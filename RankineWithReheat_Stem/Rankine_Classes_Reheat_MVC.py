#region imports
import math
from Calc_state import *
from UnitConversions import UnitConverter as UC
import numpy as np
from matplotlib import pyplot as plt
from copy import deepcopy as dc
from scipy.optimize import minimize
#these imports are necessary for drawing a matplot lib graph on my GUI
#no simple widget for this exists in QT Designer, so I have to add the widget in code.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
#endregion

#region class definitions
class rankineModel():
    def __init__(self):
        '''
        Constructor for rankine model with single reheat power cycle data (in the Model-View-Controller design pattern).
        This class is for storing data only.  The Controller class should update the model depending on input
        from the user.  The View class should display the model data depending on the desired output.
        '''
        self.p_low=None  # the low-pressure isobar
        self.p_mid=None   # JES new for project 2024, the mid-pressure isobar
        self.p_high=None  # the high-pressure isobar
        self.t_high=None  # the high temperature for the first turbine
        self.t_high_2=None   # JES new for project 2024, the high temperature for the second turbine
        self.name=None
        self.efficiency=None  # the thermal efficiency of the whole cycle NetWork/HeatIn
        self.turbine_eff=None  # the isentropic efficiency of the first turbine
        self.turbine_eff_2=None  # JES new for project 2024, the isentropic efficiency of the second turbine
        self.turbine_work=None  # the work output of the first turbine
        self.turbine_work_2=None   # JES new for project 2024, the work output of the second turbine
        self.pump_work=None  # the work input of the pump
        self.heat_added=None  # heat added by first boiler
        self.heat_added_2=None  # heat added by second turbine
        self.steam=Steam_SI()  # Instantiate a steam object for calculating state
        # Initialize the states as stateProps objects (i.e., stateProperties)
        self.state1 = stateProps()  # turbine 1 inlet
        self.state2s = stateProps()  # isentropic turbine 1 outlet (turbine 2 inlet)
        self.state2 = stateProps()  # actual turbine 1 outlet (turbine 2 inlet)
        self.state3 = stateProps()  # turbine 2 inlet
        self.state4s = stateProps()   # JES new for project 2024, isentropic turbine 2 outlet
        self.state4 = stateProps()   # JES new for project 2024, actual turbine 2 outlet
        self.state5 = stateProps()   # JES new for project 2024, pump inlet
        self.state6 = stateProps()   # JES new for project 2024, main boiler inlet
        self.SI=True  # If False, convert pressures from psi to kPa and T from F to C for inputs
        #the following are a place to store data for plotting
        self.satLiqPlotData = StateDataForPlotting()
        self.satVapPlotData = StateDataForPlotting()
        self.upperCurve = StateDataForPlotting()
        self.lowerCurve = StateDataForPlotting()
class rankineView():
    def __init__(self):
        """
        Empty constructor by design
        """
    def MakeCanvas(self):
            """
            Create a place to make graph of Rankine cycle
            Step 1:  create a Figure object called self.figure
            Step 2:  create a FigureCanvasQTAgg object called self.canvas
            Step 3:  create an axes object for making plot
            Step 4:  add self.canvas to self.widgetsVerticalLayout which is a Vertical box layout
            :return:
            """
            # Step 1.
            self.figure = Figure(figsize=(1, 1), tight_layout=True, frameon=True)
            # Step 2.
            self.canvas = FigureCanvasQTAgg(self.figure)
            # Step 3.
            self.ax = self.figure.add_subplot()
            # Step 4.
            self.Layout_Plot.addWidget(NavigationToolbar2QT(self.canvas))
            self.Layout_Plot.addWidget(self.canvas)
    def SetupCanvasInteraction(self, window):
        self.canvas.mpl_connect("motion_notify_event", window.mouseMoveEvent_Canvas)
    def setWidgets(self, *args):
        self.InputWidgets, self.DisplayWidgets = args
        self.UnpackWidgets()
        self.MakeCanvas()
    def UnpackWidgets(self):
        #Unpacks display widgets
        self.Layout_Plot, self.lbl_H1, self.lbl_H2, self.lbl_H3, self.lbl_H4, self.lbl_H5, self.lbl_H6,\
        self.lbl_TurbineWork_1, self.lbl_TurbineWork_2, self.lbl_PumpWork, self.lbl_HeatIn_1,\
        self.lbl_HeatIn_2, self.lbl_ThermalEfficiency, self.lbl_SatPropHigh, self.lbl_SatPropMid,\
        self.lbl_SatPropLow, self.cmb_XAxis, self.cmb_YAxis, self.chk_logX,\
        self.chk_logY, self.lbl_TurbineInletCondition, self.lbl_TurbineInletCondition_2,\
        self.lbl_PHigh, self.lbl_PMid, self.lbl_PLow, self.lbl_H1Units, self.lbl_H2Units,\
        self.lbl_H3Units, self.lbl_H4Units, self.lbl_H5Units, self.lbl_H6Units,\
        self.lbl_TurbineWorkUnits, self.lbl_TurbineWorkUnits_2, self.lbl_PumpWorkUnits,\
        self.lbl_HeatAddedUnits, self.lbl_HeatAddedUnits_2=self.DisplayWidgets
        #Unpacks input widgets
        self.pb_Optimize, self.rdo_Quality,  self.rdo_Quality_2, self.le_PHigh, self.le_PMid, self.le_PLow,\
        self.le_TurbineInletCondition, self.le_TurbineInletCondition_2, self.le_TurbineEff,\
        self.le_TurbineEff_2, self.rb_SI, self.chk_logX, self.chk_logY, self.cmb_XAxis,\
        self.cmb_YAxis = self.InputWidgets
        self.pb_Optimize.setToolTip("WARNING:  This function looks for an optimal P Mid to maximize Thermal Efficiency"
                                    "\nof the Rankine Cycle.  It does not change the Turbine inlet temperature(s) or "
                                    "\nquality, only the P Mid value.")
    def updateSatProps(self, Model=None):
        """
         This function checks the self.rb_SI.isChecked() to see which units are used.
         Then, it sets the text of lbl_SatPropHigh using SatPropsIsobar(float(self.le_PHigh.text())*PCF, SI=SI).txtOut
         here, PCF is the pressure conversion factor
         finally, we need to call the function self.SelectQualityOrTHigh()
         :return:
         """
        SI=self.rb_SI.isChecked()
        PCF=1 if SI else UC.psi_to_bar
        #update PHigh saturated properties
        satProp=Model.steam.getsatProps_p(float(self.le_PHigh.text())*PCF)
        self.lbl_SatPropHigh.setText(satProp.getTextOutput(SI=SI))
        #update PMid saturated properties
        satProp=Model.steam.getsatProps_p(float(self.le_PMid.text())*PCF)
        self.lbl_SatPropMid.setText(satProp.getTextOutput(SI=SI))
        #update PLow saturated prooperties
        satProp=Model.steam.getsatProps_p(float(self.le_PLow.text())*PCF)
        self.lbl_SatPropLow.setText(satProp.getTextOutput(SI=SI))
    def selectQualityOrTHigh_Turbine_1(self, Model=None):
        """
        Action to take when selecting one of the radio buttons for Quality or THigh
        :return:
        """
        SI = self.rb_SI.isChecked()
        if self.rdo_Quality.isChecked():
            self.le_TurbineInletCondition.setText("1.0")
            self.le_TurbineInletCondition.setEnabled(False)
        else:
            PCF = 1 if SI else UC.psi_to_bar
            satPHigh = Model.steam.getsatProps_p(float(self.le_PHigh.text()) * PCF)
            Tsat = satPHigh.tsat
            Tsat = Tsat if SI else UC.C_to_F(Tsat)
            CurrentT = float(self.le_TurbineInletCondition.text())
            CurrentT = max(CurrentT, Tsat)
            self.le_TurbineInletCondition.setText("{:0.2f}".format(CurrentT))
            self.le_TurbineInletCondition.setEnabled(True)
        x = self.rdo_Quality.isChecked()
        self.lbl_TurbineInletCondition.setText(
            ("Turbine Inlet: {}{} =".format('x' if x else 'THigh', '' if x else ('(C)' if SI else '(F)'))))
    def selectQualityOrTHigh_Turbine_2(self, Model=None):
        """
        Action to take when selecting one of the radio buttons for Quality or THigh
        :return:
        """
        SI = self.rb_SI.isChecked()
        if self.rdo_Quality_2.isChecked():
            self.le_TurbineInletCondition_2.setText("1.0")
            self.le_TurbineInletCondition_2.setEnabled(False)
        else:
            PCF = 1 if SI else UC.psi_to_bar
            satPMid = Model.steam.getsatProps_p(float(self.le_PMid.text()) * PCF)
            Tsat = satPMid.tsat
            Tsat = Tsat if SI else UC.C_to_F(Tsat)
            CurrentT = float(self.le_TurbineInletCondition_2.text())
            CurrentT = max(CurrentT, Tsat)
            self.le_TurbineInletCondition_2.setText("{:0.2f}".format(CurrentT))
            self.le_TurbineInletCondition_2.setEnabled(True)
        x = self.rdo_Quality_2.isChecked()
        self.lbl_TurbineInletCondition_2.setText(
            ("Turbine Inlet: {}{} =".format('x' if x else 'THigh', '' if x else ('(C)' if SI else '(F)'))))
    def updateGUI(self, Model=None):
        """
        The model should already be updated.  This only updates the GUI.
        Steps:

        :param args:
        :param Model:
        :return:
        """
        #unpack the args
        if Model.state1 is None:  # means the cycle has not been evaluated yet
            return
        steam=Model.steam
        #update the line edits and labels
        HCF=1 if Model.SI else UC.kJperkg_to_BTUperlb
        self.lbl_H1.setText("{:0.2f}".format(Model.state1.h * HCF))
        self.lbl_H2.setText("{:0.2f}".format(Model.state2.h * HCF))
        self.lbl_H3.setText("{:0.2f}".format(Model.state3.h * HCF))
        self.lbl_H4.setText("{:0.2f}".format(Model.state4.h * HCF))
        self.lbl_H5.setText("{:0.2f}".format(Model.state5.h * HCF))
        self.lbl_H6.setText("{:0.2f}".format(Model.state6.h * HCF))

        self.lbl_TurbineWork_1.setText("{:0.2f}".format(Model.turbine_work*HCF))
        self.lbl_TurbineWork_2.setText("{:0.2f}".format(Model.turbine_work_2*HCF))

        self.lbl_PumpWork.setText("{:0.2f}".format(Model.pump_work*HCF))
        self.lbl_HeatIn_1.setText("{:0.2f}".format(Model.heat_added*HCF))
        self.lbl_HeatIn_2.setText("{:0.2f}".format(Model.heat_added_2*HCF))
        self.lbl_ThermalEfficiency.setText("{:0.2f}".format(Model.efficiency))
        satPropsLow=steam.getsatProps_p(p=Model.p_low)
        satPropsMid=steam.getsatProps_p(p=Model.p_mid)
        satPropsHigh=steam.getsatProps_p(p=Model.p_high)
        self.lbl_SatPropLow.setText(satPropsLow.getTextOutput(SI=Model.SI))
        self.lbl_SatPropMid.setText(satPropsMid.getTextOutput(SI=Model.SI))
        self.lbl_SatPropHigh.setText(satPropsHigh.getTextOutput(SI=Model.SI))

        #update the plot
        self.plot_cycle_XY(Model=Model)
    def updateUnits(self, Model=None):
        """
        Updates the units on the GUI to match choice of SI or English
        :param W: A tuple of line edit widgets for input and output as well as the graph
        :param OW: A list of label widgets for units
        :param Model:  a reference to the model
        :return:
        """
        #Step 0. update the outputs
        self.updateGUI(Model=Model)
        # Update units displayed on labels
        #Step 1. Update pressures for PHigh and PLow
        pCF=1 if Model.SI else UC.bar_to_psi
        self.le_PHigh.setText("{:0.2f}".format(pCF * Model.p_high))
        self.le_PMid.setText("{:0.2f}".format(pCF * Model.p_mid))
        self.le_PLow.setText("{:0.2f}".format(pCF * Model.p_low))
        #Step 2. Update THigh for Turbine 1 if it is not None
        if not self.rdo_Quality.isChecked():
            T=float(self.le_TurbineInletCondition.text())
            T=float(self.le_TurbineInletCondition.text())
            TUnits = "C" if Model.SI else "F"
            self.le_TurbineInletCondition.setText("{:0.2f}".format(T))
            self.lbl_TurbineInletCondition.setText("Turbine 1 Inlet: THigh ({}):".format(TUnits))
        #Update THigh for Turbine 2
        if not self.rdo_Quality_2.isChecked():
            T3=float(self.le_TurbineInletCondition_2.text())
            T3=float(self.le_TurbineInletCondition_2.text())
            TUnits = "C" if Model.SI else "F"
            self.le_TurbineInletCondition_2.setText("{:0.2f}".format(T3))
            self.lbl_TurbineInletCondition_2.setText("Turbine 2 Inlet: THigh ({}):".format(TUnits))
        #Step 3. Update the units for labels
        self.lbl_PHigh.setText("P High ({})".format('bar' if Model.SI else 'psi'))
        self.lbl_PMid.setText("P Mid ({})".format('bar' if Model.SI else 'psi'))
        self.lbl_PLow.setText("P Low ({})".format('bar' if Model.SI else 'psi'))
        HUnits="kJ/kg" if Model.SI else "BTU/lb"
        self.lbl_H1Units.setText(HUnits)
        self.lbl_H2Units.setText(HUnits)
        self.lbl_H3Units.setText(HUnits)
        self.lbl_H4Units.setText(HUnits)
        self.lbl_H5Units.setText(HUnits)
        self.lbl_H6Units.setText(HUnits)
        self.lbl_TurbineWorkUnits.setText(HUnits)
        self.lbl_TurbineWorkUnits_2.setText(HUnits)
        self.lbl_PumpWorkUnits.setText(HUnits)
        self.lbl_HeatAddedUnits.setText(HUnits)
    def print_summary(self, Model=None):
        """
        Prints to CLI.
        :param Model: a rankineModel object
        :return: nothing
        """
        if Model.efficiency==None:
            Model.calc_efficiency()
        print('Cycle Summary for: ', Model.name)
        print('\tEfficiency: {:0.3f}%'.format(Model.efficiency))
        print('\tTurbine Eff:  {:0.2f}'.format(Model.turbine_eff))
        print('\tTurbine Work: {:0.3f} kJ/kg'.format(Model.turbine_work))
        print('\tPump Work: {:0.3f} kJ/kg'.format(Model.pump_work))
        print('\tHeat Added: {:0.3f} kJ/kg'.format(Model.heat_added))
        Model.state1.print()
        Model.state2.print()
        Model.state3.print()
        Model.state4.print()
    def plot_cycle_TS(self, axObj=None, Model=None):
        """
        I want to plot the Rankine cycle on T-S coordinates along with the vapor dome and shading in the cycle.
        I notice there are several lines on the plot:
        saturated liquid T(s) colored blue
        saturated vapor T(s) colored red
        The high and low isobars and lines connecting state 1 to 2, and 3 to saturated liquid at phigh
        step 1:  build data for saturated liquid line
        step 2:  build data for saturated vapor line
        step 3:  build data between state 3 and sat liquid at p_high
        step 4:  build data between sat liquid at p_high and state 1
        step 5:  build data between state 1 and state 2
        step 6:  build data between state 2 and state 3
        step 7:  put together data from 3,4,5 for top line and build bottom line
        step 8:  make and decorate plot

        Note:  will plot using pyplot if axObj is None else just returns

        :param axObj:  if None, used plt.subplot else a MatplotLib axes object.
        :return:
        """
        SI=Model.SI
        steam=Model.steam
        #region step 1&2:
        ts, ps, hfs, hgs, sfs, sgs, vfs, vgs = np.loadtxt('sat_water_table.txt', skiprows=1,
                                                          unpack=True)  # use np.loadtxt to read the saturated properties
        ax = plt.subplot() if axObj is None else axObj

        hCF = 1 if Model.SI else UC.kJperkg_to_BTUperlb
        pCF = 1 if Model.SI else UC.kpa_to_psi
        sCF = 1 if Model.SI else UC.kJperkgK_to_BTUperlbR
        vCF = 1 if Model.SI else UC.kgperm3_to_lbperft3

        sfs *= sCF
        sgs *= sCF
        hfs *= hCF
        hgs *= hCF
        vfs *= vCF
        vgs *= vCF
        ps *= pCF
        ts = [t if Model.SI else UC.C_to_F(t) for t in ts]

        xfsat = sfs
        yfsat = ts
        xgsat = sgs
        ygsat = ts

        ax.plot(xfsat, yfsat, color='blue')
        ax.plot(xgsat, ygsat, color='red')
        #endregion

        #step 3:  I'll just make a straight line between state3 and state3p
        st3p=steam.getState(Model.p_high,x=0) #saturated liquid state at p_high
        svals=np.linspace(Model.state3.s, st3p.s, 20)
        hvals=np.linspace(Model.state3.h, st3p.h, 20)
        pvals=np.linspace(Model.p_low, Model.p_high,20)
        vvals=np.linspace(Model.state3.v, st3p.v, 20)
        tvals=np.linspace(Model.state3.T, st3p.T, 20)
        line3=np.column_stack([svals, tvals])

        #step 4:
        sat_pHigh=steam.getState(Model.p_high, x=1.0)
        st1=Model.state1
        svals2p=np.linspace(st3p.s, sat_pHigh.s, 20)
        hvals2p = np.linspace(st3p.h, sat_pHigh.h, 20)
        pvals2p = [Model.p_high for i in range(20)]
        vvals2p = np.linspace(st3p.v, sat_pHigh.v, 20)
        tvals2p=[st3p.T for i in range(20)]
        line4=np.column_stack([svals2p, tvals2p])
        if st1.T>sat_pHigh.T:  #need to add data points to state1 for superheated
            svals_sh=np.linspace(sat_pHigh.s,st1.s, 20)
            tvals_sh=np.array([steam.getState(Model.p_high,s=ss).T for ss in svals_sh])
            line4 =np.append(line4, np.column_stack([svals_sh, tvals_sh]), axis=0)
        #plt.plot(line4[:,0], line4[:,1])

        #step 5:
        svals=np.linspace(Model.state1.s, Model.state2.s, 20)
        tvals=np.linspace(Model.state1.T, Model.state2.T, 20)
        line5=np.array(svals)
        line5=np.column_stack([line5, tvals])
        #plt.plot(line5[:,0], line5[:,1])

        #step 6:
        svals=np.linspace(Model.state2.s, Model.state3.s, 20)
        tvals=np.array([Model.state2.T for i in range(20)])
        line6=np.column_stack([svals, tvals])
        #plt.plot(line6[:,0], line6[:,1])

        #step 7:
        topLine=np.append(line3, line4, axis=0)
        topLine=np.append(topLine, line5, axis=0)
        xvals=topLine[:,0]
        y1=topLine[:,1]
        y2=[Model.state3.T for s in xvals]

        if not SI:
            xvals*=UC.kJperkgK_to_BTUperlbR
            for i in range(len(y1)):
                y1[i]=UC.C_to_F(y1[i])
            for i in range(len(y2)):
                y2[i]=UC.C_to_F(y2[i])

        ax.plot(xvals, y1, color='darkgreen')
        ax.plot(xvals, y2, color='black')
        # ax.fill_between(xvals, y1, y2, color='gray', alpha=0.5)

        if SI:
            ax.plot(Model.state1.s, Model.state1.T, marker='o', markeredgecolor='k', markerfacecolor='w')
            ax.plot(Model.state2.s, Model.state2.T, marker='o', markeredgecolor='k', markerfacecolor='w')
            ax.plot(Model.state3.s, Model.state3.T, marker='o', markeredgecolor='k', markerfacecolor='w')
        else:
            ax.plot(Model.state1.s * UC.kJperkgK_to_BTUperlbR, UC.C_to_F(Model.state1.T), marker='o', markeredgecolor='k', markerfacecolor='w')
            ax.plot(Model.state2.s * UC.kJperkgK_to_BTUperlbR, UC.C_to_F(Model.state2.T), marker='o', markeredgecolor='k', markerfacecolor='w')
            ax.plot(Model.state3.s * UC.kJperkgK_to_BTUperlbR, UC.C_to_F(Model.state3.T), marker='o', markeredgecolor='k', markerfacecolor='w')

        tempUnits=r'$\left(^oC\right)$' if SI else r'$\left(^oF\right)$'
        entropyUnits=r'$\left(\frac{kJ}{kg\cdot K}\right)$' if SI else r'$\left(\frac{BTU}{lb\cdot ^oR}\right)$'
        ax.set_xlabel(r's '+entropyUnits, fontsize=18)  #different than plt
        ax.set_ylabel(r'T '+tempUnits, fontsize=18)  #different than plt
        ax.set_title(Model.name)  #different than plt
        ax.grid(visible='both', alpha=0.5)
        ax.tick_params(axis='both', direction='in', labelsize=18)

        sMin=min(sfs)
        sMax=max(sgs)
        ax.set_xlim(sMin, sMax)  #different than plt

        tMin=min(ts)
        tMax=max(max(ts),st1.T)
        ax.set_ylim(tMin,tMax*1.05)  #different than plt

        energyUnits=r'$\frac{kJ}{kg}$' if SI else r'$\frac{BTU}{lb}$'
        energyCF = 1 if SI else UC.kJperkg_to_BTUperlb

        if axObj is None:  # this allows me to show plot if not being displayed on a figure
            plt.show()
    def plot_cycle_XY(self, Model=None):
        """
        I want to plot any two thermodynaimc properties on X and Y
        :param X: letter for which variable to plot on X axis
        :param Y: letter for which variable to plot on Y axis
        :return:
        """
        X = self.cmb_XAxis.currentText()
        Y = self.cmb_YAxis.currentText()
        logx = self.chk_logX.isChecked()
        logy = self.chk_logY.isChecked()
        SI=Model.SI
        if X == Y:
            return
        QTPlotting = True  # assumes we are plotting onto a QT GUI form
        if self.ax == None:
            self.ax = plt.subplot()
            QTPlotting = False  # actually, we are just using CLI and showing the plot

        self.ax.clear()
        self.ax.set_xscale('log' if logx else 'linear')
        self.ax.set_yscale('log' if logy else 'linear')
        YF = Model.satLiqPlotData.getDataCol(Y, SI=SI)
        YG = Model.satVapPlotData.getDataCol(Y, SI=SI)
        XF = Model.satLiqPlotData.getDataCol(X, SI=SI)
        XG = Model.satVapPlotData.getDataCol(X, SI=SI)
        # plot the vapor dome
        self.ax.plot(XF, YF, color='b')
        self.ax.plot(XG, YG, color='r')
        # plot the upper and lower curves
        self.ax.plot(Model.lowerCurve.getDataCol(X, SI=SI), Model.lowerCurve.getDataCol(Y, SI=SI), color='k')
        self.ax.plot(Model.upperCurve.getDataCol(X, SI=SI), Model.upperCurve.getDataCol(Y, SI=SI), color='g')
        # ax.fill_between(Model.upperCurve.getDataCol(X), Model.upperCurve.getDataCol(Y), self.lowerCurve.getDataCol(Y), color='grey', alpha=0.2)

        # add axis labels
        self.ax.set_ylabel(Model.lowerCurve.getAxisLabel(Y, SI=SI), fontsize='large' if QTPlotting else 'medium')
        self.ax.set_xlabel(Model.lowerCurve.getAxisLabel(X, SI=SI), fontsize='large' if QTPlotting else 'medium')
        # put a title on the plot
        Model.name = 'Rankine Cycle - ' + Model.state1.region + ' at Turbine Inlet'
        self.ax.set_title(Model.name, fontsize='large' if QTPlotting else 'medium')

        # modify the tick marks
        self.ax.tick_params(axis='both', which='both', direction='in', top=True, right=True,
                       labelsize='large' if QTPlotting else 'medium')  # format tick marks

        # plot the circles for states 1, 2, 3, and 4
        self.ax.plot(Model.state1.getVal(X, SI=SI), Model.state1.getVal(Y, SI=SI), marker='o', markerfacecolor='w',
                markeredgecolor='k')
        self.ax.plot(Model.state2.getVal(X, SI=SI), Model.state2.getVal(Y, SI=SI), marker='o', markerfacecolor='w',
                markeredgecolor='k')
        self.ax.plot(Model.state3.getVal(X, SI=SI), Model.state3.getVal(Y, SI=SI), marker='o', markerfacecolor='w',
                markeredgecolor='k')
        self.ax.plot(Model.state4.getVal(X, SI=SI), Model.state4.getVal(Y, SI=SI), marker='o', markerfacecolor='w',
                markeredgecolor='k')
        self.ax.plot(Model.state5.getVal(X, SI=SI), Model.state5.getVal(Y, SI=SI), marker='o', markerfacecolor='w',
                markeredgecolor='k')
        self.ax.plot(Model.state6.getVal(X, SI=SI), Model.state6.getVal(Y, SI=SI), marker='o', markerfacecolor='w',
                markeredgecolor='k')
        # set limits on x and y
        xmin = min(min(XF), min(XG), min(Model.upperCurve.getDataCol(X, SI=SI)), max(Model.lowerCurve.getDataCol(X, SI=SI)))
        xmax = max(max(XF), max(XG), max(Model.upperCurve.getDataCol(X, SI=SI)), max(Model.lowerCurve.getDataCol(X, SI=SI)))
        ymin = min(min(YF), min(YG), min(Model.upperCurve.getDataCol(Y, SI=SI)), max(Model.lowerCurve.getDataCol(Y, SI=SI)))
        ymax = max(max(YF), max(YG), max(Model.upperCurve.getDataCol(Y, SI=SI)),
                   max(Model.lowerCurve.getDataCol(Y, SI=SI))) * 1.1
        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)
        deltax = xmax - xmin
        deltay = ymax - ymin
        # add the summary text to the plot

        # show the plot
        if QTPlotting == False:
            plt.show()
        else:
            self.canvas.draw()
class rankineController():
    def __init__(self, *args):
        """
        Create rankineModel object.  The rankineController class updates the model based on user input
        and updates the rankineView as well
        :param *args: a tuple containing (inputWidgets, displayWidgets)
        """
        self.InputWidgets, self.DisplayWidgets = args
        self.UnpackInputWidgets()
        self.Model = rankineModel()
        self.View = rankineView()
        self.View.setWidgets(self.InputWidgets, self.DisplayWidgets)
        self.buildVaporDomeData()
        self.iterations=0
    def UnpackInputWidgets(self):
        #Unpacks input widgets
        self.pb_Optimize, self.rdo_Quality,  self.rdo_Quality_2, self.le_PHigh, self.le_PMid, \
        self.le_PLow, self.le_TurbineInletCondition, self.le_TurbineInletCondition_2, \
        self.le_TurbineEff, self.le_TurbineEff_2, self.rb_SI, self.chk_logX, \
        self.chk_logY, self.cmb_XAxis, self.cmb_YAxis = self.InputWidgets

    def SetupCanvasInteraction(self, window):
        self.View.SetupCanvasInteraction(window)
    def optimize(self):
        def objFn(P):
            """
            This is the objective function to be minimized by the minimize function.
            It computes the cycle thermal efficiency.
            It applies a steep penalty if the p_mid exceeds p_high or is below p_low.
            """
            p=P[0]
            self.Model.p_mid = p
            eff = self.calc_efficiency()
            #region penalty functions
            if p >= self.Model.p_high:
                pass #$JES MISSING CODE# # implement a penalty
            if p <= self.Model.p_low:
                pass  # $JES MISSING CODE# # implement a penalty
            if self.Model.heat_added_2 <0:
                pass  # $JES MISSING CODE# # implement a penalty
            #endregion
            return 100.0 - eff
        #set initial guess for p_mid as a numpy array
        ic=np.array([self.Model.p_low + 0.01 * (self.Model.p_high - self.Model.p_low)])
        #set self.Model.t_high and self.Model.t_high_2 to match view
        self.readConditionsFromGUI()
        #optimize efficiency by adjusting p_mid
        I = minimize(objFn, ic, method="Nelder-Mead")
        #Apply optimized p_mid
        self.Model.p_mid = I.x[0]
        self.le_PMid.setText("{:0.3f}".format(self.Model.p_mid))
        self.updateModel()
        pass
    def readConditionsFromGUI(self):
        """
        This reads the line edit boxes and radio buttons from the GUI
        """
        T = float(self.le_TurbineInletCondition.text())  # $UNITS$
        T_2 = float(self.le_TurbineInletCondition_2.text())
        x1 = self.rdo_Quality.isChecked()
        x2 = self.rdo_Quality_2.isChecked()

        self.Model.t_high = None if x1 else (T if self.Model.SI else UC.F_to_C(T))  # $UNITS$
        self.Model.t_high_2 = None if x2 else (T_2 if self.Model.SI else UC.F_to_C(T_2))  # $UNITS$
        self.Model.turbine_eff = float(self.le_TurbineEff.text())
        self.Model.turbine_eff_2 = float(self.le_TurbineEff_2.text())

        p_CF = 1.0 if self.Model.SI else UC.psi_to_bar
        self.Model.p_high = p_CF * float(self.le_PHigh.text())
        self.Model.p_mid = p_CF * float(self.le_PMid.text())
        self.Model.p_low = p_CF * float(self.le_PLow.text())

        self.Model.turbine_eff = float(self.le_TurbineEff.text())
        self.Model.turbine_eff_2 = float(self.le_TurbineEff_2.text())
    def updateModel(self):
        # update the model
        self.readConditionsFromGUI()
        # do the calculation
        self.calc_efficiency()
        self.updateView()
    def updateUnits(self):
        # Switching units should not change the model, but should update the view
        self.Model.SI = self.rb_SI.isChecked()
        self.View.updateUnits(Model=self.Model)
    def calc_efficiency(self):
        """
        I've modified this on 4/15/2022 to use a single SI_Steam object that is held in the model for calculating
        various states along the path of the Rankine cycle.  I use the getState function to retrieve a deep copy of
        a stateProps object.
        :return:
        """
        steam = self.Model.steam
        # calculate the 6 states

        # state 1: turbine inlet (p_high, t_high) superheated or saturated vapor
        if (self.Model.t_high == None):
            x=float(self.le_TurbineInletCondition.text())
            self.Model.state1 = steam.getState(P=self.Model.p_high, x=x, name='Turbine 1 Inlet')
        else:
            self.Model.state1 = steam.getState(P=self.Model.p_high, T=self.Model.t_high, name='Turbine 1 Inlet')

        # state 2: turbine exit (p_low, s=s_turbine inlet) two-phase
        self.Model.state2s = steam.getState(P=self.Model.p_mid, s=self.Model.state1.s, name="Turbine 1 Exit")
        if self.Model.turbine_eff < 1.0:  # eff=(h1-h2)/(h1-h2s) -> h2=h1-eff(h1-h2s)
            h2 = self.Model.state1.h - self.Model.turbine_eff * (self.Model.state1.h - self.Model.state2s.h)
            self.Model.state2 = steam.getState(P=self.Model.p_mid, h=h2, name="Turbine 1 Exit")
        else:
            self.Model.state2 = self.Model.state2s

        # state 3: turbine 2 inlet (p_mid, T=self.Model.T_high_2) 2-phase or superheated
        if (self.Model.t_high_2 == None):
            x = float(self.le_TurbineInletCondition_2.text())
            self.Model.state3 = steam.getState(P=self.Model.p_mid, x=x, name='Turbine 2 Inlet')
        else:
            self.Model.state3 = steam.getState(P=self.Model.p_mid, T=self.Model.t_high_2, name='Turbine 2 Inlet')

        # state 4: turbine 2 exit (p_low,s=s_turbine_2_inlet)
        self.Model.state4s = steam.getState(P=self.Model.p_low, s=self.Model.state3.s, name='Turbine 2 Exit')
        if self.Model.turbine_eff_2 < 1.0:  # eff=(h3-h4)/(h3-h4s) -> h4=h3-eff(h3-h42s)
            h4 = self.Model.state3.h - self.Model.turbine_eff_2 * (self.Model.state3.h - self.Model.state4s.h)
            self.Model.state4 = steam.getState(P=self.Model.p_low, h=h4, name="Turbine 2 Exit")
        else:
            self.Model.state4 = self.Model.state4s

        # state 5: pump inlet
        self.Model.state5 = steam.getState(P=self.Model.p_low, x=0.0, name="Pump Inlet")

        # state 6: boiler inlet
        self.Model.state6 = steam.getState(P=self.Model.p_high, s=self.Model.state5.s, name="Boiler Inlet")

        self.Model.turbine_work = self.Model.state1.h - self.Model.state2.h  # calculate turbine work
        self.Model.turbine_work_2 = self.Model.state3.h - self.Model.state4.h  # calculate turbine 2 work
        self.Model.pump_work = self.Model.state6.h - self.Model.state5.h  # calculate pump work
        self.Model.heat_added = self.Model.state1.h - self.Model.state6.h  # calculate heat added
        self.Model.heat_added_2 = self.Model.state3.h - self.Model.state2.h  # calculate heat added 2
        self.Model.efficiency = 100.0 * (self.Model.turbine_work + self.Model.turbine_work_2 - self.Model.pump_work) / (
                self.Model.heat_added + self.Model.heat_added_2)
        return self.Model.efficiency
    def updateSatProps(self):
        self.View.updateSatProps(self.Model)
    def selectQualityOrTHigh_Turbine1(self):
        self.View.selectQualityOrTHigh_Turbine_1(self.Model)
    def selectQualityOrTHigh_Turbine2(self):
        self.View.selectQualityOrTHigh_Turbine_2(self.Model)
    def updateView(self):
        """
        This is a pass-through function that calls and identically named function in the View, but passes along the
        Model as an argument.
        :param args: A tuple of Widgets that get unpacked and updated in the view
        :return:
        """
        self.buildDataForPlotting()
        self.View.updateGUI(Model=self.Model)
    def print_summary(self):
        """
        A pass-thrugh method for accessing View and passing Model.
        :return:
        """
        self.View.print_summary(Model=self.Model)
    def buildVaporDomeData(self, nPoints=500):
        """
        I'll build the vapor dome from just above the triple point up to the critical point
        :param nPoints:
        :return:
        """
        steam = self.Model.steam
        tp = triplePt_PT()
        cp = criticalPt_PT()
        steam.state.p = cp.p  #in bar
        steam.state.t = cp.t  #in K
        steam.calcState_1Phase()
        critProps = dc(steam.state)
        critProps.x=1.0
        P = np.logspace(math.log10(tp.p * 1.001), math.log10(cp.p * 0.99), nPoints)
        for p in P:
            sat = steam.getsatProps_p(p)
            self.Model.satLiqPlotData.addPt((sat.tsat, p, sat.uf, sat.hf, sat.sf, sat.vf))
            self.Model.satVapPlotData.addPt((sat.tsat, p, sat.uf, sat.hg, sat.sg, sat.vg))
        self.Model.satLiqPlotData.addPt((critProps.t, critProps.p, critProps.u, critProps.h, critProps.s, critProps.v))
        self.Model.satVapPlotData.addPt((critProps.t, critProps.p, critProps.u, critProps.h, critProps.s, critProps.v))
    def buildDataForPlotting(self):
        """
        I want to create data for plotting the Rankine cycle.  The key states are:
        State 1.  Entrance to Turbine (either saturated vapor or superheated steam at p_High)
        State 2.  Entrance to reheat
        State 3.  Entrance to Turbine 2
        State 4.  Entrance to Condenser (probably two-phase at p_Low)
        State 5.  Entrance to the pump (saturated liquid at p_Low)
        State 6.  Entrance to the boiler (sub-cooled liquid at p_High, isentropic pump)

        I want to create h, s, v, p, T data between states 1-2, 2-3, 3-4, 4-5, 5-6, 6-1
        I'll piece together an upperCurve data set from 5-6 + 6-1 + 1-2 + 2-3 + 3-4
        The lowerCurve data set is 4-5
        :return:
        """
        # clear out any old data
        self.Model.upperCurve.clear()
        self.Model.lowerCurve.clear()

        # get saturated properties at PHigh, PMid,  and PLow
        satPLow = self.Model.steam.getsatProps_p(self.Model.p_low)
        satPMid = self.Model.steam.getsatProps_p(self.Model.p_mid)
        satPHigh = self.Model.steam.getsatProps_p(self.Model.p_high)

        steam = self.Model.steam

        # region build upperCurve
        # region states from 5-6
        nPts = 15
        DeltaP = (satPHigh.psat - satPLow.psat)
        state = steam.getState(P=satPLow.psat, s=satPLow.sf)
        for n in range(nPts):
            z = n * 1.0 / (nPts - 1)
            state = steam.getState(P=(satPLow.psat + z * DeltaP), s=satPLow.sf)
            self.Model.upperCurve.addStatePt(state)
        # endregion

        # region states from 6-1
        # first from T6 to T7 where T7 is the saturated liquid at p_High
        T6 = state.t
        T7 = satPHigh.tsat
        DeltaT = (T7 - T6)
        nPts = 20
        P = satPHigh.psat
        for n in range(nPts - 1):
            z = n * 1.0 / (nPts - 2)
            T = T6 + z * DeltaT
            if T < T7:
                state = steam.getState(P=P, T=T)
                self.Model.upperCurve.addStatePt(state)
        for n in range(nPts):
            z = n * 1.0 / (nPts - 1)
            state = steam.getState(satPHigh.psat, x=z)
            self.Model.upperCurve.addStatePt(state)
        if self.Model.state1.t > (satPHigh.tsat + 1):
            DeltaT = self.Model.state1.t - T7
            for n in range(0, nPts):
                z = n * 1.0 / (nPts - 1)
                if z > 0:
                    state = steam.getState(satPHigh.psat, T=T7 + z * DeltaT)
                    self.Model.upperCurve.addStatePt(state)
        # endregion

        # region states between 1 and 2
        # I'm assuming a linear change in Pressure from P1 to P2, along with linear change in s,
        # but not sure of details inside the turbine, so this is just a guess.
        s1 = self.Model.state1.s
        s2 = self.Model.state2.s
        P1 = self.Model.state1.p
        P2 = self.Model.state2.p
        Deltas = s2 - s1
        DeltaP = P2 - P1
        for n in range(nPts):
            z = n * 1.0 / (nPts - 1)
            state = steam.getState(P=P1 + z * DeltaP, s=s1 + z * Deltas)
            self.Model.upperCurve.addStatePt(state)
        # endregion

        # region states between 2 and 3
        if self.Model.state3.t >= satPMid.tsat:  # state 3 is superheated or 2-phase
            # region case i) states from 2-2Sat
            if self.Model.state2.s <= satPMid.sg:  # if needed, add points up to saturated vapor
                nPts = 5
                x2 = self.Model.state2.x
                Deltax = round(self.Model.state3.x - x2,3)
                if Deltax > 0.0:
                    for n in range(nPts):
                        z = n * 1.0 / (nPts - 1)
                        x=UC.clamp(x2 + z * Deltax,0.0, 1.0)
                        state = steam.getState(P=satPMid.psat, x=x)
                        self.Model.upperCurve.addStatePt(state)
                #state = steam.getState(P=satPMid.psat, x=1.0)
                #self.Model.upperCurve.addStatePt(state)
                T2 = satPMid.tsat
            else:
                T2 = self.Model.state2.t
            # endregion
            T3 = self.Model.state3.t
            if (T3 != T2):
                DeltaT = (T3 - T2)
                nPts = 20
                P = satPMid.psat
                for n in range(nPts):
                    z = n * 1.0 / (nPts - 1)
                    T = T2 + z * DeltaT
                    state = steam.getState(P=P, T=T)
                    self.Model.upperCurve.addStatePt(state)
            pass
        # endregion

        # region states between 3 and 4
        state = self.Model.state3
        s3 = self.Model.state3.s
        P3 = self.Model.state3.p
        s4 = self.Model.state4.s
        P4 = self.Model.state4.p
        nPts = 15
        DeltaP = P4 - P3
        Deltas = s4 - s3
        for n in range(nPts):
            z = n * 1.0 / (nPts - 1)
            state = steam.getState(P=P3 + z * DeltaP, s=s3 + z * Deltas)
            self.Model.upperCurve.addStatePt(state)
        # endregion
        # endregion
        # region build lowerCurve
        # region states between 4 and 5
        if self.Model.state4.t > satPLow.tsat:
            t4p = satPLow.tsat
            t4 = self.Model.state4.t
            deltaT = t4p - t4
            nPts = 15
            for n in range(nPts):
                z = n * 1.0 / (nPts - 1)
                state = steam.getState(P=satPLow.psat, T=t4 + z * deltaT)
                self.Model.lowerCurve.addStatePt(state)
        x5 = 0
        x4 = self.Model.state4.x
        deltax = x5 - x4
        nPts = 15
        for n in range(nPts):
            z = n * 1.0 / (nPts - 1)
            state = steam.getState(P=satPLow.psat, x=UC.clamp(x4 + z * deltax,0.0,1.0))
            self.Model.lowerCurve.addStatePt(state)
        # endregion
        # endregion
        pass
    def updatePlot(self):
        self.View.plot_cycle_XY(Model=self.Model)
#endregion