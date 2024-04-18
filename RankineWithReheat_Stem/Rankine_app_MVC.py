#region imports
import sys
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from Rankine_GUI import Ui_Form
from Rankine_Classes_Reheat_MVC import rankineController
from UnitConversions import UnitConverter as UC
#these imports are necessary for drawing a matplot lib graph on my GUI
#no simple widget for this exists in QT Designer, so I have to add the widget in code.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
#endregion

#region class definitions
class MainWindow(qtw.QWidget, Ui_Form):
    def __init__(self):
        """
        MainWindow constructor
        """
        super().__init__()  #if you inherit, you generally should run the parent constructor first.
        # Main UI code goes here
        self.setupUi(self)
        self.MakeCanvas()

        #Catagorize the widgets as display or input
        self.DisplayWidgets = (self.figure, self.canvas, self.ax, self.lbl_H1, self.lbl_H2, self.lbl_H3, self.lbl_H4, self.lbl_H5, self.lbl_H6,
                               self.lbl_TurbineWork_1, self.lbl_TurbineWork_2, self.lbl_PumpWork, self.lbl_HeatIn_1,
                               self.lbl_HeatIn_2, self.lbl_ThermalEfficiency, self.lbl_SatPropHigh, self.lbl_SatPropMid,
                               self.lbl_SatPropLow, self.cmb_XAxis, self.cmb_YAxis, self.chk_logX,
                               self.chk_logY, self.lbl_TurbineInletCondition, self.lbl_TurbineInletCondition_2,
                               self.lbl_PHigh, self.lbl_PMid, self.lbl_PLow, self.lbl_H1Units, self.lbl_H2Units,
                               self.lbl_H3Units, self.lbl_H4Units, self.lbl_H5Units, self.lbl_H6Units,
                               self.lbl_TurbineWorkUnits, self.lbl_TurbineWorkUnits_2, self.lbl_PumpWorkUnits,
                               self.lbl_HeatAddedUnits, self.lbl_HeatAddedUnits_2)
        self.InputWidgets = (self.pb_Optimize, self.rdo_Quality,  self.rdo_Quality_2, self.le_PHigh, self.le_PMid, self.le_PLow,
                             self.le_TurbineInletCondition, self.le_TurbineInletCondition_2, self.le_TurbineEff,
                             self.le_TurbineEff_2, self.rb_SI, self.chk_logX, self.chk_logY, self.cmb_XAxis,
                             self.cmb_YAxis)
        #Instantiate a new rankine controller
        self.RC=rankineController(self.InputWidgets, self.DisplayWidgets)
        #Note:  all slots should be functions of the rankine controller
        self.AssignSlots()
        self.RC.updateModel()
        # a place to store coordinates from last position on graph
        self.oldXData=0.0
        self.oldYData=0.0
        # End main ui code
        self.show()

    def AssignSlots(self):
        """
        Setup signals and slots for my program.  You should notice that all the slots
        are a part of the Rankine Controller.  This is part of good MVC design pattern.
        :return:
        """
        self.btn_Calculate.clicked.connect(self.RC.updateModel)
        self.pb_Optimize.clicked.connect(self.RC.optimize)
        self.rdo_Quality.clicked.connect(self.RC.selectQualityOrTHigh_Turbine1)
        self.rdo_THigh.clicked.connect(self.RC.selectQualityOrTHigh_Turbine1)
        self.rdo_Quality_2.clicked.connect(self.RC.selectQualityOrTHigh_Turbine2)
        self.rdo_THigh_2.clicked.connect(self.RC.selectQualityOrTHigh_Turbine2)
        self.le_PHigh.editingFinished.connect(self.RC.updateSatProps)
        self.le_PMid.editingFinished.connect(self.RC.updateSatProps)
        self.le_PLow.editingFinished.connect(self.RC.updateSatProps)
        self.rb_SI.clicked.connect(self.RC.updateUnits)
        self.rb_English.clicked.connect(self.RC.updateUnits)
        self.cmb_XAxis.currentIndexChanged.connect(self.RC.updatePlot)
        self.cmb_YAxis.currentIndexChanged.connect(self.RC.updatePlot)
        self.chk_logX.toggled.connect(self.RC.updatePlot)
        self.chk_logY.toggled.connect(self.RC.updatePlot)

    #since my main window is a widget, I can customize its events by overriding the default event
    def mouseMoveEvent_Canvas(self, event):
        self.oldXData=event.xdata if event.xdata is not None else self.oldXData
        self.oldYData=event.ydata if event.ydata is not None else self.oldYData
        sUnit='kJ/(kg*K)' if self.rb_SI.isChecked() else 'BTU/(lb*R)'
        TUnit='C' if self.rb_SI.isChecked() else 'F'
        self.setWindowTitle('s:{:0.2f} {}, T:{:0.2f} {}'.format(self.oldXData,sUnit, self.oldYData,TUnit))

    def MakeCanvas(self):
        """
        Create a place to make graph of Rankine cycle
        Step 1:  create a Figure object called self.figure
        Step 2:  create a FigureCanvasQTAgg object called self.canvas
        Step 3:  create an axes object for making plot
        Step 4:  add self.canvas to self.widgetsVerticalLayout which is a Vertical box layout
        :return:
        """
        #Step 1.
        self.figure=Figure(figsize=(1,1),tight_layout=True, frameon=True)
        #Step 2.
        self.canvas=FigureCanvasQTAgg(self.figure)
        #Step 3.
        self.ax = self.figure.add_subplot()
        #Step 4.
        self.Layout_Plot.addWidget(NavigationToolbar2QT(self.canvas,self))
        self.Layout_Plot.addWidget(self.canvas)
        self.canvas.mpl_connect("motion_notify_event", self.mouseMoveEvent_Canvas)


#endregion

#if this module is being imported, this won't run. If it is the main module, it will run.
if __name__== '__main__':
    app = qtw.QApplication(sys.argv)
    mw = MainWindow()
    mw.setWindowTitle('Rankine calculator')
    sys.exit(app.exec())