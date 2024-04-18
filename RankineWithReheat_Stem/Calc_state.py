#region imports
import numpy as np
from copy import deepcopy as dc
from scipy.optimize import fsolve
from UnitConversions import UnitConverter as UC
from pyXSteam.XSteam import XSteam
from pyXSteam import Constants
#endregion

#region class definitions
class triplePt_PT():
    def __init__(self):
        """
        get the triple point pressure in bar and temperature in C
        """
        self.t=Constants.__TRIPLE_POINT_TEMPERATURE__+Constants.__ABSOLUTE_ZERO_CELSIUS__
        self.p=Constants.__TRIPLE_POINT_PRESSURE__*1000*UC.kpa_to_bar

class criticalPt_PT():
    def __init__(self):
        """
        get the triple point pressure in bar and temperature in C
        In PyXSteam the critical pressure is returned in MPa and critical pressure in K
        https://pyxsteam.readthedocs.io/en/latest/pyXSteam.html?highlight=Constants.__CRITICAL_PRESSURE__#pyxsteam-constants-module
        """
        self.t=Constants.__CRITICAL_TEMPERATURE__+Constants.__ABSOLUTE_ZERO_CELSIUS__
        self.p=Constants.__CRITICAL_PRESSURE__*1000*UC.kpa_to_bar

class satProps():
    """
    For storage and retrieval of saturated properties at a given P&T
    """
    def __init__(self):
        self.tsat = None
        self.psat = None  # storage is in bar b/c pyXSteam
        self.hf = None
        self.hg = None
        self.hgf = None
        self.uf = None
        self.ug = None
        self.ugf = None
        self.sf = None
        self.sg = None
        self.sgf = None
        self.vf = None
        self.vg = None
        self.vgf = None

    def set(self, vals):
        self.tsat, self.psat, self.uf, self.ug, self.hf, self.hg, self.sf, self.sg, self.vf, self.vg = vals
        self.ugf = self.ug - self.uf
        self.hgf = self.hg - self.hf
        self.sgf = self.sg - self.sf
        self.vgf = self.vg - self.vf

    def get(self):
        return (self.tsat, self.psat, self.hf, self.hg, self.hgf, self.sf, self.sg, self.sgf, self.vf, self.vg, self.vgf)

    def getTextOutput(self, SI=True):
        """
        Sets the self.txtOut string for display.
        :param SI:
        :return: self.txtOut
        """
        if SI is False:
            P = self.psat*UC.bar_to_psi
            PUnits="psi"
            T =UC.C_to_F(self.tsat)
            TUnits="F"
            hf=self.hf*UC.kJperkg_to_BTUperlb
            hg=self.hg*UC.kJperkg_to_BTUperlb
            HUnits="BTU/lb"
            sf=self.sf*UC.kJperkgK_to_BTUperlbR
            sg=self.sg*UC.kJperkgK_to_BTUperlbR
            SUnits="BTU/lb*R"
            vf=self.vf*UC.m3perkg_to_ft3perlb
            vg=self.vg*UC.m3perkg_to_ft3perlb
            VUnits="ft^3/lb"
        else:
            P = self.psat
            PUnits = "bar"
            T = self.tsat
            TUnits = "C"
            hf = self.hf
            hg = self.hg
            HUnits = "kJ/kg"
            sf = self.sf
            sg = self.sg
            SUnits = "kJ/kg*K"
            vf = self.vf
            vg = self.vg
            VUnits = "m^3/kg"

        self.txtOut = "psat = {:0.2f} {}, tsat = {:0.2f} {}".format(P, PUnits, T, TUnits)
        self.txtOut += "\nhf = {:0.2f} {}, hg = {:0.2f} {}".format(hf, HUnits,hg,HUnits)
        self.txtOut += "\nsf = {:0.2f} {}, sg = {:0.2f} {}".format(sf,SUnits,sg,SUnits)
        self.txtOut += "\nvf = {:0.4f} {}, vg = {:0.4f} {}".format(vf,VUnits,vg,VUnits)
        return self.txtOut
    
class stateProps():
    """
    for storage and retrieval of a thermodynamic state
    T (C), P (kPa), u(kJ/kg), h (kJ/kg), s (kJ/kg*K), v (m^3/kg), x (dimensionless)
    """
    def __init__(self):
        self.name = None
        self.t = None
        self.p = None
        self.u = None
        self.h = None
        self.s = None
        self.v = None
        self.x = None
        self.region = None

    def getVal(self, name='T', SI=True ):
        if SI:
            uCF=1
            hCF=1
            sCF=1
            vCF=1
            pCF=1
        else:
            uCF=UC.kJperkg_to_BTUperlb
            hCF=UC.kJperkg_to_BTUperlb
            sCF=UC.kJperkgK_to_BTUperlbR
            vCF=UC.m3perkg_to_ft3perlb
            pCF=UC.kpa_to_psi

        n=name.lower()
        if n == 't':
            return self.t if SI else UC.C_to_F(self.t)
        if n == 'h':
            return self.h*hCF
        if n == 's':
            return self.s*sCF
        if n == 'v':
            return self.v*vCF
        if n == 'p':
            return self.p*pCF

    def print(self):
        if self.name is not None:
            print(self.name)
        if self.x is None or self.x < 0.0:
            print('Region: compressed liquid')
            print('p = {:0.2f} bar'.format(self.p))
            print('h = {:0.2f} kJ/kg'.format(self.h))
        else:
            print('Region: ', self.region)
            print('p = {:0.2f} bar'.format(self.p))
            print('T = {:0.1f} degrees C'.format(self.t))
            print('h = {:0.2f} kJ/kg'.format(self.h))
            print('s = {:0.4f} kJ/(kg K)'.format(self.s))
            print('v = {:0.6f} m^3/kg'.format(self.v))
            print('x = {:0.4f}'.format(self.x))
        print()

class StateDataForPlotting:
    """
    I'm making this class for easy storage of data for plotting.
    """
    def __init__(self):
        self.t = []
        self.p = []
        self.u = []
        self.h = []
        self.s = []
        self.v = []

    def clear(self):
        self.t.clear()
        self.p.clear()
        self.u.clear()
        self.h.clear()
        self.s.clear()
        self.v.clear()

    def addPt(self, vals):
        """
        adds a thermodynamic state point to the list
        :param vals: a list or tuple with T, P, u, h, s, v in that order
        :return:
        """
        T, P, u, h, s, v = vals
        self.t.append(T)
        self.p.append(P)
        self.u.append(u)
        self.h.append(h)
        self.s.append(s)
        self.v.append(v)

    def addStatePt(self, st):
        self.addPt((st.t, st.p,st.u,st.h,st.s,st.v))

    def getAxisLabel(self, W='T', SI=True):
        w=W.lower()
        if w == 't':
            return r'T $\left(^oC\right)$' if SI else r'T $\left(^oF\right)$'
        if w == 'h':
            return r'h $\left(\frac{kJ}{kg}\right)$' if SI else r'h $\left(\frac{BTU}{lb}\right)$'
        if w == 's':
            return r's $\left(\frac{kJ}{kg\cdot K}\right)$' if SI else r's $\left(\frac{BTU}{lb\cdot ^oR}\right)$'
        if w == 'v':
            return r'v $\left(\frac{m^3}{kg}\right)$' if SI else r'v $\left(\frac{ft^3}{lb}\right)$'
        if w == 'p':
            return r'P $\left(kPa\right)$' if SI else r'P (psi)'

    def getDataCol(self, W='T', SI=True):

        if SI:
            uCF=1
            hCF=1
            sCF=1
            vCF=1
            pCF=1
        else:
            uCF=UC.kJperkg_to_BTUperlb
            hCF=UC.kJperkg_to_BTUperlb
            sCF=UC.kJperkgK_to_BTUperlbR
            vCF=UC.m3perkg_to_ft3perlb
            pCF=UC.kpa_to_psi

        w=W.lower()
        if w=='t':
            return self.t if SI else [UC.C_to_F(t) for t in self.t]
        if w == 'u':
            return np.array(self.u) * uCF
        if w=='h':
            return np.array(self.h)*hCF
        if w=='s':
            return np.array(self.s)*sCF
        if w=='v':
            return np.array(self.v)*vCF
        if w=='p':
            return np.array(self.p)*pCF

class Steam_SI:
    def __init__(self, P=None, T=None, x=None, v=None, h=None, u=None, s=None, name=None):
        """
        This is a general steam class for sub-critical (i.e., superheated, subcooled and saturated) properties of steam.
        The user may specify any two properties to calculate all other properties of the steam.
        Note: we have 7 properties, but only can specify two of them.  Combinations=7!/(2!5!)=42
        But, since order of specifying properties does not matter, I really get only 21

        I handle all cases in self.getState

        :param P: Pressure (bar)
        :param T: Temperature (C)
        :param x: Quality
        :param v: Specific Volume (kg/m^3)
        :param h: Enthalpy (kJ/kg)
        :param u: Internal Energy (kJ/kg)
        :param s: Entropy (kJ/(kg*K))
        :param name:
        """
        self.steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
        self.satProps = satProps()
        self.state = stateProps()
        self.state.p = P  # pressure - bar
        self.state.t = T  # Temperature - degrees C
        self.state.x = x  # quality (a value between 0 and 1)
        self.state.v = v  # specific volume - m^3/kg
        self.state.h = h  # enthalpy - kJ/kg
        self.state.u = u  # internal energy - kJ/kg
        self.state.s = s  # entropy - kJ/(kg K)
        self.state.name = name  # a useful identifier
        self.state.region = None  # 'superheated' or 'saturated'
        self.RW=UC.R/UC.MW_Water #water gas constant kJ/(kg*K)
        self.getState(P=P,T=T,x=x,v=v,h=h, u=u ,s=s)

    def getsatProps_p(self, p):
        """
        Given a pressure, calculate the saturated properties for that isobar
        :param p:
        :return:
        """
        self.satProps.tsat = self.steamTable.tsat_p(p)
        self.satProps.psat = p
        self.satProps.vf = self.steamTable.vL_p(p)
        self.satProps.vg = self.steamTable.vV_p(p)
        self.satProps.hf = self.steamTable.hL_p(p)
        self.satProps.hg = self.steamTable.hV_p(p)
        self.satProps.uf = self.steamTable.uL_p(p)
        self.satProps.ug = self.steamTable.uV_p(p)
        self.satProps.sf = self.steamTable.sL_p(p)
        self.satProps.sg = self.steamTable.sV_p(p)
        self.satProps.vgf = self.satProps.vg-self.satProps.vf
        self.satProps.hgf = self.satProps.hg-self.satProps.hf
        self.satProps.sgf = self.satProps.sg-self.satProps.sf
        self.satProps.ugf = self.satProps.ug-self.satProps.uf
        return dc(self.satProps)

    def getsatProps_t(self, t):
        """
        Given a temperature, calculate the saturation pressure and then
        calculate all other saturated properties
        :param t:
        :return:
        """
        self.satProps.tsat=t
        self.satProps.psat=self.steamTable.psat_t(t)
        self.getsatProps_p(self.psat)
        return dc(self.satProps)

    def calcState_1Phase(self):
        """
        Given that I have T&P, find the other properties
        p should be in bar
        t should be in C
        :return:
        """
        self.state.u = self.steamTable.u_pt(self.state.p, self.state.t)
        self.state.h = self.steamTable.h_pt(self.state.p, self.state.t)
        self.state.s = self.steamTable.s_pt(self.state.p, self.state.t)
        self.state.v = self.steamTable.v_pt(self.state.p, self.state.t)
        self.state.x = 1.0 if self.state.t > self.steamTable.tsat_p(self.state.p) else 0.0
        pass

    def calcState_2Phase(self):
        """
        Given P and x, find all other properties
        :return:
        """
        if self.state.x == 0.0:
            self.state.region = "saturated liquid"
        elif self.state.x ==1.0:
            self.state.region = "saturated vapor"
        else:
            self.state.region = "two-phase"

        self.state.u = self.satProps.uf + self.state.x * self.satProps.ugf
        self.state.h = self.satProps.hf + self.state.x * self.satProps.hgf
        self.state.s = self.satProps.sf + self.state.x * self.satProps.sgf
        self.state.v = self.satProps.vf + self.state.x * self.satProps.vgf
        
    def between(self, x, xLow,xHigh):
        """
        This is just for convenience when finding if subcooled, superheated or saturated
        :param x: a thermodynamic property
        :param xLow: low value
        :param xHigh: high value
        :return: boolean
        """
        if x<xLow: return False
        if x>xHigh: return False
        return True
    
    def clamp(self, x, xLow, xHigh):
        """
        A convenience function to ensure a varible is within bounds
        :param x:
        :param xLow:
        :param xHigh:
        :return:
        """
        if x<xLow: return xLow
        if x>xHigh: return xHigh
        return x
    
    def getState(self, P=None, T=None, x=None, v=None, u=None, h=None, s=None, name=None):
        """
        Calculates the thermodynamic state variables based on specified input values.
        I have thermodynamic variables:  P, T, v, h, u, s and x (7 things) for which I am choosing two.
        Possible number of permutations:  7!/5! =42.
        But, order of the two things does not matter, so 42/2=21
        PT, Pv, Ph, Pu, Ps, Px (6)
        Tv, Th, Tu, Ts, Tx (5)
        vh, vu, vs, vx (4)
        hu, hs, hx (3)
        us, ux (2)
        sx (1)
        Total of 21 cases to deal with.  I will attack them in the order shown above.  In general,
        I find P&T for each case and then calculate all other properties from these.
        """
        if name is not None:
            self.state.name = name

        # Step 1: select the proper case from the 21 and encode it in a two letter string
        #region select case
        case=None
        if P is not None:  # pressure is specified
            if T is not None:  # 1
                case="pt"
            elif v is not None:  # 2
                case="pv"
            elif h is not None:  # 3
                case="ph"
            elif u is not None:  # 4
                case="pu"
            elif s is not None:  # 5
                case="ps"
            elif x is not None:  # 6
                case="px"
        elif case is None and T is not None:   #temperature is specified
            if x is not None:  # 7
                case="tx"
            elif v is not None:  # 8
                case="tv"
            elif u is not None:  # 9
                case = "tu"
            elif h is not None:  # 10
                case="th"
            elif s is not None:  # 11
                case="ts"
        elif case is None and x is not None:  # quality is specified
            if v is not None:  # 12
                case ="xv"
            elif u is not None:  # 13
                case = "xu"
            elif h is not None:  # 14
                case="xh"
            elif s is not None:  # 15
                case="xs"
        elif case is None and v is not None:  # quality is specified
            if h is not None:  # 16
                case="vh"
            elif u is not None:  # 17
                case = "vu"
            elif s is not None:  # 18
                case="vs"
        elif case is None and h is not None:  # enthalpy is specified
            if s is not None:  # 19
                case="hs"
            elif u is not None:  # 20
                case = "hu"
        elif case is None and s is not None:  # enthalpy is specified
            case = "su"  # 21
            
        if case is None:
            return self.state
        #endregion
        
        case=case.lower()  # probably unnecessary

        # Step 2: select the proper case from the 21.  Note that PT is the same as TP etc.
        if case.__contains__("p"):
            self.state.p = P
            self.getsatProps_p(self.state.p)
            # case 1:  pt or tp
            if case.__contains__("t"):
                self.state.t=T
                self.satProps.tsat = round(self.satProps.tsat, 3)  # I will compare at 3 three decimal places
                # compare T to tsat
                if T < self.satProps.tsat or T > self.satProps.tsat or self.satProps.tsat is "nan":
                    self.region = "sub-cooled liquid" if T < self.satProps.tsat else "super-heated vapor"
                    self.calcState_1Phase()
                else:  # this is ambiguous since at saturated temperature
                    self.state.x = 1.0
                    self.calcState_2Phase()
            # case 2: pv or vp
            elif case.__contains__("v"):
                self.state.v = v
                self.satProps.vf = round(self.satProps.vf, 5)
                self.satProps.vg = round(self.satProps.vg, 3)
                # compare v to vf and vg
                if self.state.v < self.satProps.vf or self.state.v > self.satProps.vg:
                    self.state.region = "sub-cooled liquid" if self.state.v < self.satProps.vf else "super-heated vapor"
                    # since I can't find properties using v, I will use fsolve to find T
                    dt = 1.0 if self.state.v > self.satProps.vg else -1.0
                    fn1 = lambda T: self.state.v - self.steamTable.v_pt(self.state.p, T)
                    self.state.t = fsolve(fn1, [self.satProps.tsat + dt])[0]
                    # now use P and T
                    self.calcState_1Phase()
                else:  # two-phase
                    # first calculate quality
                    self.state.x = (self.state.v - self.satProps.vf) / (self.satProps.vgf)
                    self.state.t = self.satProps.tsat
                    self.calcState_2Phase()
            # case 3 pu or up
            elif case.__contains__('u'):
                self.state.u = u
                # compare u to uf and ug
                if self.state.u < self.satProps.uf or self.state.u > self.satProps.ug:
                    self.state.region = "sub-cooled liquid" if self.state.u < self.satProps.uf else "super-heated vapor"
                    # since I can't find properties using u, I will use fsolve to find T
                    dt = 1.0 if self.state.u > self.satProps.ug else -1.0
                    fn3 = lambda T: self.state.u - self.steamTable.u_pt(self.state.p, T)
                    self.state.t = fsolve(fn3, [self.satProps.tsat + dt])[0]
                    # now use P and T
                    self.calcState_1Phase()
                else:  # two-phase
                    # first calculate quality
                    self.state.x = (self.state.u - self.satProps.uf) / (self.satProps.ugf)
                    self.state.t = self.satProps.tsat
                    self.calcState_2Phase()
            # case 4 ph or hp
            elif case.__contains__('h'):
                self.state.h = h
                self.getsatProps_p(self.state.p)
                # compare h to hf and hg
                if self.state.h < self.satProps.hf or self.state.h > self.satProps.hg:
                    self.state.region = "sub-cooled liquid" if self.state.h < self.satProps.hf else "super-heated vapor"
                    self.state.t = self.steamTable.t_ph(self.state.p, self.state.h)
                    # now use P and T
                    self.calcState_1Phase()
                else:  # two-phase
                    # first calculate quality
                    self.state.x = (self.state.h - self.satProps.hf) / (self.satProps.hgf)
                    self.state.t = self.satProps.tsat
                    self.calcState_2Phase()
            # case 5 ps or sp
            elif case.__contains__('s'):
                self.state.s = s
                self.getsatProps_p(self.state.p)
                # compare s to sf and sg
                if self.state.s < self.satProps.sf or self.state.s > self.satProps.sg:
                    self.state.region = "sub-cooled liquid" if self.state.s < self.satProps.sf else "super-heated vapor"
                    self.state.t = self.steamTable.t_ps(self.state.p, self.state.s)
                    # now use P and T
                    self.calcState_1Phase()
                else:  # two-phase
                    # first calculate quality
                    self.state.t=self.satProps.tsat
                    self.state.x = (self.state.s - self.satProps.sf) / (self.satProps.sgf)
                    self.calcState_2Phase()
            # case 6 px or xp
            elif case.__contains__('x'):
                self.state.x = x
                self.getsatProps_p(self.state.p)
                self.state.t = self.satProps.tsat
                self.state.x = self.clamp(self.state.x, 0.0, 1.0)
                self.calcState_2Phase()
        elif case.__contains__('t'):
            self.state.t=T
            self.satProps=self.getsatProps_t(self.state.t)
            # case 7:  tv or vt
            if case.__contains__('v'):
                self.state.v = v
                self.satProps.vf = round(self.satProps.vf, 5)
                self.satProps.vg = round(self.satProps.vg, 3)
                # compare v to vf and vg
                if self.state.v < self.satProps.vf or self.state.v > self.satProps.vg:
                    self.state.region = "sub-cooled liquid" if self.state.v < self.satProps.vf else "super-heated vapor"
                    # since I can't find properties using v, I will use fsolve to find P
                    dp = -0.1 if self.state.v > self.satProps.vg else 0.1
                    fn3 = lambda P: [self.state.v - self.steamTable.v_pt(P, self.state.t)]
                    self.state.p = fsolve(fn3, [self.state.satProps.psat + dp])[0]
                    # now use P and T
                    self.calcState_1Phase()
                else:  # two-phase
                    # first calculate quality
                    self.state.x = (self.state.v - self.satProps.vf) / (self.satProps.vgf)
                    self.state.p = self.state.satProps.psat
                    self.calcState_2Phase()
            # case 8:  tu or ut
            elif case.__contains__('u'):
                self.state.u = u
                self.getsatProps_t(self.state.t)
                # compare u to uf and ug
                if self.state.u < self.satProps.uf or self.state.u > self.satProps.ug:
                    self.state.region = "sub-cooled liquid" if self.state.u < self.satProps.uf else "super-heated vapor"
                    # since I can't find properties using u, I will use fsolve to find P
                    dp = 0.1 if self.state.u > self.satProps.ug else -0.1
                    fn8 = lambda P: self.state.u - self.steamTable.u_pt(self.state.t, P)
                    self.state.p = fsolve(fn8, [self.state.satProps.psat + dp])[0]
                    # now use P and T
                    self.calcState_1Phase()
                else:  # two-phase
                    # first calculate quality
                    self.state.x = (self.state.u - self.satProps.uf) / (self.satProps.ugf)
                    self.state.p = self.state.satProps.psat
                    self.calcState_2Phase()
            # case 9:  th or ht
            elif case.__contains__('h'):
                self.state.h = h
                self.getsatProps_t(self.state.t)
                # compare h to hf and hg
                if self.state.h < self.satProps.hf or self.state.h > self.satProps.hg:
                    self.state.region = "sub-cooled liquid" if self.state.h < self.satProps.hf else "super-heated vapor"
                    self.state.p = self.steamTable.p_th(self.state.t, self.state.h)
                    # now use P and T
                    self.calcState_1Phase()
                else:  # two-phase
                    # first calculate quality
                    self.state.p = self.state.satProps.psat
                    self.state.x = (self.state.h - self.satProps.hf) / (self.satProps.hgf)
                    self.calcState_2Phase()
            # case 10:  ts or st
            elif case.__contains__('s'):
                self.state.s = s
                self.getsatProps_t(self.state.t)
                # compare s to sf and sg
                if self.state.s < self.satProps.sf or self.state.s > self.satProps.sg:
                    self.state.region = "sub-cooled liquid" if self.state.s < self.satProps.sf else "super-heated vapor"
                    self.state.p = self.steamTable.p_ts(self.state.t, self.state.s)
                    # now use P and T
                    self.calcState_1Phase()
                else:  # two-phase
                    # first calculate quality
                    self.state.p = self.state.satProps.psat
                    self.state.x = (self.state.s - self.satProps.sf) / (self.satProps.sgf)
                    self.calcState_2Phase()
            # case 11:  tx or xt
            elif case.__contains__('x'):
                self.state.x = x
                self.getsatProps_t(self.state.t)
                self.state.p = self.state.satProps.psat
                self.state.x = self.clamp(self.state.x, 0.0, 1.0)
                self.calcState_2Phase()
        elif case.__contains__('v') :
            self.state.v=v
            # case 12:  vh or hv
            if case.__contains__('h'):
                self.state.h = h
                def fn12(P):
                    # could be single phase or two-phase, but both v&h have to match at same x
                    self.getsatProps_p(P)
                    if self.between(self.state.h, self.satProps.hf, self.satProps.hg):
                        self.state.x = (self.state.h - self.satProps.hf) / self.satProps.hgf
                        return self.state.v - (self.satProps.vf + self.state.x * self.satProps.vgf)
                    # could be single phase
                    return self.state.v - self.steamTable.v_ph(P, self.state.h)

                self.state.p = fsolve(fn12, [1.0])[0]
                self.state.t = self.steamTable.t_ph(self.state.p, self.state.h)
                self.getsatProps_p(self.state.p)
                self.satProps.vf = round(self.satProps.vf, 5)
                self.satProps.vg = round(self.satProps.vg, 3)
                # compare v to vf and vg
                if self.state.v < self.satProps.vf or self.state.v > self.satProps.vg:
                    self.state.region = "sub-cooled liquid" if self.state.v < self.satProps.vf else "super-heated vapor"
                    # now use P and T
                    self.calcState_1Phase()
                else:  # two-phase
                    # first calculate quality
                    self.state.x = (self.state.v - self.satProps.vf) / (self.satProps.vgf)
                    self.calcState_2Phase()
            # case 13:  vu or uv
            elif case.__contains__('u'):
                self.state.u = u
                # use fsolve to find P&T at this v & u
                def fn13(PT):
                    self.getsatProps_p(PT[0])
                    if self.between(self.state.u, self.satProps.uf, self.satProps.ug):
                        self.state.t = self.satProps.tsat
                        self.state.x = (self.state.u - self.satProps.uf) / self.satProps.ugf
                        return [self.state.v - (self.satProps.vf + self.state.x * self.satProps.vgf), 0]
                    return [self.state.v - self.steamTable.v_pt(PT[0], PT[1]),
                            self.state.u - self.steamTable.u_pt(PT[0], PT[1])]
                props = fsolve(fn13, [1, 100])
                self.state.p = props[0]
                self.state.t = props[1]
                self.getsatProps_p(self.state.p)
                # compare u to uf and ug
                if self.state.u < self.satProps.uf or self.state.u > self.satProps.ug:
                    self.state.region = "sub-cooled liquid" if self.state.u < self.satProps.uf else "super-heated vapor"
                    # now use P and T
                    self.calcState_1Phase()
                else:  # two-phase
                    # first calculate quality
                    self.state.x = (self.state.u - self.satProps.uf) / (self.satProps.ugf)
                    self.state.p = self.satProps.psat
                    self.state.t = self.satProps.tsat
                    self.calcState_2Phase()
            # case 14:  vs os sv
            elif case.__contains__('s'):
                self.state.s = s
                def fn14(PT):
                    self.getsatProps_p(PT[0])
                    if self.between(self.state.s, self.satProps.sf, self.satProps.sg):
                        self.state.x = (self.state.s - self.satProps.sf) / self.satProps.sgf
                        return [self.state.v - self.satProps.vf - self.state.x * self.satProps.vgf, 0.0]
                    return [self.state.v - self.steamTable.v_pt(PT[0], PT[1]),
                            self.state.s - self.steamTable.s_pt(PT[0], PT[1])]
                props = fsolve(fn14, [1, 100])
                self.state.p = props[0]
                self.state.t = props[1]
                self.getsatProps_p(self.state.p)
                # compare s to sf and sg
                if self.state.s < self.satProps.sf or self.state.s > self.satProps.sg:
                    self.state.region = "sub-cooled liquid" if self.state.s < self.satProps.sf else "super-heated vapor"
                    # now use P and T
                    self.calcState_1Phase()
                else:  # two-phase
                    # first calculate quality
                    self.state.p = self.state.satProps.psat
                    self.state.t = self.satProps.tsat
                    self.state.x = (self.state.s - self.satProps.sf) / (self.satProps.sgf)
                    self.calcState_2Phase()
            # case 15:  vx or xv
            elif case.__contains__('x'):
                self.state.x = self.clamp(x,0.0,1.0)
                def fn15(p):
                    self.getsatProps_p(p)
                    return self.state.v - (self.satProps.vf + self.state.x * self.satProps.vgf)
                props = fsolve(fn15, [1])
                self.state.p = props[0]
                self.state.t = self.satProps.tsat
                self.calcState_2Phase()
        elif case.__contains__('h'):
            self.state.h=h
            # case 16:  hu or uh
            if case.__contains__('u'):
                self.state.u = u
                # use fsolve to find P&T at this h & u
                def fn16(PT):
                    self.getsatProps_p(PT[0])
                    if self.between(self.state.u, self.satProps.uf, self.satProps.ug):
                        self.state.x = (self.state.u - self.satProps.uf) / self.satProps.ugf
                        return [self.state.h - self.satProps.hf - self.state.x * self.satProps.hgf, 0.0]
                    return [self.state.h - self.steamTable.h_pt(PT[0], PT[1]),
                            self.state.u - self.steamTable.u_pt(PT[0], PT[1])]
                props = fsolve(fn16, [1, 100])
                self.state.p = props[0]
                self.state.t = props[1]
                self.getsatProps_p(self.state.p)
                # compare u to uf and ug
                if self.state.u < self.satProps.uf or self.state.u > self.satProps.ug:
                    self.state.region = "sub-cooled liquid" if self.state.u < self.satProps.uf else "super-heated vapor"
                    # now use P and T
                    self.calcState_1Phase()
                else:  # two-phase
                    # first calculate quality
                    self.state.x = (self.state.u - self.satProps.uf) / (self.satProps.ugf)
                    self.state.p = self.satProps.psat
                    self.state.t = self.satProps.tsat
                    self.calcState_2Phase()
            # case 17:  hs or sh
            elif case.__contains__('s'):
                self.state.s = s
                def fn17(PT):
                    self.getsatProps_p(PT[0])
                    if self.between(self.state.s, self.satProps.sf, self.satProps.sg):
                        self.state.x = (self.state.s - self.satProps.sf) / self.satProps.sgf
                        return [self.state.h - self.satProps.hf - self.state.x * self.satProps.hgf, 0.0]
                    return [self.state.h - self.steamTable.h_pt(PT[0], PT[1]),
                            self.state.s - self.steamTable.s_pt(PT[0], PT[1])]
                props = fsolve(fn17, [1, 100])
                self.state.p = props[0]
                self.state.t = props[1]
                self.getsatProps_p(self.state.p)
                # compare s to sf and sg
                if self.state.s < self.satProps.sf or self.state.s > self.satProps.sg:
                    self.state.region = "sub-cooled liquid" if self.state.s < self.satProps.sf else "super-heated vapor"
                    # now use P and T
                    self.calcState_1Phase()
                else:  # two-phase
                    # first calculate quality
                    self.state.p = self.state.satProps.psat
                    self.state.t = self.satProps.tsat
                    self.state.x = (self.state.s - self.satProps.sf) / (self.satProps.sgf)
                    self.calcState_2Phase()
            # case 18:  hx or xh
            elif case.__contains__('x'):
                self.state.x = self.clamp(x,0.0,1.0)
                def fn18(p):
                    self.getsatProps_p(p)
                    return self.state.h - (self.satProps.hf + self.state.x * self.satProps.hgf)
                props = fsolve(fn18, [1])
                self.state.p = props[0]
                self.getsatProps_p(self.state.p)
                self.state.t = self.satProps.tsat
                self.calcState_2Phase()
        elif case.__contains__("u"):
            self.state.u = u
            # case 19:  us or su
            if case.__contains__('s'):
                self.state.s = s
                def fn19(PT):
                    self.getsatProps_p(PT[0])
                    if self.between(self.state.s, self.satProps.sf, self.satProps.sg):
                        self.state.x = (self.state.s - self.satProps.sf) / self.satProps.sgf
                        return [self.state.u - self.satProps.uf - self.state.x * self.satProps.ugf, 0.0]
                    return [self.state.u - self.steamTable.u_pt(PT[0], PT[1]),
                            self.state.s - self.steamTable.s_pt(PT[0], PT[1])]
                props = fsolve(fn19, [1, 100])
                self.state.p = props[0]
                self.state.t = props[1]
                self.getsatProps_p(self.state.p)
                # compare s to sf and sg
                if self.state.s < self.satProps.sf or self.state.s > self.satProps.sg:
                    self.state.region = "sub-cooled liquid" if self.state.s < self.satProps.sf else "super-heated vapor"
                    # now use P and T
                    self.calcState_1Phase()
                else:  # two-phase
                    # first calculate quality
                    self.state.p = self.satProps.psat
                    self.state.t = self.satProps.tsat
                    self.state.x = (self.state.s - self.satProps.sf) / (self.satProps.sgf)
                    self.calcState_2Phase()
            # case 20:  ux or xu
            elif case.__contains__('x'):
                self.state.x = self.clamp(x,0.0,1.0)
                self.state.region = "two-phase"
                def fn20(p):
                    self.getsatProps_p(p)
                    return self.state.h - (self.satProps.hf + self.state.x * self.satProps.hgf)
                props = fsolve(fn20, [1])
                self.state.p = props[0]
                self.getsatProps_p(self.state.p)
                self.state.t = self.satProps.tsat
                self.calcState_2Phase()
        elif case.__contains__('s'):
            self.state.s=s
            # case 21:  sx or xs
            if case.__contains__('x'):
                self.state.x = self.clamp(x, 0.0,1.0)
                def fn21(p):
                    self.getsatProps_p(p)
                    return self.state.h - (self.satProps.hf + self.state.x * self.satProps.hgf)
                props = fsolve(fn21, [1])
                self.state.p = props[0]
                self.getsatProps_p(self.state.p)
                self.state.t = self.satProps.tsat
                self.calcState_2Phase()
        return dc(self.state)

    def igl_v(self):
        # ideal gas law V=RT/P-> v=(R/MW)*(T+273)/(P)
        # calculate a column of specific volume for the superheated water table using ideal gas
        return self.RW * (self.t + 273) / self.p

    def print(self):
        if self.state.name is not None:
            print('Name: {}'.format(self.state.name))
        if self.state.region is not None:
            print('Region: {}'.format(self.state.region))
        if self.state.p is not None:
            print('p = {:.2f} bar'.format(self.state.p))
        if self.state.t is not None:
            print('T = {:.1f} degrees C'.format(self.state.t))
        if self.state.u is not None:
            print('u = {:.2f} kJ/kg'.format(self.state.u))
        if self.state.h is not None:
            print('h = {:.2f} kJ/kg'.format(self.state.h))
        if self.state.s is not None:
            print('s = {:.4f} kJ/(kg K)'.format(self.state.s))
        if self.state.v is not None:
            print('v = {:.6f} m^3/kg'.format(self.state.v))
        if self.state.x is not None:
            print('x = {:.4f}'.format(self.state.x))
        print('')
#endregion

#region function definitions
def main():

    inlet=Steam_SI(P=80, x=1.0, name="Turbine Inlet")
    inlet.print()
    h1 = inlet.state.h
    s1 = inlet.state.s
    print(h1,s1,'\n')
    steam=XSteam(XSteam.UNIT_SYSTEM_MKS)

    outlet=Steam_SI(P=1, s=inlet.state.s, name='Turbine Exit')
    #notice -  s=inlet.s
    outlet.print()

    another=Steam_SI(P=85.75, h=2050, name='State 3')
    another.print()

    yetanother=Steam_SI(P=89, h=4125, name='State 4')
    yetanother.print()

    g1=Steam_SI(P=8, x=1.0, name='Gap1')
    g1.print()
    g2=Steam_SI(P=g1.state.p, s=g1.state.s * 1.0012, name='Gap2')
    g2.print()
    g2=Steam_SI(P=g1.state.p, s=6.6699, name='Gap3')
    g2.print()




    #uncommenting the next two lines will cause a ValueError error
    #final1 = State(7000, T=250, name='State 5')
    #final1.print()
#endregion

#region function calls
if __name__ == "__main__":
   main()
#endregion
