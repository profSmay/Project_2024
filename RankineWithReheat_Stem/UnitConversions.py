
class UnitConverter():
    def __init__(self):
        """
        This unit converter class is useful for the pipe network and perhaps other problems.
        The strategy is (number in current units)*(conversion factor)=(number desired units), for instance:
            1(ft)*(self.ft_to_m) = 1/3.28084 (m)
            1(in^2)*(self.in2_to_m2) = 1*(1/(12*3.28084))**2 (m^2)
        """
        pass

    """
    These constants can be used directly from the class without instantiating an object
    """
    #length conversions
    m_to_ft = 3.28084
    ft_to_m = 1 / m_to_ft
    in_to_m = ft_to_m / 12
    m_to_in = 1 / in_to_m
    mm_to_in = m_to_in/1000
    in_to_mm = 1/mm_to_in

    #area conversions
    m2_to_ft2 = m_to_ft**2
    ft2_to_m2 = ft_to_m ** 2
    in2_to_m2 = in_to_m ** 2
    m2_to_in2 = 1 / in2_to_m2

    #volume conversions
    m3_to_ft3 = m_to_ft**3
    ft3_to_m3 = ft_to_m ** 3
    m3_to_L = 1000
    ft3_to_L = ft3_to_m3 * m3_to_L
    L_to_ft3 = 1 / ft3_to_L

    g_SI = 9.80665  # m/s^2
    g_EN = 32.174  # 32.174 ft/s^2
    gc_EN = 32.174  # lbm*ft/lbf*s^2
    gc_SI = 1.0  # kg*m/N*s^2

    #mass/force
    lbf_to_kg = 1 / 2.20462
    kg_to_lbf = 1/lbf_to_kg
    lbf_to_N = lbf_to_kg * g_SI

    #pressure/head
    pa_to_psi = (1 / (lbf_to_N)) * in2_to_m2
    kpa_to_psi = 1000*pa_to_psi
    kpa_to_bar = 1/100.0
    bar_to_kpa = 1/kpa_to_bar
    psi_to_pa = 1/pa_to_psi
    psi_to_kpa = psi_to_pa/1000
    psi_to_bar = psi_to_kpa*kpa_to_bar
    bar_to_psi = 1/psi_to_bar
    mh2o_to_psi = 1000*g_SI*1.0*pa_to_psi #1000(kg/m^3)*9.81(m/s^2)*1(m) -> 9810Pa

    #viscosity
    PaS_to_slbPerSqft = pa_to_psi*144
    slbPerSqFt_to_PaS = 1/PaS_to_slbPerSqft

    #density
    kgperm3_to_lbperft3 = kg_to_lbf/m3_to_ft3
    lbperft3_to_kgperm3 = 1/kgperm3_to_lbperft3

    #delta T
    deltaK_to_deltaR = 9/5*(1.0)

    # Energy
    BTU_to_J=1055.06
    kJ_to_BTU = 1000/BTU_to_J
    BTU_to_kJ = 1/kJ_to_BTU
    kJperkg_to_BTUperlb = kJ_to_BTU/kg_to_lbf
    m3perkg_to_ft3perlb = m3_to_ft3/kg_to_lbf

    # Entropy
    kJperkgK_to_BTUperlbR = kJperkg_to_BTUperlb/deltaK_to_deltaR

    #Universal ideal gas constant
    R=8.31446261815324  #J/(mol*K). kJ/(kmol*K)

    #Water molecular weight
    MW_Water = 18.01528  #kg/kmol

    @classmethod  # a classmethod can be used directly from a class without needing to instantiate an object
    def viscosityEnglishToSI(cls, mu, toSI=True):
        """
        Converts between lb*s/ft^2 and Pa*s
        :param mu: the viscosity in english units
        :param toSI:  True assumes english in, False assumes SI in
        :return: the viscosity in Pa*s if toSI=True, lb*s/ft^2 if toSI=False
        """
        # (lb*s)/ft^2*((3.3 ft/m)^2)*(1kg/2.2lb)*(9.81m/s^2)->(Pa*s)
        cf = (1 / cls.ft2_to_m2) * (cls.lbf_to_kg) * cls.g_SI
        return mu * cf if toSI else mu / cf

    @classmethod  # a classmethod can be used directly from a class without needing to instantiate an object
    def densityEnglishToSI(cls, rho, toSI=True):
        """
        Converts between lb/ft^3 and kg/m^3
        :param rho: specific weight or density
        :param toSI:  True assumes english in, False assumes SI in
        :return: density in SI or EN
        """
        # (lb/ft^3)*((3.3ft/m)^3)*(1kg/2.2lb) -> kg/m^3
        cf = cls.lbf_to_kg / cls.ft3_to_m3
        return rho * cf if toSI else rho / cf

    @classmethod  # a classmethod can be used directly from a class without needing to instantiate an object
    def head_to_pressure(cls, h, rho, SI=True):
        """
        Convert from height of column of fluid to pressure in consistent units
        :param h: head in height of fluid (in or m)
        :return: pressure in (psi or Pa)
        """
        if SI:  # p = rho*g*h = g*cf
            cf = rho * cls.g_SI / cls.gc_SI  # kg*m/m^3*s^2
            return h * cf
        else:  # p = rho*g*h = g*cf (h in in)
            cf = rho * cls.g_EN / cls.gc_EN * (1 / 12) ** 2  # (lbm*ft/ft^3*s^2)(lbf*s^2/lbm*ft)(ft^2/in^2)
            return h * cf
        # convert m of water to psi
        # (m)*(3.3*12in/m)*rho(kg/m^3)*(2.2lb/kg)*(1m/(3.3*12in))^3
        psi = p * cls.rho * 2.2 / ((3.3 * 12) ** 2)
        return psi

    @classmethod  # a classmethod can be used directly from a class without needing to instantiate an object
    def m_to_psi(cls, h, rho):
        """
        For converting from height of fluid to psi
        :param h: height of fluid in m
        :param rho: density of fluid in kg/m^3
        :return: pressure in psi
        """
        return cls.head_to_pressure(h, rho) * cls.pa_to_psi

    @classmethod  # a classmethod can be used directly from a class without needing to instantiate an object
    def psi_to_m(cls, p, rho):
        """
        For converting from psi to height of fluid.
        first convert psi to pa
        :param p: pressure in psi
        :param rho: density of fluid in kg/m^3
        :return: height of fluid in m
        """
        pa = p / cls.pa_to_psi
        h = pa / (rho * cls.g_SI)
        return h

    @classmethod # a classmethod can be used directly from a class without needing to instantiate an object
    def C_to_F(cls, T=0):
        return 9/5*(T)+32

    @classmethod # a classmethod can be used directly from a class without needing to instantiate an object
    def F_to_C(cls, T=0):
        return 5/9*(T-32)

    @classmethod # a classmethod can be used directly from a class without needing to instantiate an object
    def K_to_R(cls, T=0):
        return cls.C_to_F(T-273.15)+459.67

    @classmethod
    def clamp(cls, x, low=0.0, high=1.0):
        if x < low: return low
        if x > high: return high
        return x