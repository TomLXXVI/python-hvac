from .. import Quantity
from ..fluids import Fluid, FluidState

Q_ = Quantity


# noinspection PyCallingNonCallable
class ReciprocatingCompressor:
    """
    Model of an ideal reciprocating compressor.
    The clearance volume efficiency and polytropic work relations assume that
    there are no pressure drops across the valves and that the compression and
    expansion processes are ideal and polytropic. These relations give an
    upper limit for the performance.
    """
    def __init__(
        self,
        V_dis: Quantity,
        C: Quantity,
        speed: Quantity,
        n: float,
        refrigerant_type: Fluid
    ) -> None:
        """
        Parameters
        ----------
        V_dis: Quantity
            Displacement volume.
        C: Quantity
            Clearance volume fraction.
        speed: Quantity
            Compressor speed.
        n: float
            Polytropic exponent.
        refrigerant_type: Fluid
            Type of refrigerant.
        """
        self.V_dis = V_dis
        self.C = C
        self.speed = speed
        self.n = n
        self.Refrigerant = refrigerant_type

        self._Pc: Quantity = None
        self._Pe: Quantity = None
        self._v_suc: Quantity = None

    @property
    def Pc(self) -> Quantity:
        """Get condenser pressure."""
        return self._Pc

    @Pc.setter
    def Pc(self, v: Quantity) -> None:
        """Set condenser pressure."""
        self._Pc = v

    @property
    def Pe(self) -> Quantity:
        """Get evaporator pressure."""
        return self._Pe

    @Pe.setter
    def Pe(self, v: Quantity) -> None:
        """Set evaporator pressure."""
        self._Pe = v

    @property
    def v_suc(self) -> Quantity:
        """Get specific volume of suction gas at compressor inlet."""
        return self._v_suc

    @v_suc.setter
    def v_suc(self, v: Quantity) -> None:
        """Set specific volume of suction gas at compressor inlet."""
        self._v_suc = v

    @property
    def eta_vol(self) -> Quantity:
        """Get clearance volumetric efficiency at given working conditions."""
        eta_vol = 1 + self.C - self.C * (self.Pc / self.Pe) ** (1 / self.n)
        return eta_vol

    @property
    def m_dot(self) -> Quantity:
        """Get mass flow rate of refrigerant at given working conditions."""
        m_dot = self.eta_vol * self.speed * self.V_dis / self.v_suc
        return m_dot

    @property
    def Wc_dot(self) -> Quantity:
        """Get compressor power at given working conditions."""
        e = self.n / (self.n - 1)
        a = self.eta_vol * self.speed * self.V_dis * self.Pe * e
        b = (self.Pc / self.Pe) ** (1 / e) - 1
        Wc_dot = a * b
        return Wc_dot

    @property
    def suction_gas(self) -> FluidState:
        """Get state of suction gas at compressor inlet."""
        return self.Refrigerant(P=self.Pe, rho=1 / self.v_suc)

    @suction_gas.setter
    def suction_gas(self, refrigerant: FluidState) -> None:
        """Set state of suction gas at compressor inlet."""
        self.Pe = refrigerant.P
        self.v_suc = 1 / refrigerant.rho

    @property
    def v_dis(self) -> Quantity:
        """Get specific volume of discharge gas at compressor outlet."""
        # v_dis = (self.Pe * (self.v_suc ** self.n) / self.Pc) ** (1 / self.n)
        Pe = self.Pe.to('bar').m
        Pc = self.Pc.to('bar').m
        v_dis = (Pe / Pc) ** (1 / self.n) * self.v_suc
        return v_dis

    @property
    def discharge_gas(self) -> FluidState:
        """Get state of discharge gas at compressor outlet."""
        return self.Refrigerant(P=self.Pc, rho=1 / self.v_dis)

    @property
    def Qc_dot(self) -> Quantity:
        """Get cooling capacity of compressor at given working conditions."""
        condenser_out = self.Refrigerant(P=self.Pc, x=Q_(0, 'frac'))
        evaporator_in = self.Refrigerant(P=self.Pe, h=condenser_out.h)
        evaporator_out = self.Refrigerant(P=self.Pe, rho=1 / self.v_suc)
        Qc_dot = self.m_dot * (evaporator_out.h - evaporator_in.h)
        return Qc_dot

    @property
    def COP(self) -> Quantity:
        """Get COP of compressor under given working conditions."""
        return self.Qc_dot / self.Wc_dot
