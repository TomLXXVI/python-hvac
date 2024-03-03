import numpy as np
from scipy import interpolate
from hvac import Quantity, UNITS
from hvac.fluid_flow.conduit import Duct
from ..abstract_fitting import AbstractFitting

Q_ = Quantity
u = UNITS


def convert_zeta(zeta1: float, duct1: Duct, duct2: Duct) -> float:
    # convert `zeta_1` referred to `duct_1` to `zeta` referred to `duct_2`
    pv1 = duct1.velocity_pressure.to('Pa').magnitude
    pv2 = duct2.velocity_pressure.to('Pa').magnitude
    zeta2 = zeta1 * pv1 / pv2
    return zeta2


class ReynoldsNumberCorrectionFactor:
    _Re = [1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 14.0]
    _k_R_on_W_050 = [1.40, 1.26, 1.19, 1.14, 1.09, 1.06, 1.04, 1.0]
    _k_R_on_W_075 = [2.0, 1.77, 1.64, 1.56, 1.46, 1.38, 1.30, 1.15]

    def __init__(self, R_on_W: float, Re: float):
        self.Re = Re * 1e-4
        if R_on_W <= 0.5:
            self._R_on_W = self._k_R_on_W_050
        elif R_on_W >= 0.75:
            self._R_on_W = self._k_R_on_W_075
        else:
            raise ValueError("R / W must be <= 0.5 or >= 0.75")

    def __call__(self) -> float:
        if self.Re >= 20.0:
            return 1.0
        else:
            k = np.interp(self.Re, self._Re, self._R_on_W)
            return k


class ElbowA7A(AbstractFitting):
    description = "Elbow, smooth radius, round"""
    _R_on_D = [0.5, 0.75, 1.0, 1.5, 2.0, 2.5]
    _zeta = [0.71, 0.33, 0.22, 0.15, 0.13, 0.12]
    _theta = [0, 20, 30, 45, 60, 75, 90, 110, 130, 150, 180]
    _k_theta = [0, 0.31, 0.45, 0.60, 0.78, 0.90, 1.00, 1.13, 1.20, 1.28, 1.40]

    def __init__(self, duct: Duct, R_on_D: Quantity, theta: Quantity = Q_(90.0, 'deg'), ID: str = ''):
        super().__init__(ID)
        self.R_on_D = R_on_D.to('frac').magnitude
        self.theta = theta.to('deg').magnitude
        self.Pv = duct.velocity_pressure.to('Pa').magnitude

    @property
    def zeta(self) -> float:
        zeta = np.interp(self.R_on_D, self._R_on_D, self._zeta)
        k_theta = np.interp(self.theta, self._theta, self._k_theta)
        return zeta * k_theta

    @property
    def pressure_drop(self) -> Quantity:
        return Q_(self.zeta * self.Pv, 'Pa')


class ElbowA7B(AbstractFitting):
    description = "Elbow, round, 3 to 5 piece, 90°"
    _R_on_D = [0.5, 0.75, 1.0, 1.5, 2.0]
    _zeta_5p = [float('nan'), 0.46, 0.33, 0.24, 0.19]
    _zeta_4p = [float('nan'), 0.5, 0.37, 0.27, 0.24]
    _zeta_3p = [0.98, 0.54, 0.42, 0.34, 0.33]

    def __init__(self, duct: Duct, radius: Quantity, no_pieces: int, ID: str = ''):
        super().__init__(ID)
        if 3 <= no_pieces <= 5:
            D = duct.cross_section.equivalent_diameter.to('m').magnitude
            radius = radius.to('m').magnitude
            self.R_on_D = radius / D
            self.no_pieces = no_pieces,
            self.Pv = duct.velocity_pressure.to('Pa').magnitude
        else:
            raise ValueError('number of pieces must be between 3..5')

    @property
    def zeta(self) -> float:
        if self.no_pieces == 5:
            zeta = np.interp(self.R_on_D, self._R_on_D, self._zeta_5p)
        elif self.no_pieces == 4:
            zeta = np.interp(self.R_on_D, self._R_on_D, self._zeta_4p)
        else:
            zeta = np.interp(self.R_on_D, self._R_on_D, self._zeta_3p)
        return zeta

    @property
    def pressure_drop(self) -> Quantity:
        return Q_(self.zeta * self.Pv, 'Pa')


class ElbowA7C(AbstractFitting):
    description = "Elbow, round, mitered"
    _theta = [20.0, 30.0, 45.0, 60.0, 75.0, 90.0]
    _zeta = [0.08, 0.16, 0.34, 0.55, 0.81, 1.2]

    def __init__(self, duct: Duct, theta: Quantity, ID: str = ''):
        super().__init__(ID)
        self.theta = theta.to('deg').magnitude
        self.Pv = duct.velocity_pressure.to('Pa').magnitude
        self.k_Re = ReynoldsNumberCorrectionFactor(R_on_W=0.0, Re=duct.reynolds_number)

    @property
    def zeta(self) -> float:
        zeta = np.interp(self.theta, self._theta, self._zeta)
        k_Re = self.k_Re()
        return zeta * k_Re

    @property
    def pressure_drop(self) -> Quantity:
        return Q_(self.zeta * self.Pv, 'Pa')


class ElbowA7D(AbstractFitting):
    description = "Elbow, rectangular, mitered"
    _H_on_W = np.array([0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 4.0, 6.0, 8.0])
    _theta = np.array([20.0, 30.0, 45.0, 60.0, 75.0, 90.0])
    _zeta = np.array([
        [0.08, 0.08, 0.08, 0.07, 0.07, 0.07, 0.06, 0.05, 0.05],
        [0.18, 0.17, 0.17, 0.16, 0.15, 0.15, 0.13, 0.12, 0.11],
        [0.38, 0.37, 0.36, 0.34, 0.33, 0.31, 0.27, 0.25, 0.24],
        [0.60, 0.59, 0.57, 0.55, 0.52, 0.49, 0.43, 0.39, 0.38],
        [0.89, 0.87, 0.84, 0.81, 0.77, 0.73, 0.63, 0.58, 0.57],
        [1.30, 1.30, 1.20, 1.20, 1.10, 1.10, 0.92, 0.85, 0.83]
    ])

    def __init__(self, duct: Duct, theta: Quantity, ID: str = ''):
        super().__init__(ID)
        H = duct.cross_section.height.to('m').magnitude
        W = duct.cross_section.width.to('m').magnitude
        self.Pv = duct.velocity_pressure.to('Pa').magnitude
        self.H_on_W = H / W
        self.theta = theta.to('deg').magnitude
        self.interp = interpolate.interp2d(self._H_on_W, self._theta, self._zeta)
        self.k_Re = ReynoldsNumberCorrectionFactor(R_on_W=0.0, Re=duct.reynolds_number)

    @property
    def zeta(self) -> float:
        zeta = self.interp(self.H_on_W, self.theta)
        k_Re = self.k_Re()
        return zeta[0] * k_Re

    @property
    def pressure_drop(self) -> Quantity:
        return Q_(self.zeta * self.Pv, 'Pa')


class ElbowA7E(AbstractFitting):
    description = "Elbow, rectangular, mitered with converging or diverging flow"
    _W1_on_W = np.array([0.6, 0.8, 1.2, 1.4, 1.6, 2.0])
    _H_on_W = np.array([0.25, 1.0, 4.0, 1e6])
    _zeta = np.array([
        [1.8, 1.4, 1.1, 1.1, 1.1, 1.1],
        [1.7, 1.1, 1.0, 0.95, 0.90, 0.84],
        [1.5, 1.0, 0.81, 0.76, 0.72, 0.66],
        [1.5, 1.0, 0.69, 0.63, 0.60, 0.55]
    ])

    def __init__(self, duct_down: Duct, duct_up: Duct, ID: str = ''):
        super().__init__(ID)
        W1 = duct_up.cross_section.width.to('m').magnitude
        W = duct_down.cross_section.width.to('m').magnitude
        self.Pv = duct_up.velocity_pressure.to('Pa').magnitude
        self.W1_on_W = W1 / W
        H = duct_down.cross_section.height.to('m').magnitude
        self.H_on_W = H / W
        self.interp = interpolate.interp2d(self._W1_on_W, self._H_on_W, self._zeta)
        self.k_Re = ReynoldsNumberCorrectionFactor(R_on_W=0.0, Re=duct_up.reynolds_number)

    @property
    def zeta(self) -> float:
        zeta = self.interp(self.W1_on_W, self.H_on_W)
        return zeta[0]

    @property
    def pressure_drop(self) -> Quantity:
        return Q_(self.zeta * self.Pv, 'Pa')


class ElbowA7F(AbstractFitting):
    description = "Elbow, rectangular, smooth radius without vanes"
    _H_on_W = np.array([0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0])
    _R_on_W = np.array([0.5, 0.75, 1.0, 1.5, 2.0])
    _zeta = np.array([
        [1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 1.0, 1.1, 1.1, 1.2, 1.2],
        [0.57, 0.52, 0.48, 0.44, 0.40, 0.39, 0.39, 0.40, 0.42, 0.43, 0.44],
        [0.27, 0.25, 0.23, 0.21, 0.19, 0.18, 0.18, 0.19, 0.20, 0.27, 0.21],
        [0.22, 0.20, 0.19, 0.17, 0.15, 0.14, 0.14, 0.15, 0.16, 0.17, 0.17],
        [0.20, 0.18, 0.16, 0.15, 0.14, 0.13, 0.13, 0.14, 0.14, 0.15, 0.15]
    ])

    def __init__(self, duct: Duct, radius: Quantity, ID: str = ''):
        super().__init__(ID)
        radius = radius.to('m').magnitude
        W = duct.cross_section.width.to('m').magnitude
        H = duct.cross_section.height.to('m').magnitude
        self.R_on_W = radius / W
        self.H_on_W = H / W
        self.Pv = duct.velocity_pressure.to('Pa').magnitude
        self.interp = interpolate.interp2d(self._H_on_W, self._R_on_W, self._zeta)
        self.k_Re = ReynoldsNumberCorrectionFactor(self.R_on_W, duct.reynolds_number)

    @property
    def zeta(self) -> float:
        zeta = self.interp(self.H_on_W, self.R_on_W)
        k_Re = self.k_Re()
        return zeta[0] * k_Re

    @property
    def pressure_drop(self) -> Quantity:
        return Q_(self.zeta * self.Pv, 'Pa')


class ElbowA7G(AbstractFitting):
    description = "Elbow, rectangular, smooth radius with splitter vanes"
    _H_on_W = np.array([0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
    _R_on_W = np.array([0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50])
    _zeta_v1 = np.array([
        [0.52, 0.40, 0.43, 0.49, 0.55, 0.66, 0.75, 0.84, 0.93, 1.0, 1.1],
        [0.36, 0.27, 0.25, 0.28, 0.30, 0.35, 0.39, 0.42, 0.46, 0.49, 0.52],
        [0.28, 0.21, 0.18, 0.19, 0.20, 0.22, 0.25, 0.26, 0.28, 0.30, 0.32],
        [0.22, 0.16, 0.14, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21],
        [0.18, 0.13, 0.11, 0.11, 0.11, 0.12, 0.13, 0.14, 0.14, 0.15, 0.15],
        [0.15, 0.11, 0.09, 0.09, 0.09, 0.09, 0.10, 0.10, 0.11, 0.11, 0.12],
        [0.13, 0.09, 0.08, 0.07, 0.07, 0.08, 0.08, 0.08, 0.08, 0.09, 0.09],
        [0.11, 0.08, 0.07, 0.06, 0.06, 0.06, 0.06, 0.07, 0.07, 0.07, 0.07],
        [0.10, 0.07, 0.06, 0.05, 0.05, 0.05, 0.05, 0.05, 0.06, 0.06, 0.06],
        [0.09, 0.06, 0.05, 0.05, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05, 0.05]
    ])
    _zeta_v2 = np.array([
        [0.26, 0.20, 0.22, 0.25, 0.28, 0.33, 0.37, 0.41, 0.45, 0.48, 0.51],
        [0.17, 0.13, 0.11, 0.12, 0.13, 0.15, 0.16, 0.17, 0.19, 0.20, 0.21],
        [0.12, 0.09, 0.08, 0.08, 0.08, 0.09, 0.10, 0.10, 0.11, 0.11, 0.11],
        [0.09, 0.07, 0.06, 0.05, 0.06, 0.06, 0.06, 0.06, 0.07, 0.07, 0.07],
        [0.08, 0.05, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05, 0.05, 0.05],
        [0.06, 0.04, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04, 0.04, 0.04],
        [0.05, 0.04, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03],
        [0.05, 0.03, 0.03, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02],
        [0.04, 0.92, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02],
        [0.03, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
    ])
    _zeta_v3 = np.array([
        [0.11, 0.10, 0.12, 0.13, 0.14, 0.16, 0.18, 0.19, 0.21, 0.22, 0.23],
        [0.07, 0.05, 0.06, 0.06, 0.06, 0.07, 0.07, 0.08, 0.08, 0.08, 0.09],
        [0.05, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05],
        [0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03],
        [0.03, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02],
        [0.03, 0.02, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
        [0.02, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
        [0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
        [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
        [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
    ])
    _theta = [0, 20, 30, 45, 60, 75, 90, 110, 130, 150, 180]
    _k_theta = [0, 0.31, 0.45, 0.60, 0.78, 0.90, 1.00, 1.13, 1.20, 1.28, 1.40]

    def __init__(self, duct: Duct, radius: Quantity, no_vanes: int, theta: Quantity = Q_(90.0, 'deg'), ID: str = ''):
        super().__init__(ID)
        radius = radius.to('m').magnitude
        W = duct.cross_section.width.to('m').magnitude
        H = duct.cross_section.height.to('m').magnitude
        self.R_on_W = radius / W
        self.H_on_W = H / W
        self.Pv = duct.velocity_pressure.to('Pa').magnitude
        self.theta = theta.to('deg').magnitude
        if no_vanes == 1:
            self.interp = interpolate.interp2d(self._H_on_W, self._R_on_W, self._zeta_v1)
        elif no_vanes == 2:
            self.interp = interpolate.interp2d(self._H_on_W, self._R_on_W, self._zeta_v2)
        elif no_vanes == 3:
            self.interp = interpolate.interp2d(self._H_on_W, self._R_on_W, self._zeta_v3)
        else:
            raise ValueError("number of splitter vanes must be 1, 2 or 3")

    @property
    def zeta(self) -> float:
        zeta = self.interp(self.H_on_W, self.R_on_W)
        k_theta = np.interp(self.theta, self._theta, self._k_theta)
        return zeta[0] * k_theta

    @property
    def pressure_drop(self) -> Quantity:
        return Q_(self.zeta * self.Pv, 'Pa')


class ElbowA7H(AbstractFitting):
    description = "Elbow, rectangular, mitered with turning vanes"
    _velocity = np.array([5, 7.5, 10.0, 12.5])  # m/s
    _zeta_single_R50S38 = np.array([0.24, 0.23, 0.22, 0.20])
    _zeta_single_R114S83 = np.array([0.26, 0.24, 0.23, 0.22])
    _zeta_double_R50S38 = np.array([0.43, 0.42, 0.41, 0.40])
    _zeta_double_R50S56 = np.array([0.53, 0.53, 0.50, 0.49])
    _zeta_double_R114S83 = np.array([0.27, 0.25, 0.24, 0.23])

    def __init__(
            self,
            duct: Duct,
            vane_thickness: str = 'single',
            radius: Quantity = Q_(50, 'mm'),
            spacing: Quantity = Q_(38, 'mm'),
            ID: str = ''
    ):
        super().__init__(ID)
        self.v = duct.velocity.to('m / s').magnitude
        self.Pv = duct.velocity_pressure.to('Pa').magnitude
        if vane_thickness == 'single':
            if 50.0 <= radius('mm') < 114.0:
                self.interp = interpolate.interp1d(self._velocity, self._zeta_single_R50S38)
            elif radius('mm') >= 114.0:
                self.interp = interpolate.interp1d(self._velocity, self._zeta_single_R114S83)
            else:
                raise ValueError('vane radius must be equal to or greater than 50 mm')
        elif vane_thickness == 'double':
            if 50.0 <= radius('mm') < 114.0:
                if 38.0 <= spacing('mm') < 56.0:
                    self.interp = interpolate.interp1d(self._velocity, self._zeta_double_R50S38)
                elif spacing('mm') >= 56.0:
                    self.interp = interpolate.interp1d(self._velocity, self._zeta_double_R50S56)
                else:
                    raise ValueError('vane spacing must be equal to or greater than 38 mm')
            elif radius('mm') >= 114.0:
                self.interp = interpolate.interp1d(self._velocity, self._zeta_double_R114S83)
            else:
                raise ValueError('vane radius must be equal to or greater than 50 mm')

    @property
    def zeta(self) -> float:
        zeta = self.interp(self.v)
        return zeta

    @property
    def pressure_drop(self) -> Quantity:
        return Q_(self.zeta * self.Pv, 'Pa')


class ElbowA7I(AbstractFitting):
    description = "Elbow, 90°, rectangular, z-shaped"
    _L_on_H = np.array(
        [0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 4.0, 5.0, 6.0, 7.0, 9.0, 10.0, 1e6]
    )
    _zeta = np.array(
        [0.0, 0.62, 0.90, 1.6, 2.6, 3.6, 4.0, 4.2, 4.2, 4.2, 3.7, 3.3, 3.2, 3.1, 2.9, 2.8, 2.7, 2.6, 2.5, 2.3]
    )
    _W_on_H = np.array(
        [0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0]
    )
    _k_W_on_H = np.array(
        [1.10, 1.07, 1.04, 1.0, 0.95, 0.90, 0.83, 0.78, 0.72, 0.70]
    )

    def __init__(self, duct: Duct, height: Quantity, ID: str = ''):
        super().__init__(ID)
        self.Pv = duct.velocity_pressure.to('Pa').magnitude
        H = duct.cross_section.height.to('m').magnitude
        W = duct.cross_section.width.to('m').magnitude
        self.L_on_H = height.to('m').m / H
        self.W_on_H = W / H
        self.interp_zeta = interpolate.interp1d(self._L_on_H, self._zeta)
        self.interp_k_W_on_H = interpolate.interp1d(self._W_on_H, self._k_W_on_H)
        self.k_Re = ReynoldsNumberCorrectionFactor(R_on_W=0.0, Re=duct.reynolds_number)

    @property
    def zeta(self) -> float:
        zeta = self.interp_zeta(self.L_on_H)
        k_W_on_H = self.interp_k_W_on_H(self.W_on_H)
        k_Re = self.k_Re()
        return zeta * k_W_on_H * k_Re

    @property
    def pressure_drop(self) -> Quantity:
        return Q_(self.zeta * self.Pv, 'Pa')


class ElbowA7J(AbstractFitting):
    description = "Elbows, 90°, rectangular in different planes"
    _L_on_W = np.array([
        0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 4.0, 5.0, 6.0, 7.0, 9.0, 10.0, 1e6
    ])
    _zeta = np.array([
        1.2, 2.4, 2.9, 3.3, 3.4, 3.4, 3.4, 3.3, 3.2, 3.1, 3.2, 3.2, 3.2, 3.0, 2.9, 2.8, 2.7, 2.5, 2.4, 2.3
    ])
    _W_on_H = np.array(
        [0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0]
    )
    _k_W_on_H = np.array(
        [1.10, 1.07, 1.04, 1.0, 0.95, 0.90, 0.83, 0.78, 0.72, 0.70]
    )

    def __init__(self, duct: Duct, height: Quantity, ID: str = ''):
        super().__init__(ID)
        self.Pv = duct.velocity_pressure.to('Pa').magnitude
        W = duct.cross_section.width.to('m').magnitude
        H = duct.cross_section.height.to('m').magnitude
        self.L_on_W = height.to('m').magnitude / W
        self.W_on_H = W / H
        self.interp_zeta = interpolate.interp1d(self._L_on_W, self._zeta)
        self.interp_k_W_on_H = interpolate.interp1d(self._W_on_H, self._k_W_on_H)
        self.k_Re = ReynoldsNumberCorrectionFactor(R_on_W=0.0, Re=duct.reynolds_number)

    @property
    def zeta(self) -> float:
        zeta = self.interp_zeta(self.L_on_W)
        k_W_on_H = self.interp_k_W_on_H(self.W_on_H)
        k_Re = self.k_Re()
        return zeta * k_W_on_H * k_Re

    @property
    def pressure_drop(self) -> Quantity:
        return Q_(self.zeta * self.Pv, 'Pa')


class ElbowA7K(AbstractFitting):
    description = "Elbow, 30°, round, offset"
    _L_on_D = np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    _zeta = np.array([0.0, 0.15, 0.15, 0.16, 0.16, 0.16, 0.16])

    def __init__(self, duct: Duct, height: Quantity, ID: str = ''):
        super().__init__(ID)
        self.Pv = duct.velocity_pressure.to('Pa').magnitude
        D = duct.cross_section.equivalent_diameter.to('m').magnitude
        self.L_on_D = height.to('m').m / D
        self.interp = interpolate.interp1d(self._L_on_D, self._zeta)
        self.k_Re = ReynoldsNumberCorrectionFactor(R_on_W=0.0, Re=duct.reynolds_number)

    @property
    def zeta(self) -> float:
        zeta = self.interp(self.L_on_D)
        k_Re = self.k_Re()
        return zeta * k_Re

    @property
    def pressure_drop(self) -> Quantity:
        return Q_(self.zeta * self.Pv, 'Pa')


class DivergingTransitionA8A(AbstractFitting):
    description = "A. Transition, round, conical"
    _theta = np.array([16.0, 20.0, 30.0, 45.0, 60.0, 90.0, 120.0, 180.0])
    _A1_on_A = np.array([2.0, 4.0, 6.0, 10.0, 16])
    _zeta_Re05 = np.array([
        [0.14, 0.19, 0.32, 0.33, 0.33, 0.32, 0.31, 0.30],
        [0.23, 0.30, 0.46, 0.61, 0.68, 0.64, 0.63, 0.62],
        [0.27, 0.33, 0.48, 0.66, 0.77, 0.74, 0.73, 0.72],
        [0.29, 0.38, 0.59, 0.76, 0.80, 0.83, 0.84, 0.83],
        [0.31, 0.38, 0.60, 0.84, 0.88, 0.88, 0.88, 0.88]
    ])
    _zeta_Re2 = np.array([
        [0.07, 0.12, 0.23, 0.28, 0.27, 0.27, 0.27, 0.26],
        [0.15, 0.18, 0.36, 0.55, 0.59, 0.59, 0.58, 0.57],
        [0.19, 0.28, 0.44, 0.90, 0.70, 0.71, 0.71, 0.69],
        [0.20, 0.24, 0.43, 0.76, 0.80, 0.81, 0.81, 0.81],
        [0.21, 0.28, 0.52, 0.76, 0.87, 0.87, 0.87, 0.87]
    ])
    _zeta_Re6 = np.array([
        [0.05, 0.07, 0.12, 0.27, 0.27, 0.27, 0.27, 0.27],
        [0.17, 0.24, 0.38, 0.51, 0.56, 0.58, 0.58, 0.57],
        [0.16, 0.29, 0.46, 0.60, 0.69, 0.71, 0.70, 0.70],
        [0.21, 0.33, 0.52, 0.60, 0.76, 0.83, 0.84, 0.83],
        [0.21, 0.34, 0.56, 0.72, 0.79, 0.85, 0.87, 0.89]
    ])

    def __init__(self, duct_small: Duct, duct_large: Duct, theta: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct_small = duct_small
        self.duct_large = duct_large
        Re = duct_small.reynolds_number
        A1 = duct_large.cross_section.area.to('m ** 2').magnitude
        A = duct_small.cross_section.area.to('m ** 2').magnitude
        self.A1_on_A = A1 / A
        self.theta = theta.to('deg').magnitude
        if Re < 2e5:
            self.interp = interpolate.interp2d(self._theta, self._A1_on_A, self._zeta_Re05)
        elif 2e5 <= Re < 6e5:
            self.interp = interpolate.interp2d(self._theta, self._A1_on_A, self._zeta_Re2)
        elif Re >= 6e5:
            self.interp = interpolate.interp2d(self._theta, self._A1_on_A, self._zeta_Re6)

    @property
    def zeta_small(self) -> float:
        zeta = self.interp(self.theta, self.A1_on_A)[0]
        return zeta

    @property
    def zeta_large(self) -> float:
        zeta = convert_zeta(self.zeta_small, self.duct_small, self.duct_large)
        return zeta

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct_small.velocity_pressure.to('Pa')
        return self.zeta_small * Pv


class DivergingTransitionA8B(AbstractFitting):
    description = "Transition, rectangular, pyramidal"
    _theta = np.array([16.0, 20.0, 30.0, 45.0, 60.0, 90.0, 120.0, 180.0])
    _A1_on_A = np.array([2.0, 4.0, 6.0, 10.0])
    _zeta = np.array([
        [0.18, 0.22, 0.25, 0.29, 0.31, 0.32, 0.33, 0.30],
        [0.36, 0.43, 0.50, 0.56, 0.61, 0.63, 0.63, 0.63],
        [0.42, 0.47, 0.58, 0.68, 0.72, 0.76, 0.76, 0.75],
        [0.42, 0.49, 0.59, 0.70, 0.80, 0.87, 0.85, 0.86]
    ])

    def __init__(self, duct_small: Duct, duct_large: Duct, theta: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct_small = duct_small
        self.duct_large = duct_large
        A1 = duct_large.cross_section.area.to('m ** 2').magnitude
        A = duct_small.cross_section.area.to('m ** 2').magnitude
        self.A1_on_A = A1 / A
        self.theta = theta.to('deg').magnitude
        self.interp = interpolate.interp2d(self._theta, self._A1_on_A, self._zeta)

    @property
    def zeta_small(self) -> float:
        zeta = self.interp(self.theta, self.A1_on_A)[0]
        return zeta

    @property
    def zeta_large(self) -> float:
        zeta = convert_zeta(self.zeta_small, self.duct_small, self.duct_large)
        return zeta

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct_small.velocity_pressure.to('Pa')
        return self.zeta_small * Pv


class DivergingTransitionA8C(DivergingTransitionA8B):
    description = "Transition, round to rectangular"

    def __init__(self, duct_small: Duct, duct_large: Duct, length: Quantity, ID: str = ''):
        H = duct_large.cross_section.height.to('m').magnitude
        W = duct_large.cross_section.width.to('m').magnitude
        D = duct_small.cross_section.equivalent_diameter.to('m').magnitude
        length = length.to('m').magnitude
        theta = 2.0 * np.arctan((1.13 * (H * W) ** 0.5 - D) / (2 * length))
        super().__init__(duct_small, duct_large, Q_(theta, 'rad'), ID)


class DivergingTransitionA8D(DivergingTransitionA8B):
    description = "Transition, rectangular to round"

    def __init__(self, duct_small: Duct, duct_large: Duct, length: Quantity, ID: str = ''):
        H = duct_large.cross_section.height.to('m').magnitude
        W = duct_large.cross_section.width.to('m').magnitude
        D = duct_small.cross_section.equivalent_diameter.to('m').magnitude
        length = length.to('m').magnitude
        theta = 2.0 * np.arctan((D - 1.13 * (H * W) ** 0.5) / (2 * length))
        super().__init__(duct_small, duct_large, Q_(theta, 'rad'), ID)


class DivergingTransitionA8E(AbstractFitting):
    description = "Transition, rectangular sides straight"
    _theta = np.array([14.0, 20.0, 30.0, 45.0, 60.0, 90.0, 180.0])
    _A1_on_A = np.array([2.0, 4.0, 6.0])
    _zeta = np.array([
        [0.09, 0.12, 0.20, 0.34, 0.37, 0.38, 0.35],
        [0.16, 0.25, 0.42, 0.60, 0.68, 0.70, 0.66],
        [0.19, 0.30, 0.48, 0.65, 0.76, 0.83, 0.80]
    ])

    def __init__(self, duct_small: Duct, duct_large: Duct, theta: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct_small = duct_small
        self.duct_large = duct_large
        A1 = duct_large.cross_section.area.to('m ** 2').magnitude
        A = duct_small.cross_section.area.to('m ** 2').magnitude
        self.A1_on_A = A1 / A
        self.theta = theta.to('deg').magnitude
        self.interp = interpolate.interp2d(self._theta, self._A1_on_A, self._zeta)

    @property
    def zeta_small(self) -> float:
        zeta = self.interp(self.theta, self.A1_on_A)[0]
        return zeta

    @property
    def zeta_large(self) -> float:
        zeta = convert_zeta(self.zeta_small, self.duct_small, self.duct_large)
        return zeta

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct_small.velocity_pressure.to('Pa')
        return self.zeta_small * Pv


class DivergingTransitionA8F(AbstractFitting):
    description = "Transition, symmetric at fan with duct sides straight"
    _theta = np.array([10.0, 15.0, 20.0, 25.0, 30.0, 35.0])
    _A1_on_A = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
    _zeta = np.array([
        [0.05, 0.07, 0.09, 0.10, 0.11, 0.11],
        [0.06, 0.09, 0.11, 0.13, 0.13, 0.14],
        [0.07, 0.10, 0.13, 0.15, 0.16, 0.16],
        [0.08, 0.13, 0.16, 0.19, 0.21, 0.23],
        [0.16, 0.24, 0.29, 0.32, 0.34, 0.35],
        [0.24, 0.34, 0.39, 0.44, 0.48, 0.50]
    ])

    def __init__(self, duct_small: Duct, duct_large: Duct, theta: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct_large = duct_large
        self.duct_small = duct_small
        A = duct_small.cross_section.area.to('m ** 2').magnitude
        A1 = duct_large.cross_section.area.to('m ** 2').magnitude
        self.A1_on_A = A1 / A
        self.theta = theta.to('deg').magnitude
        self.interp = interpolate.interp2d(self._A1_on_A, self._theta, self._zeta)

    @property
    def zeta_small(self) -> float:
        zeta = self.interp(self.A1_on_A, self.theta)[0]
        return zeta

    @property
    def zeta_large(self) -> float:
        zeta = convert_zeta(self.zeta_small, self.duct_small, self.duct_large)
        return zeta

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct_small.velocity_pressure.to('Pa')
        return self.zeta_small * Pv


class DivergingTransitionA8G(DivergingTransitionA8F):
    description = "Transition, asymmetric at fan with duct sides straight, top level"
    _zeta = np.array([
        [0.08, 0.09, 0.10, 0.10, 0.11, 0.11],
        [0.10, 0.11, 0.12, 0.13, 0.14, 0.15],
        [0.12, 0.14, 0.15, 0.16, 0.17, 0.18],
        [0.15, 0.18, 0.21, 0.23, 0.25, 0.26],
        [0.18, 0.25, 0.30, 0.33, 0.35, 0.35],
        [0.21, 0.31, 0.38, 0.41, 0.43, 0.44]
    ])


class DivergingTransitionA8H(DivergingTransitionA8F):
    description = "Transition, asymmetric at fan with duct sides straight, top 10° down"
    _zeta = np.array([
        [0.11, 0.13, 0.14, 0.14, 0.14, 0.14],
        [0.13, 0.15, 0.16, 0.17, 0.18, 0.18],
        [0.19, 0.22, 0.24, 0.26, 0.28, 0.30],
        [0.29, 0.32, 0.35, 0.37, 0.39, 0.40],
        [0.36, 0.42, 0.46, 0.49, 0.51, 0.51],
        [0.44, 0.54, 0.61, 0.64, 0.66, 0.66]
    ])


class DivergingTransitionA8I(DivergingTransitionA8F):
    description = "Transition, asymmetric at fan with duct sides straight, top 10° up"
    _zeta = np.array([
        [0.05, 0.08, 0.11, 0.13, 0.13, 0.14],
        [0.06, 0.10, 0.12, 0.14, 0.15, 0.15],
        [0.07, 0.11, 0.14, 0.15, 0.16, 0.16],
        [0.09, 0.14, 0.18, 0.20, 0.21, 0.22],
        [0.13, 0.18, 0.23, 0.26, 0.28, 0.29],
        [0.15, 0.23, 0.28, 0.33, 0.35, 0.36]
    ])


class DivergingTransitionA8J(DivergingTransitionA8F):
    description = "Transition, pyramidal at fan with duct"
    _theta = np.array([10.0, 15.0, 20.0, 25.0, 30.0])
    _zeta = np.array([
        [0.10, 0.18, 0.21, 0.23, 0.24, 0.25],
        [0.23, 0.33, 0.38, 0.40, 0.42, 0.44],
        [0.31, 0.43, 0.48, 0.53, 0.56, 0.58],
        [0.36, 0.49, 0.55, 0.58, 0.62, 0.64],
        [0.42, 0.53, 0.59, 0.64, 0.67, 0.69]
    ])


class ConvergingTransitionA9A(AbstractFitting):
    description = "Contraction, round and rectangular, gradual to abrupt"
    _theta = np.array([10, 15, 20, 30, 40, 50, 60, 90, 120, 150, 180])
    _A1_on_A = np.array([2, 4, 6, 10])
    _zeta = np.array([
        [0.05, 0.05, 0.05, 0.05, 0.05, 0.06, 0.06, 0.12, 0.18, 0.24, 0.26],
        [0.05, 0.04, 0.04, 0.04, 0.04, 0.07, 0.07, 0.17, 0.27, 0.35, 0.41],
        [0.05, 0.04, 0.04, 0.04, 0.04, 0.07, 0.07, 0.18, 0.28, 0.36, 0.42],
        [0.05, 0.05, 0.05, 0.05, 0.05, 0.08, 0.08, 0.19, 0.29, 0.37, 0.43]
    ])

    def __init__(self, duct_small: Duct, duct_large: Duct, theta: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct_small = duct_small
        self.duct_large = duct_large
        A1 = duct_large.cross_section.area.to('m ** 2').magnitude
        A = duct_small.cross_section.area.to('m ** 2').magnitude
        self.A1_on_A = A1 / A
        self.theta = theta.to('deg').magnitude
        self.interp = interpolate.interp2d(self._theta, self._A1_on_A, self._zeta)

    @property
    def zeta_small(self) -> float:
        zeta = self.interp(self.theta, self.A1_on_A)[0]
        return zeta

    @property
    def zeta_large(self) -> float:
        zeta = convert_zeta(self.zeta_small, self.duct_small, self.duct_large)
        return zeta

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct_small.velocity_pressure.to('Pa')
        return self.zeta_small * Pv


class ConvergingTransitionA9B(AbstractFitting):
    description = "Contraction, conical, round and rectangular"
    _theta = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 60.0, 100.0, 140.0, 180.0])
    _A_on_A1 = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0])
    _k_A = np.array([1.0, 0.85, 0.68, 0.50, 0.30, 0.18, 0.0])
    _L_on_D = np.array([0.025, 0.05, 0.075, 0.1, 0.15, 0.60])
    _zeta = np.array([
        [0.50, 0.47, 0.45, 0.43, 0.41, 0.40, 0.42, 0.45, 0.50],
        [0.50, 0.45, 0.41, 0.36, 0.33, 0.30, 0.35, 0.42, 0.50],
        [0.50, 0.42, 0.35, 0.30, 0.26, 0.23, 0.30, 0.40, 0.50],
        [0.50, 0.39, 0.32, 0.25, 0.22, 0.18, 0.27, 0.38, 0.50],
        [0.50, 0.37, 0.27, 0.20, 0.16, 0.15, 0.25, 0.37, 0.50],
        [0.50, 0.27, 0.18, 0.13, 0.11, 0.12, 0.23, 0.36, 0.50]
    ])

    def __init__(self, duct_small: Duct, duct_large: Duct, theta: Quantity, length: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct_small = duct_small
        self.duct_large = duct_large
        A1 = duct_large.cross_section.area.to('m ** 2').magnitude
        A = duct_small.cross_section.area.to('m ** 2').magnitude
        D = duct_small.cross_section.hydraulic_diameter.to('m').magnitude
        self.A_on_A1 = A / A1
        self.theta = theta.to('deg').magnitude
        self.L_on_D = length.to('m').magnitude / D
        self.interp_k_A = interpolate.interp1d(self._A_on_A1, self._k_A)
        self.interp_zeta = interpolate.interp2d(self._theta, self._L_on_D, self._zeta)

    @property
    def zeta_small(self) -> float:
        zeta = self.interp_zeta(self.theta, self.L_on_D)[0]
        k_A = self.interp_k_A(self.A_on_A1)
        return zeta * k_A

    @property
    def zeta_large(self) -> float:
        zeta = convert_zeta(self.zeta_small, self.duct_small, self.duct_large)
        return zeta

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct_small.velocity_pressure.to('Pa')
        return self.zeta_small * Pv


class ConvergingTransitionA9C(AbstractFitting):
    description = "Contraction, rectangular slot to round"
    _Re = np.array([1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0, 40.0])
    _zeta = np.array([0.27, 0.25, 0.20, 0.17, 0.14, 0.11, 0.04, 0])

    def __init__(self, duct_small: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_down = duct_small
        self.Re = duct_small.reynolds_number * 1e-4
        self.interp = interpolate.interp1d(self._Re, self._zeta)

    @property
    def zeta_small(self) -> float:
        zeta = self.interp(self.Re)
        return zeta

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct_down.velocity_pressure.to('Pa')
        return self.zeta_small * Pv


class ConvergingJunctionA10A(AbstractFitting):
    description = "Converging wye 45°, round"
    _Ab_on_Ac = np.array([0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0])
    _Vb_on_Vc = np.array([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0])
    _Vs_on_Vc = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    _zeta_b = np.array([
        [-0.56, -0.44, -0.35, -0.28, -0.15, -0.04, 0.05],
        [-0.48, -0.37, -0.28, -0.21, -0.09, 0.02, 0.11],
        [-0.38, -0.27, -0.19, -0.12, 0.0, 0.10, 0.18],
        [-0.26, -0.16, -0.08, -0.01, 0.10, 0.20, 0.28],
        [-0.21, -0.02, 0.05, 0.12, 0.23, 0.32, 0.40],
        [0.04, 0.13, 0.21, 0.27, 0.37, 0.46, 0.53],
        [0.22, 0.31, 0.38, 0.44, 0.53, 0.62, 0.69],
        [1.4, 1.5, 1.5, 1.6, 1.7, 1.7, 1.8],
        [3.1, 3.2, 3.2, 3.2, 3.3, 3.3, 3.3],
        [5.3, 5.3, 5.3, 5.4, 5.4, 5.4, 5.4],
        [8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0]
    ])
    _zeta_c = np.array([
        [-8.6, -4.1, -2.5, -1.7, -0.97, -0.58, -0.34],
        [-6.7, -3.1, -1.9, -1.3, -0.67, -0.36, -0.18],
        [-5.0, -2.2, -1.3, -0.88, -0.42, -0.19, -0.05],
        [-3.5, -1.5, -0.88, -0.55, -0.21, -0.05, 0.05],
        [-2.3, -0.95, -0.51, -0.28, -0.06, 0.06, 0.13],
        [-1.3, -0.50, -0.22, -0.09, 0.05, 0.12, 0.17],
        [-0.63, -0.18, -0.03, 0.04, 0.12, 0.16, 0.18],
        [-0.18, 0.01, 0.07, 0.10, 0.13, 0.15, 0.17],
        [0.03, 0.07, 0.08, 0.09, 0.10, 0.11, 0.13],
        [-0.01, 0.0, 0.0, 0.10, 0.02, 0.04, 0.05]
    ])

    def __init__(self, duct_s: Duct, duct_b: Duct, duct_c: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        Vs = duct_s.volume_flow_rate.to('m ** 3 / s').magnitude
        Vb = duct_b.volume_flow_rate.to('m ** 3 / s').magnitude
        Vc = duct_c.volume_flow_rate.to('m ** 3 / s').magnitude
        Ab = duct_b.cross_section.area.to('m ** 2').magnitude
        Ac = duct_c.cross_section.area.to('m ** 2').magnitude
        self.Vb_on_Vc = Vb / Vc
        self.Vs_on_Vc = Vs / Vc
        self.Ab_on_Ac = Ab / Ac
        self.interp_zeta_b = interpolate.interp2d(self._Ab_on_Ac, self._Vb_on_Vc, self._zeta_b)
        self.interp_zeta_c = interpolate.interp2d(self._Ab_on_Ac, self._Vs_on_Vc, self._zeta_c)

    @property
    def zeta_b(self) -> float:
        zeta = self.interp_zeta_b(self.Ab_on_Ac, self.Vb_on_Vc)[0]
        zeta = convert_zeta(zeta, self.duct_c, self.duct_b)   # refer `zeta_b` to velocity in branch duct
        return zeta

    @property
    def zeta_c(self) -> float:
        zeta = self.interp_zeta_c(self.Ab_on_Ac, self.Vs_on_Vc)[0]
        return zeta

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class ConvergingJunctionA10B(AbstractFitting):
    description = "Converging tee, 90°, round"
    _Ab_on_Ac = np.array([0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0])
    _Vb_on_Vc = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    _zeta_b = np.array([
        [0.40, -0.37, -0.51, -0.46, -0.50, -0.51, -0.52],
        [3.8, 0.72, 0.17, -0.02, -0.14, -0.18, -0.24],
        [9.2, 2.3, 1.0, 0.44, 0.21, 0.11, -0.08],
        [16.0, 4.3, 2.1, 0.94, 0.54, 0.40, 0.32],
        [26.0, 6.8, 3.2, 1.1, 0.66, 0.49, 0.42],
        [37.0, 9.7, 4.7, 1.6, 0.92, 0.69, 0.57],
        [43.0, 13.0, 6.3, 2.1, 1.2, 0.88, 0.72],
        [65.0, 17.0, 7.9, 2.7, 1.5, 1.1, 0.86],
        [82.0, 21.0, 9.7, 3.4, 1.8, 1.2, 0.99],
        [101.0, 26.0, 12.0, 4.0, 2.1, 1.4, 1.1]
    ])
    _zeta_c = np.array(
        [0.16, 0.27, 0.38, 0.46, 0.53, 0.57, 0.59, 0.60, 0.59, 0.55]
    )

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        Ab = duct_b.cross_section.area.to('m ** 2').magnitude
        Ac = duct_c.cross_section.area.to('m ** 2').magnitude
        Vb = duct_b.volume_flow_rate.to('m ** 3 / s').magnitude
        Vc = duct_c.volume_flow_rate.to('m ** 3 / s').magnitude
        self.Ab_on_Ac = Ab / Ac
        self.Vb_on_Vc = Vb / Vc
        self.interp_zeta_b = interpolate.interp2d(self._Ab_on_Ac, self._Vb_on_Vc, self._zeta_b)
        self.interp_zeta_c = interpolate.interp1d(
            self._Vb_on_Vc, self._zeta_c, bounds_error=False, fill_value='extrapolate'
        )

    @property
    def zeta_b(self) -> float:
        zeta = self.interp_zeta_b(self.Ab_on_Ac, self.Vb_on_Vc)[0]
        zeta = convert_zeta(zeta, self.duct_c, self.duct_b)  # refer `zeta_b` to velocity in branch duct
        return zeta

    @property
    def zeta_c(self) -> float:
        zeta = self.interp_zeta_c(self.Vb_on_Vc)
        return zeta

    @property
    def zeta_s(self) -> float:
        # refer `zeta_c` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class ConvergingJunctionA10C(AbstractFitting):
    description = "Converging tee, round branch to rectangular main"
    _Vb_on_Vc = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    _zeta_b = np.array([
        [-0.63, -0.55, -0.13, 0.23, 0.78, 1.30, 1.93, 3.10, 4.88, 5.60],
        [-0.49, -0.21, -0.23, 0.60, 1.27, 2.06, 2.75, 3.70, 4.93, 5.95]
    ])
    _zeta_c = np.array(
        [0.16, 0.27, 0.38, 0.46, 0.53, 0.57, 0.59, 0.60, 0.59, 0.55]
    )

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        Vb = duct_b.volume_flow_rate.to('m ** 3 / s').magnitude
        Vc = duct_c.volume_flow_rate.to('m ** 3 / s').magnitude
        vc = duct_c.velocity.to('ft / min').magnitude
        self.Vb_on_Vc = Vb / Vc
        self.interp_zeta_c = interpolate.interp1d(self._Vb_on_Vc, self._zeta_c)
        if vc < 1200:
            zeta_b = self._zeta_b[0]
        else:
            zeta_b = self._zeta_b[1]
        self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, zeta_b)

    @property
    def zeta_b(self) -> float:
        zeta = self.interp_zeta_b(self.Vb_on_Vc)
        zeta = convert_zeta(zeta, self.duct_c, self.duct_b)  # refer `zeta_b` to velocity in branch duct
        return zeta

    @property
    def zeta_c(self) -> float:
        zeta = self.interp_zeta_c(self.Vb_on_Vc)
        return zeta

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class ConvergingJunctionA10D(AbstractFitting):
    description = "Converging tee, rectangular main and branch"
    _Vb_on_Vc = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    _zeta_b = np.array([
        [-0.75, -0.53, -0.03, 0.33, 1.03, 1.10, 2.15, 2.93, 4.18, 4.78],
        [-0.69, -0.21, -0.23, 0.67, 1.17, 1.66, 2.67, 3.36, 3.93, 5.13]
    ])
    _zeta_c = np.array(
        [0.16, 0.27, 0.38, 0.46, 0.53, 0.57, 0.59, 0.60, 0.59, 0.55]
    )

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        Vb = duct_b.volume_flow_rate.to('m ** 3 / s').magnitude
        Vc = duct_c.volume_flow_rate.to('m ** 3 / s').magnitude
        vc = duct_c.velocity.to('ft / min').magnitude
        self.Vb_on_Vc = Vb / Vc
        self.interp_zeta_c = interpolate.interp1d(self._Vb_on_Vc, self._zeta_c)
        if vc < 1200:
            zeta_b = self._zeta_b[0]
        else:
            zeta_b = self._zeta_b[1]
        self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, zeta_b)

    @property
    def zeta_b(self) -> float:
        zeta = self.interp_zeta_b(self.Vb_on_Vc)
        zeta = convert_zeta(zeta, self.duct_c, self.duct_b)  # refer `zeta_b` to velocity in branch duct
        return zeta

    @property
    def zeta_c(self) -> float:
        zeta = self.interp_zeta_c(self.Vb_on_Vc)
        return zeta

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class ConvergingJunctionA10E(AbstractFitting):
    description = "Converging wye 45°, conical round"
    _Vb_on_Vs = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    _zeta_b = np.array([
        [-2.4, -0.01, 2.0, 3.8, 5.3, 6.6, 7.8, 8.9, 9.8, 11],
        [-2.8, -1.2, 0.12, 1.1, 1.9, 2.6, 3.2, 3.7, 4.2, 4.6],
        [-1.2, 0.93, 2.8, 4.5, 5.9, 7.2, 8.4, 9.5, 10, 11],
        [-1.6, -0.27, 0.81, 1.7, 2.4, 3.0, 3.6, 4.1, 4.5, 4.9],
        [-1.8, -0.72, 0.07, 0.66, 1.1, 1.5, 1.8, 2.1, 2.3, 2.5],
        [-0.046, 1.5, 3.3, 4.9, 6.4, 7.7, 8.8, 9.9, 11, 12],
        [-0.094, 0.25, 1.2, 2.0, 2.7, 3.3, 3.8, 4.2, 4.7, 5.0],
        [-1.1, -0.24, 0.42, 0.92, 1.3, 1.6, 1.9, 2.1, 2.3, 2.5],
        [-1.2, -0.38, 0.18, 0.58, 0.88, 1.1, 1.3, 1.5, 1.6, 1.7],
        [-0.55, 1.3, 3.1, 4.7, 6.1, 7.4, 8.6, 9.6, 11, 12],
        [-1.1, 0, 0.88, 1.6, 2.3, 2.8, 3.3, 3.7, 4.1, 4.5],
        [-1.2, -0.48, 0.1, 0.54, 0.89, 1.2, 1.4, 1.6, 1.8, 2.0],
        [-1.3, -0.62, -0.14, 0.21, 0.47, 0.68, 0.85, 0.99, 1.1, 1.2],
        [-1.3, -0.69, -0.26, 0.04, 0.26, 0.42, 0.57, 0.66, 0.75, 0.82],
        [0.06, 1.8, 3.5, 5.1, 6.5, 7.8, 8.9, 10, 11, 12],
        [-0.52, 0.35, 1.1, 1.7, 2.3, 2.8, 3.2, 3.6, 3.9, 4.2],
        [-0.67, -0.05, 0.43, 0.8, 1.1, 1.4, 1.6, 1.8, 1.9, 2.1],
        [-0.73, -0.19, 0.18, 0.46, 0.68, 0.85, 0.99, 1.1, 1.2, 1.3],
        [-0.75, -0.27, 0.05, 0.28, 0.45, 0.58, 0.68, 0.76, 0.83, 0.88],
        [-0.77, -0.31, -0.02, 0.18, 0.32, 0.43, 0.50, 0.56, 0.61, 0.65],
        [-0.78, -0.34, -0.07, 0.12, 0.24, 0.33, 0.39, 0.44, 0.47, 0.50],
        [9999.9, 2.1, 3.7, 5.2, 6.6, 7.8, 9.0, 11, 11, 12],
        [9999.9, 0.54, 1.2, 1.8, 2.3, 2.7, 3.1, 3.7, 3.7, 4.0],
        [9999.9, 0.21, 0.62, 0.96, 1.2, 1.5, 1.7, 2.0, 2.0, 2.1],
        [9999.9, 0.05, 0.37, 0.60, 0.79, 0.93, 1.1, 1.2, 1.2, 1.3],
        [9999.9, -0.02, 0.23, 0.42, 0.55, 0.66, 0.73, 0.80, 0.85, 0.89],
        [9999.9, -0.1, 0.11, 0.24, 0.33, 0.39, 0.43, 0.46, 0.47, 0.48],
        [9999.9, -0.14, 0.05, 0.16, 0.23, 0.27, 0.29, 0.30, 0.30, 0.29]
    ])

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        As = duct_s.cross_section.area.to('m ** 2').magnitude
        Ac = duct_c.cross_section.area.to('m ** 2').magnitude
        Ab = duct_b.cross_section.area.to('m ** 2').magnitude
        As_on_Ac = round(As / Ac, 1)
        self.Ab_on_Ac = Ab / Ac
        Vb = duct_b.volume_flow_rate.to('m ** 3 / s').magnitude
        Vs = duct_s.volume_flow_rate.to('m ** 3 / s').magnitude
        self.Vb_on_Vs = Vb / Vs
        if As_on_Ac == 0.3:
            _Ab_on_Ac = np.array([0.2, 0.3])
            _zeta_b = self._zeta_b[0:2]
            self.interp_zeta_b = interpolate.interp2d(_Ab_on_Ac, self._Vb_on_Vs, _zeta_b)
        elif As_on_Ac == 0.4:
            _Ab_on_Ac = np.array([0.2, 0.3, 0.4])
            _zeta_b = self._zeta_b[2:5]
            self.interp_zeta_b = interpolate.interp2d(_Ab_on_Ac, self._Vb_on_Vs, _zeta_b)
        elif As_on_Ac == 0.5:
            _Ab_on_Ac = np.array([0.2, 0.3, 0.4, 0.5])
            _zeta_b = self._zeta_b[5:9]
            self.interp_zeta_b = interpolate.interp2d(_Ab_on_Ac, self._Vb_on_Vs, _zeta_b)
        elif As_on_Ac == 0.6:
            _Ab_on_Ac = np.array([0.2, 0.3, 0.4, 0.5, 0.6])
            _zeta_b = self._zeta_b[9:14]
            self.interp_zeta_b = interpolate.interp2d(_Ab_on_Ac, self._Vb_on_Vs, _zeta_b)
        elif As_on_Ac == 0.8:
            _Ab_on_Ac = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
            _zeta_b = self._zeta_b[14:21]
            self.interp_zeta_b = interpolate.interp2d(_Ab_on_Ac, self._Vb_on_Vs, _zeta_b)
        elif As_on_Ac == 1.0:
            _Ab_on_Ac = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0])
            _zeta_b = self._zeta_b[21:]
            self.interp_zeta_b = interpolate.interp2d(_Ab_on_Ac, self._Vb_on_Vs, _zeta_b)

    @property
    def zeta_b(self) -> float:
        zeta = self.interp_zeta_b(self.Ab_on_Ac, self.Vb_on_Vs)[0]
        zeta = convert_zeta(zeta, self.duct_c, self.duct_b)  # refer `zeta_b` to velocity in branch duct
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv


class ConvergingJunctionA10F(AbstractFitting):
    description = "Converging tee, 45° entry branch to rectangular main"
    _Vb_on_Vc = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    _Vb_on_Vs = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    _zeta_b = np.array([
        [-0.83, -0.68, -0.30, 0.28, 0.55, 1.03, 1.50, 1.93, 2.50, 3.03],
        [-0.72, -0.52, -0.23, 0.34, 0.76, 1.14, 1.83, 2.01, 2.90, 3.63]
    ])
    _zeta_c = np.array([
        [5.3, -0.01, 2.0, 1.1, 0.34, -0.20, -0.61, -0.93, -1.2, -1.4],
        [5.4, 3.7, 2.5, 1.6, 1.0, 0.53, 0.16, -0.14, -0.38, -0.58],
        [1.9, 1.1, 0.46, -0.07, -0.49, -0.83, -1.1, -1.3, -1.5, -1.7],
        [2.0, 1.4, 0.81, 0.42, 0.08, -0.20, -0.43, -0.62, -0.78, -0.92],
        [2.0, 1.5, 1.0, 0.68, 0.39, 0.16, -0.04, -0.21, -0.35, -0.47],
        [0.77, 0.34, -0.09, -0.48, -0.81, -1.1, -1.3, -1.5, -1.7, -1.8],
        [0.85, 0.56, 0.25, -0.03, -0.27, -0.48, -0.67, -0.82, -0.96, -1.1],
        [0.88, 0.66, 0.43, 0.21, 0.02, -0.15, -0.30, -0.42, -0.54, -0.64],
        [0.91, 0.73, 0.54, 0.36, 0.21, 0.06, -0.06, -0.17, -0.26, -0.35],
        [0.3, 0, -0.34, -0.67, -0.96, -1.2, -1.4, -1.6, -1.8, -1.9],
        [0.37, 0.21, -0.02, -0.24, -0.44, -0.63, -0.79, -0.93, -1.1, -1.2],
        [0.40, 0.31, 0.16, -0.01, -0.16, -0.30, -0.43, -0.54, -0.64, -0.73],
        [0.43, 0.37, 0.26, 0.14, 0.02, -0.09, -0.20, -0.29, -0.37, -0.45],
        [0.44, 0.41, 0.33, 0.24, 0.14, 0.05, -0.03, -0.11, -0.18, -0.25],
        [-0.06, -0.27, -0.57, -0.86, -1.1, -1.4, -1.6, -1.7, -1.9, -2.0],
        [0, -0.08, -0.25, -0.43, -0.62, -0.78, -0.93, -1.1, -1.2, -1.3],
        [0.04, 0.02, -0.08, -0.21, -0.34, -0.46, -0.57, -0.67, -0.77, -0.85],
        [0.06, 0.08, 0.02, -0.06, -0.16, -0.25, -0.34, -0.42, -0.50, -0.57],
        [0.07, 0.12, 0.09, 0.03, -0.04, -0.11, -0.18, -0.25, -0.31, -0.37],
        [0.08, 0.15, 0.14, 0.10, 0.05, -0.01, -0.07, -0.12, -0.17, -0.22],
        [0.09, 0.17, 0.18, 0.16, 0.11, 0.07, 0.02, -0.02, -0.07, -0.11],
        [9999.9, -0.39, -0.67, -0.96, -1.2, -1.5, -1.6, -1.8, -2.0, -2.1],
        [9999.9, -0.19, -0.35, -0.54, -0.71, -0.87, -1.0, -1.2, -1.3, -1.4],
        [9999.9, -0.10, -0.19, -0.31, -0.43, -0.55, -0.66, -0.77, -0.86, -0.94],
        [9999.9, -0.04, -0.09, -0.17, -0.26, -0.35, -0.44, -0.52, -0.59, -0.66],
        [9999.9, 0, -0.02, -0.07, -0.14, -0.21, -0.28, -0.34, -0.40, -0.46],
        [9999.9, 0.06, 0.07, 0.05, 0.02, -0.03, -0.07, -0.12, -0.16, -0.20],
        [9999.9, 0.09, 0.13, 0.13, 0.11, 0.08, 0.06, 0.03, -0.01, -0.03]
    ])
        
    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        Vb = duct_b.volume_flow_rate.to('m ** 3 / s').magnitude
        Vc = duct_c.volume_flow_rate.to('m ** 3 / s').magnitude
        Vs = duct_s.volume_flow_rate.to('m ** 3 / s').magnitude
        vc = duct_c.velocity.to('ft / min').magnitude
        self.Vb_on_Vc = Vb / Vc
        self.Vb_on_Vs = Vb / Vs
        As = duct_s.cross_section.area.to('m ** 2').magnitude
        Ac = duct_c.cross_section.area.to('m ** 2').magnitude
        Ab = duct_b.cross_section.area.to('m ** 2').magnitude
        As_on_Ac = round(As / Ac, 1)
        self.Ab_on_Ac = Ab / Ac
        if vc < 1200:
            zeta_b = self._zeta_b[0]
        else:
            zeta_b = self._zeta_b[1]
        self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, zeta_b)
        if As_on_Ac == 0.3:
            _Ab_on_Ac = np.array([0.2, 0.3])
            _zeta_c = self._zeta_c[0:2]
            self.interp_zeta_c = interpolate.interp2d(_Ab_on_Ac, self._Vb_on_Vs, _zeta_c)
        elif As_on_Ac == 0.4:
            _Ab_on_Ac = np.array([0.2, 0.3, 0.4])
            _zeta_c = self._zeta_c[2:5]
            self.interp_zeta_c = interpolate.interp2d(_Ab_on_Ac, self._Vb_on_Vs, _zeta_c)
        elif As_on_Ac == 0.5:
            _Ab_on_Ac = np.array([0.2, 0.3, 0.4, 0.5])
            _zeta_c = self._zeta_c[5:9]
            self.interp_zeta_c = interpolate.interp2d(_Ab_on_Ac, self._Vb_on_Vs, _zeta_c)
        elif As_on_Ac == 0.6:
            _Ab_on_Ac = np.array([0.2, 0.3, 0.4, 0.5, 0.6])
            _zeta_c = self._zeta_c[9:14]
            self.interp_zeta_c = interpolate.interp2d(_Ab_on_Ac, self._Vb_on_Vs, _zeta_c)
        elif As_on_Ac == 0.8:
            _Ab_on_Ac = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
            _zeta_c = self._zeta_c[14:21]
            self.interp_zeta_c = interpolate.interp2d(_Ab_on_Ac, self._Vb_on_Vs, _zeta_c)
        elif As_on_Ac == 1.0:
            _Ab_on_Ac = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0])
            _zeta_c = self._zeta_c[21:]
            self.interp_zeta_c = interpolate.interp2d(_Ab_on_Ac, self._Vb_on_Vs, _zeta_c)
        else:
            self.interp_zeta_c = None

    @property
    def zeta_b(self) -> float:
        if self.interp_zeta_b:
            zeta = self.interp_zeta_b(self.Vb_on_Vc)
            zeta = convert_zeta(zeta, self.duct_c, self.duct_b)  # refer `zeta_b` to velocity in branch duct
            return zeta
        else:
            return float('nan')

    @property
    def zeta_c(self) -> float:
        if self.interp_zeta_c:
            zeta = self.interp_zeta_c(self.Ab_on_Ac, self.Vb_on_Vs)[0]
            return zeta
        else:
            return float('nan')

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class ConvergingJunctionA10G(AbstractFitting):
    description = "Symmetrical wye, dovetail, rectangular"

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        Ab = duct_b.cross_section.area.to('m ** 2').magnitude
        Ac = duct_c.cross_section.area.to('m ** 2').magnitude
        Ab_on_Ac = round(Ab / Ac, 1)
        if Ab_on_Ac == 0.5:
            self._zeta = 0.23
        elif Ab_on_Ac == 1.0:
            self._zeta = 0.07
        else:
            self._zeta = float('nan')

    @property
    def zeta_b(self) -> float:
        zeta = convert_zeta(self._zeta, self.duct_c, self.duct_b)
        return zeta

    @property
    def zeta_c(self) -> float:
        return self._zeta

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class ConvergingJunctionA10H(AbstractFitting):
    description = "Converging wye, rectangular"
    _Vb_on_Vc = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    _zeta_b = np.array([
        [-0.50, 0, 0.50, 1.2, 2.2, 3.7, 5.8, 8.4, 11],
        [-1.2, -0.40, 0.40, 1.6, 3.0, 4.8, 6.8, 8.9, 11],
        [-0.50, -0.20, 0, 0.25, 0.45, 0.70, 1.0, 1.5, 2.0],
        [-1.0, -0.60, -0.20, 0.10, 0.30, 0.60, 1.0, 1.5, 2.0],
        [-2.2, -1.5, -0.95, -0.50, 0, 0.40, 0.80, 1.3, 1.9],
        [-0.60, -0.30, -0.10, -0.04, 0.13, 0.21, 0.29, 0.36, 0.42],
        [-1.2, -0.80, -0.40, -0.20, 0, 0.16, 0.24, 0.32, 0.38],
        [-2.1, -1.4, -0.90, -0.5, -0.20, 0, 0.20, 0.25, 0.30]
    ])
    _zeta_m = np.array([
        [0.30, 0.30, 0.20, -0.1, -0.45, -0.92, -1.5, -0.20, -0.26],
        [0.17, 0.16, 0.10, 0, -0.08, -0.18, -0.27, -0.37, -0.46],
        [0.27, 0.35, 0.32, 0.25, 0.12, -0.03, -0.23, -0.42, -0.58],
        [1.2, 1.1, 0.90, 0.65, 0.35, 0, -0.40, -0.80, -1.3],
        [0.18, 0.24, 0.27, 0.26, 0.23, 0.18, 0.10, 0, -0.12],
        [0.75, 0.36, 0.38, 0.35, 0.27, 0.18, 0.05, -0.08, -0.22],
        [0.80, 0.87, 0.80, 0.68, 0.55, 0.40, 0.25, 0.08, -0.10]
    ])

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        Vb = duct_b.volume_flow_rate.to('m ** 3 / s').magnitude
        Vc = duct_c.volume_flow_rate.to('m ** 3 / s').magnitude
        self.Vb_on_Vc = Vb / Vc
        Ab = duct_b.cross_section.area.to('m ** 2').magnitude
        Ac = duct_c.cross_section.area.to('m ** 2').magnitude
        As = duct_s.cross_section.area.to('m ** 2').magnitude
        Ab_on_As = round(Ab / As, 2)
        Ab_on_Ac = round(Ab / Ac, 2)
        if Ab_on_As == 0.25 and Ab_on_Ac == 0.25:
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[0])
        elif Ab_on_As == 0.33 and Ab_on_Ac == 0.25:
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[1])
        elif Ab_on_As == 0.50 and Ab_on_Ac == 0.50:
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[2])
        elif Ab_on_As == 0.67 and Ab_on_Ac == 0.50:
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[3])
        elif Ab_on_As == 1.00 and Ab_on_Ac == 0.50:
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[4])
        elif Ab_on_As == 1.00 and Ab_on_Ac == 1.00:
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[5])
        elif Ab_on_As == 1.33 and Ab_on_Ac == 1.00:
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[6])
        elif Ab_on_As == 2.00 and Ab_on_Ac == 1.00:
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[7])
        else:
            self.interp_zeta_b = None
        As_on_Ac = round(As / Ac, 2)
        if As_on_Ac == 0.75 and Ab_on_Ac == 0.25:
            self.interp_zeta_m = interpolate.interp1d(self._Vb_on_Vc, self._zeta_m[0])
        elif As_on_Ac == 1.05 and Ab_on_Ac == 0.50:
            self.interp_zeta_m = interpolate.interp1d(self._Vb_on_Vc, self._zeta_m[1])
        elif As_on_Ac == 0.75 and Ab_on_Ac == 0.50:
            self.interp_zeta_m = interpolate.interp1d(self._Vb_on_Vc, self._zeta_m[2])
        elif As_on_Ac == 0.50 and Ab_on_Ac == 0.50:
            self.interp_zeta_m = interpolate.interp1d(self._Vb_on_Vc, self._zeta_m[3])
        elif As_on_Ac == 1.00 and Ab_on_Ac == 1.00:
            self.interp_zeta_m = interpolate.interp1d(self._Vb_on_Vc, self._zeta_m[4])
        elif As_on_Ac == 0.75 and Ab_on_Ac == 1.00:
            self.interp_zeta_m = interpolate.interp1d(self._Vb_on_Vc, self._zeta_m[5])
        elif As_on_Ac == 0.50 and Ab_on_Ac == 1.00:
            self.interp_zeta_m = interpolate.interp1d(self._Vb_on_Vc, self._zeta_m[6])
        else:
            self.interp_zeta_m = None

    @property
    def zeta_b(self) -> float:
        if self.interp_zeta_b:
            zeta = self.interp_zeta_b(self.Vb_on_Vc)
            zeta = convert_zeta(zeta, self.duct_c, self.duct_b)  # refer zeta to branch duct
            return zeta
        else:
            return float('nan')

    @property
    def zeta_c(self) -> float:
        if self.interp_zeta_m:
            zeta = self.interp_zeta_m(self.Vb_on_Vc)
            return zeta
        else:
            return float('nan')

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class ConvergingJunctionA10I(AbstractFitting):
    description = "Wye, rectangular and round"
    _theta = np.array([15.0, 30.0, 45.0])
    _Vb_on_Vc = np.array([0.0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00])
    _zeta = np.array([
        [-2.6, -0.19, -0.13, -0.77, -0.30, 0.10, 0.41, 0.67, 0.85, 0.97, 1.0],
        [-2.1, -0.15, -0.10, -0.53, -0.1, 0.28, 0.69, 0.91, 1.1, 1.4, 1.6],
        [-1.3, -0.93, -0.55, -0.16, 0.20, 0.56, 0.92, 1.26, 1.60, 2.00, 2.3]
    ])

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, theta: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        Vb = duct_b.volume_flow_rate.to('m ** 3 / s').magnitude
        Vc = duct_c.volume_flow_rate.to('m ** 3 / s').magnitude
        self.Vb_on_Vc = Vb / Vc
        self.theta = theta.to('deg').magnitude
        self.interp_zeta = interpolate.interp2d(self._theta, self._Vb_on_Vc, self._zeta)

    @property
    def zeta_b(self) -> float:
        zeta = self.interp_zeta(self.theta, self.Vb_on_Vc)[0]
        zeta = convert_zeta(zeta, self.duct_c, self.duct_b)
        return zeta

    @property
    def zeta_c(self) -> float:
        zeta = self.interp_zeta(self.theta, self.Vb_on_Vc)
        return zeta

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class DivergingJunctionA11A(AbstractFitting):
    description = "Tee or wye, 30° to 90°, round"
    _Ab_on_Ac = np.array([0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1])
    _Vb_on_Vc = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    _vs_on_vc = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0])
    _zeta_b_wye30 = np.array([
        [0.75, 0.55, 0.40, 0.28, 0.21, 0.16, 0.15, 0.16, 0.19],
        [0.72, 0.51, 0.36, 0.25, 0.18, 0.15, 0.16, 0.20, 0.26],
        [0.69, 0.46, 0.31, 0.21, 0.17, 0.16, 0.20, 0.28, 0.39],
        [0.65, 0.41, 0.26, 0.19, 0.18, 0.22, 0.32, 0.47, 0.67],
        [0.59, 0.33, 0.21, 0.20, 0.27, 0.40, 0.62, 0.92, 1.30],
        [0.55, 0.28, 0.24, 0.38, 0.76, 1.30, 2.00, 2.00, 2.00],
        [0.40, 0.26, 0.58, 1.30, 2.50, 2.50, 2.50, 2.50, 2.50],
        [0.28, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50]
    ])
    _zeta_b_wye45 = np.array([
        [0.78, 0.62, 0.49, 0.40, 0.34, 0.31, 0.32, 0.35, 0.40],
        [0.77, 0.59, 0.47, 0.38, 0.34, 0.32, 0.35, 0.41, 0.50],
        [0.74, 0.56, 0.44, 0.37, 0.35, 0.36, 0.43, 0.54, 0.68],
        [0.71, 0.52, 0.41, 0.38, 0.40, 0.45, 0.59, 0.78, 1.00],
        [0.66, 0.47, 0.40, 0.43, 0.54, 0.69, 0.95, 1.30, 1.70],
        [0.66, 0.48, 0.52, 0.73, 1.20, 1.80, 2.70, 2.70, 2.70],
        [0.56, 0.56, 1.00, 1.80, 1.80, 1.80, 1.80, 1.80, 1.80],
        [0.60, 2.10, 2.10, 2.10, 2.10, 2.10, 2.10, 2.10, 2.10]
    ])
    _zeta_b_wye60 = np.array([
        [0.83, 0.71, 0.62, 0.56, 0.52, 0.50, 0.53, 0.60, 0.68],
        [0.82, 0.69, 0.61, 0.56, 0.54, 0.54, 0.60, 0.70, 0.82],
        [0.81, 0.68, 0.60, 0.58, 0.58, 0.61, 0.72, 0.87, 1.10],
        [0.79, 0.66, 0.61, 0.62, 0.68, 0.76, 0.94, 1.20, 1.50],
        [0.76, 0.65, 0.65, 0.74, 0.89, 1.10, 1.40, 1.80, 2.30],
        [0.80, 0.75, 0.89, 1.20, 1.80, 2.60, 3.50, 3.50, 3.50],
        [0.77, 0.96, 1.60, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50],
        [1.00, 2.90, 2.90, 2.90, 2.90, 2.90, 2.90, 2.90, 2.90]
    ])
    _zeta_b_wye90 = np.array([
        [0.95, 0.92, 0.92, 0.93, 0.94, 0.95, 1.10, 1.20, 1.40],
        [0.95, 0.94, 0.95, 0.98, 1.00, 1.10, 1.20, 1.40, 1.60],
        [0.96, 0.97, 1.00, 1.10, 1.10, 1.20, 1.40, 1.70, 2.00],
        [0.97, 1.00, 1.10, 1.20, 1.40, 1.50, 1.80, 2.10, 2.50],
        [0.99, 1.10, 1.30, 1.50, 1.70, 2.00, 2.40, 2.40, 2.40],
        [1.10, 1.40, 1.80, 2.30, 2.30, 2.30, 2.30, 2.30, 2.30],
        [1.30, 1.90, 2.90, 2.90, 2.90, 2.90, 2.90, 2.90, 2.90],
        [2.10, 2.10, 2.10, 2.10, 2.10, 2.10, 2.10, 2.10, 2.10]
    ])
    _zeta_c = np.array([0.35, 0.28, 0.22, 0.17, 0.13, 0.09, 0.06, 0.02, 0.00])

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, theta: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        vs = duct_s.velocity.to('m / s').magnitude
        vc = duct_c.velocity.to('m / s').magnitude
        self.vs_on_vc = vs / vc
        Ab = duct_b.cross_section.area.to('m ** 2').magnitude
        Ac = duct_c.cross_section.area.to('m ** 2').magnitude
        self.Ab_on_Ac = Ab / Ac
        Vb = duct_b.volume_flow_rate.to('m ** 3 / s').magnitude
        Vc = duct_c.volume_flow_rate.to('m ** 3 / s').magnitude
        self.Vb_on_Vc = Vb / Vc
        theta = theta.to('deg')
        self.interp_zeta_c = interpolate.interp1d(
            self._vs_on_vc,
            self._zeta_c,
            bounds_error=False,
            fill_value=(self._zeta_c[0], self._zeta_c[-1])
        )
        if theta == 30 * u.deg:
            self.interp_zeta_b = interpolate.interp2d(self._Vb_on_Vc, self._Ab_on_Ac, self._zeta_b_wye30)
        elif theta == 45 * u.deg:
            self.interp_zeta_b = interpolate.interp2d(self._Vb_on_Vc, self._Ab_on_Ac, self._zeta_b_wye45)
        elif theta == 60 * u.deg:
            self.interp_zeta_b = interpolate.interp2d(self._Vb_on_Vc, self._Ab_on_Ac, self._zeta_b_wye60)
        elif theta == 90 * u.deg:
            self.interp_zeta_b = interpolate.interp2d(self._Vb_on_Vc, self._Ab_on_Ac, self._zeta_b_wye90)
        else:
            self.interp_zeta_b = None

    @property
    def zeta_b(self) -> float:
        zeta = self.interp_zeta_b(self.Vb_on_Vc, self.Ab_on_Ac)[0]
        zeta = convert_zeta(zeta, self.duct_c, self.duct_b)
        return zeta

    @property
    def zeta_c(self) -> float:
        zeta = self.interp_zeta_c(self.vs_on_vc)
        return zeta

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class DivergingJunctionA11B(AbstractFitting):
    description = "90° conical tee, round"
    _vb_on_vc = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    _zeta_b = np.array([1.0, 0.85, 0.74, 0.62, 0.52, 0.42, 0.36, 0.32, 0.32, 0.37, 0.52])
    _vs_on_vc = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0])
    _zeta_c = np.array([0.35, 0.28, 0.22, 0.17, 0.13, 0.09, 0.06, 0.02, 0.00])

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        vb = duct_b.velocity.to('m / s').magnitude
        vc = duct_c.velocity.to('m / s').magnitude
        vs = duct_s.velocity.to('m / s').magnitude
        self.vb_on_vc = vb / vc
        self.vs_on_vc = vs / vc
        self.interp_zeta_b = interpolate.interp1d(self._vb_on_vc, self._zeta_b)
        self.interp_zeta_c = interpolate.interp1d(self._vs_on_vc, self._zeta_c)

    @property
    def zeta_b(self) -> float:
        zeta = self.interp_zeta_b(self.vb_on_vc)
        zeta = convert_zeta(zeta, self.duct_c, self.duct_b)
        return zeta

    @property
    def zeta_c(self) -> float:
        zeta = self.interp_zeta_c(self.vs_on_vc)
        return zeta

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class DivergingJunctionA11C(DivergingJunctionA11B):
    description = "45° conical wye, round"
    _zeta_b = np.array([1.0, 0.84, 0.61, 0.41, 0.27, 0.17, 0.12, 0.12, 0.14, 0.18, 0.27])


class DivergingJunctionA11D(DivergingJunctionA11B):
    description = "90° tee, round, rolled 45° with 45° elbow, branch 90° to main"
    _zeta_b = np.array([1.0, 1.32, 1.51, 1.60, 1.65, 1.74, 1.87, 2.0, 2.2, 2.5, 2.7])


class DivergingJunctionA11E(DivergingJunctionA11B):
    description = "90° tee, round, with 90° elbow, branch 90° to main"
    _zeta_b = np.array([1.0, 1.03, 1.08, 1.18, 1.33, 1.56, 1.86, 2.20, 2.60, 3.00, 3.40])


class DivergingJunctionA11F(DivergingJunctionA11B):
    description = "90° tee, round, rolled 45° with 60° elbow, branch 45° to main"
    _zeta_b = np.array([1.0, 1.06, 1.15, 1.29, 1.45, 1.65, 1.89, 2.20, 2.50, 2.90, 3.30])


class DivergingJunctionA11G(DivergingJunctionA11B):
    description = "90° conical tee, round, rolled 45° with 45° elbow, branch 90° to main"
    _zeta_b = np.array([1.0, 0.94, 0.88, 0.84, 0.80, 0.82, 0.84, 0.87, 0.90, 0.95, 1.02])


class DivergingJunctionA11H(DivergingJunctionA11B):
    description = "90° conical tee, round, rolled 45° with 60° elbow, branch 45° to main"
    _zeta_b = np.array([1.0, 0.95, 0.90, 0.86, 0.81, 0.79, 0.79, 0.81, 0.86, 0.96, 1.10])


class DivergingJunctionA11I(DivergingJunctionA11B):
    description = "45° wye, round, rolled 45° with 60° elbow, branch 90° to main"
    _zeta_b = np.array([1.0, 0.88, 0.77, 0.68, 0.65, 0.69, 0.73, 0.88, 1.14, 1.54, 2.20])


class DivergingJunctionA11J(DivergingJunctionA11B):
    description = "45° conical wye, round, rolled 45° with 60° elbow, branch 90° to main"
    _zeta_b = np.array([1.0, 0.82, 0.63, 0.52, 0.45, 0.42, 0.41, 0.40, 0.41, 0.45, 0.56])


class DivergingJunctionA11K(DivergingJunctionA11B):
    description = "45° wye, round, rolled 45° with 30° elbow, branch 45° to main"
    _zeta_b = np.array([1.0, 0.84, 0.72, 0.62, 0.54, 0.50, 0.56, 0.71, 0.92, 1.22, 1.66])


class DivergingJunctionA11L(DivergingJunctionA11B):
    description = "45° conical wye, round, rolled 45° with 30° elbow, branch 45° to main"
    _zeta_b = np.array([1.0, 0.93, 0.71, 0.55, 0.44, 0.42, 0.42, 0.44, 0.47, 0.54, 0.62])


class DivergingJunctionA11M(AbstractFitting):
    description = "45° conical main and branch with 45° elbow, branch 90° to main"
    _vb_on_vc = np.array([0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
                          1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0])
    _zeta_b = np.array([0.76, 0.60, 0.52, 0.50, 0.51, 0.52, 0.56, 0.61, 0.68,
                        0.86, 1.1, 1.4, 1.8, 2.2, 2.6, 3.1, 3.7, 4.2])
    _vs_on_vc = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    _zeta_c = np.array([0.14, 0.06, 0.05, 0.09, 0.18, 0.30, 0.46, 0.64, 0.84, 1.0])

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        vb = duct_b.velocity.to('m / s').magnitude
        vc = duct_c.velocity.to('m / s').magnitude
        vs = duct_s.velocity.to('m / s').magnitude
        self.vb_on_vc = vb / vc
        self.vs_on_vc = vs / vc
        self.interp_zeta_b = interpolate.interp1d(self._vb_on_vc, self._zeta_b)
        self.interp_zeta_c = interpolate.interp1d(self._vs_on_vc, self._zeta_c)

    @property
    def zeta_b(self) -> float:
        zeta = self.interp_zeta_b(self.vb_on_vc)
        zeta = convert_zeta(zeta, self.duct_c, self.duct_b)
        return zeta

    @property
    def zeta_c(self) -> float:
        zeta = self.interp_zeta_c(self.vs_on_vc)
        return zeta

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class DivergingJunctionA11N(AbstractFitting):
    description = "tee, 45° rectangular main and branch"
    _vb_on_vc = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8])
    _Vb_on_Vc = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    _zeta_b = np.array([
        [0.91, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.81, 0.79, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.77, 0.72, 0.70, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.78, 0.73, 0.69, 0.66, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.78, 0.98, 0.85, 0.79, 0.74, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.90, 1.11, 1.16, 1.23, 1.03, 0.86, 9999.9, 9999.9, 9999.9],
        [1.19, 1.22, 1.26, 1.29, 1.54, 1.25, 0.92, 9999.9, 9999.9],
        [1.35, 1.42, 1.55, 1.59, 1.63, 1.50, 1.31, 1.09, 9999.9],
        [1.44, 1.50, 1.75, 1.74, 1.72, 2.24, 1.63, 1.40, 1.17]
    ])
    _vs_on_vc = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0])
    _zeta_c = np.array([0.35, 0.28, 0.22, 0.17, 0.13, 0.09, 0.06, 0.02, 0.00])

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        vb = duct_b.velocity.to('m / s').magnitude
        vc = duct_c.velocity.to('m / s').magnitude
        vs = duct_s.velocity.to('m / s').magnitude
        Vb = duct_b.volume_flow_rate.to('m ** 3 / s').magnitude
        Vc = duct_c.volume_flow_rate.to('m ** 3 / s').magnitude
        self.vb_on_vc = vb / vc
        self.Vb_on_Vc = Vb / Vc
        self.interp_zeta_b = interpolate.interp2d(self._Vb_on_Vc, self._vb_on_vc, self._zeta_b)
        self.vs_on_vc = vs / vc
        self.interp_zeta_c = interpolate.interp1d(self._vs_on_vc, self._zeta_c)

    @property
    def zeta_b(self) -> float:
        zeta = self.interp_zeta_b(self.Vb_on_Vc, self.vb_on_vc)[0]
        if zeta < 100.0:
            zeta = convert_zeta(zeta, self.duct_c, self.duct_b)
            return zeta
        else:
            return np.nan

    @property
    def zeta_c(self) -> float:
        zeta = self.interp_zeta_c(self.vs_on_vc)
        return zeta

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class DivergingJunctionA11R(AbstractFitting):
    description = "tee, rectangular main and branch with extractor"
    _vb_on_vc = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8])
    _Vb_on_Vc = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    _zeta_b = np.array([
        [0.60, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.62, 0.69, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.74, 0.80, 0.82, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.99, 1.10, 0.95, 0.90, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [1.48, 1.12, 1.41, 1.24, 1.21, 9999.9, 9999.9, 9999.9, 9999.9],
        [1.91, 1.33, 1.43, 1.52, 1.55, 1.64, 9999.9, 9999.9, 9999.9],
        [2.47, 1.67, 1.70, 2.04, 1.86, 1.98, 2.47, 9999.9, 9999.9],
        [3.17, 2.40, 2.33, 2.53, 2.31, 2.51, 3.13, 3.25, 9999.9],
        [3.85, 3.37, 2.89, 3.23, 3.09, 3.03, 3.30, 3.74, 4.11]
    ])
    _zeta_c = np.array([0.03, 0.04, 0.07, 0.12, 0.13, 0.14, 0.27, 0.30, 0.25])

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        vb = duct_b.velocity.to('m / s').magnitude
        vc = duct_c.velocity.to('m / s').magnitude
        Vb = duct_b.volume_flow_rate.to('m ** 3 / s').magnitude
        Vc = duct_c.volume_flow_rate.to('m ** 3 / s').magnitude
        self.vb_on_vc = vb / vc
        self.Vb_on_Vc = Vb / Vc
        self.interp_zeta_b = interpolate.interp2d(self._Vb_on_Vc, self._vb_on_vc, self._zeta_b)
        self.interp_zeta_c = interpolate.interp1d(self._vb_on_vc, self._zeta_c)

    @property
    def zeta_b(self) -> float:
        zeta = self.interp_zeta_b(self.Vb_on_Vc, self.vb_on_vc)[0]
        if zeta < 100:
            zeta = convert_zeta(zeta, self.duct_c, self.duct_b)
            return zeta
        else:
            return np.nan

    @property
    def zeta_c(self) -> float:
        zeta = self.interp_zeta_c(self.vb_on_vc)
        return zeta

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class DivergingJunctionA11O(DivergingJunctionA11R):
    description = "tee, 45° entry, rectangular main and branch with damper"
    _zeta_b = np.array([
        [0.61, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.46, 0.61, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.43, 0.50, 0.54, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.39, 0.43, 0.62, 0.53, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.34, 0.57, 0.77, 0.73, 0.68, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.37, 0.64, 0.85, 0.98, 1.07, 0.83, 9999.9, 9999.9, 9999.9],
        [0.57, 0.71, 1.04, 1.16, 1.54, 1.36, 1.18, 9999.9, 9999.9],
        [0.89, 1.08, 1.28, 1.30, 1.69, 2.09, 1.81, 1.47, 9999.9],
        [1.33, 1.34, 2.04, 1.78, 1.90, 2.40, 2.77, 2.23, 1.92]
    ])


class DivergingJunctionA11P(DivergingJunctionA11N):
    description = "tee, rectangular main and branch"
    _zeta_b = np.array([
        [1.03, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [1.04, 1.01, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [1.11, 1.03, 1.05, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [1.16, 1.21, 1.17, 1.12, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [1.38, 1.40, 1.30, 1.36, 1.27, 9999.9, 9999.9, 9999.9, 9999.9],
        [1.52, 1.61, 1.68, 1.91, 1.47, 1.66, 9999.9, 9999.9, 9999.9],
        [1.79, 2.01, 1.90, 2.31, 2.28, 2.20, 1.95, 9999.9, 9999.9],
        [2.07, 2.28, 2.13, 2.71, 2.99, 2.81, 2.09, 2.20, 9999.9],
        [2.32, 2.54, 2.64, 3.09, 3.72, 3.48, 2.21, 2.29, 2.57]
    ])


class DivergingJunctionA11Q(DivergingJunctionA11R):
    description = "tee, rectangular main and branch with damper"
    _zeta_b = np.array([
        [0.58, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.67, 0.64, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.78, 0.76, 0.75, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [0.88, 0.98, 0.81, 1.01, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [1.12, 1.05, 1.08, 1.18, 1.29, 9999.9, 9999.9, 9999.9, 9999.9],
        [1.49, 1.48, 1.40, 1.51, 1.70, 1.91, 9999.9, 9999.9, 9999.9],
        [2.10, 2.21, 2.25, 2.29, 2.32, 2.48, 2.53, 9999.9, 9999.9],
        [2.72, 3.30, 2.84, 3.09, 3.30, 3.19, 3.29, 3.16, 9999.9],
        [3.42, 4.58, 3.65, 3.92, 4.20, 4.15, 4.14, 4.10, 4.05]
    ])


class DivergingJunctionA11S(DivergingJunctionA11N):
    description = "tee, rectangular main to round branch"
    _zeta_b = np.array([
        [1.00, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [1.01, 1.07, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [1.14, 1.10, 1.08, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [1.18, 1.31, 1.12, 1.13, 9999.9, 9999.9, 9999.9, 9999.9, 9999.9],
        [1.30, 1.38, 1.20, 1.23, 1.26, 9999.9, 9999.9, 9999.9, 9999.9],
        [1.46, 1.58, 1.45, 1.31, 1.39, 1.48, 9999.9, 9999.9, 9999.9],
        [1.70, 1.82, 1.65, 1.51, 1.56, 1.64, 1.71, 9999.9, 9999.9],
        [1.93, 2.06, 2.00, 1.85, 1.70, 1.76, 1.80, 1.88, 9999.9],
        [2.06, 2.17, 2.20, 2.13, 2.06, 1.98, 1.99, 2.00, 2.07]
    ])


class DivergingJunctionA11T(AbstractFitting):
    description = "wye, rectangular"
    _vb_on_vc = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    _vs_on_vc = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    _theta = np.array([15, 30, 45, 60, 90])
    _As_on_Ac = np.array([0.4, 0.5, 0.6, 0.7, 0.8])
    _zeta_b = np.array([
        [0.81, 0.65, 0.51, 0.38, 0.28, 0.20, 0.11, 0.06, 0.14, 0.30, 0.51, 0.76, 1.00],
        [0.84, 0.69, 0.56, 0.44, 0.34, 0.26, 0.19, 0.15, 0.15, 0.30, 0.51, 0.76, 1.00],
        [0.87, 0.74, 0.63, 0.54, 0.45, 0.38, 0.29, 0.24, 0.23, 0.30, 0.51, 0.76, 1.00],
        [0.90, 0.82, 0.79, 0.66, 0.59, 0.53, 0.43, 0.36, 0.33, 0.39, 0.51, 0.76, 1.00],
        [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]
    ])
    _zeta_c_th15_60 = np.array([1.0, 0.81, 0.64, 0.50, 0.36, 0.25, 0.16, 0.04, 0.00, 0.07, 0.39, 0.90, 1.80, 3.20])
    _zeta_c_th90 = np.array([
        [1.00, 1.00, 1.00, 1.00, 1.00],
        [0.81, 0.81, 0.81, 0.81, 0.81],
        [0.64, 0.64, 0.64, 0.64, 0.64],
        [0.50, 0.52, 0.52, 0.50, 0.50],
        [0.36, 0.40, 0.38, 0.37, 0.36],
        [0.25, 0.30, 0.28, 0.27, 0.25],
        [0.16, 0.23, 0.20, 0.18, 0.16],
        [0.04, 0.17, 0.10, 0.07, 0.04],
        [0.00, 0.20, 0.10, 0.05, 0.00],
        [0.07, 0.36, 0.21, 0.14, 0.07],
        [0.39, 0.79, 0.59, 0.39, 9999.9],
        [0.90, 1.40, 1.20, 9999.9, 9999.9],
        [1.80, 2.40, 9999.9, 9999.9, 9999.9],
        [3.20, 4.00, 9999.9, 9999.9, 9999.9]
    ])

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, theta: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        vb = duct_b.velocity.to('m / s').magnitude
        vc = duct_c.velocity.to('m / s').magnitude
        vs = duct_s.velocity.to('m / s').magnitude
        As = duct_s.cross_section.area.to('m ** 2').magnitude
        Ac = duct_c.cross_section.area.to('m ** 2').magnitude
        self.As_on_Ac = As / Ac
        self.vb_on_vc = vb / vc
        self.vs_on_vc = vs / vc
        self.theta = theta.to('deg').magnitude
        self.interp_zeta_b = interpolate.interp2d(self._vb_on_vc, self._theta, self._zeta_b)
        if 15 <= self.theta <= 60:
            self.interp_zeta_c_1d = interpolate.interp1d(self._vs_on_vc, self._zeta_c_th15_60)
        if self.theta == 90:
            self.interp_zeta_c_2d = interpolate.interp2d(self._As_on_Ac, self._vs_on_vc, self._zeta_c_th90)

    @property
    def zeta_b(self) -> float:
        zeta = self.interp_zeta_b(self.vb_on_vc, self.theta)[0]
        zeta = convert_zeta(zeta, self.duct_c, self.duct_b)
        return zeta

    @property
    def zeta_c(self) -> float:
        if 15 <= self.theta <= 60:
            zeta = self.interp_zeta_c_1d(self.vs_on_vc)
        else:
            zeta = self.interp_zeta_c_2d(self.As_on_Ac, self.vs_on_vc)[0]
        return zeta

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class DivergingJunctionA11U(AbstractFitting):
    description = "tee, rectangular main to conical branch"
    _vb_on_vc = np.array([0.40, 0.50, 0.75, 1.00, 1.30, 1.50])
    _zeta_b = np.array([0.80, 0.83, 0.90, 1.00, 1.10, 1.40])
    _vs_on_vc = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0])
    _zeta_c = np.array([0.35, 0.28, 0.22, 0.17, 0.13, 0.09, 0.06, 0.02, 0.00])

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        vb = duct_b.velocity.to('m / s').magnitude
        vc = duct_c.velocity.to('m / s').magnitude
        vs = duct_s.velocity.to('m / s').magnitude
        self.vb_on_vc = vb / vc
        self.interp_zeta_b = interpolate.interp1d(self._vb_on_vc, self._zeta_b)
        self.vs_on_vc = vs / vc
        self.interp_zeta_c = interpolate.interp1d(self._vs_on_vc, self._zeta_c)

    @property
    def zeta_b(self) -> float:
        zeta = self.interp_zeta_b(self.vb_on_vc)
        zeta = convert_zeta(zeta, self.duct_c, self.duct_b)
        return zeta

    @property
    def zeta_c(self) -> float:
        zeta = self.interp_zeta_c(self.vs_on_vc)
        return zeta

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class DivergingJunctionA11V(AbstractFitting):
    description = "wye, rectangular"
    _Vb_on_Vc = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    _zeta_b = np.array([
        [0.55, 0.50, 0.60, 0.85, 1.20, 1.80, 3.10, 4.40, 6.00],
        [0.35, 0.35, 0.50, 0.80, 1.30, 2.00, 2.80, 3.80, 5.00],
        [0.62, 0.48, 0.40, 0.40, 0.48, 0.60, 0.78, 1.10, 1.50],
        [0.52, 0.40, 0.32, 0.30, 0.34, 0.44, 0.62, 0.92, 1.40],
        [0.44, 0.38, 0.38, 0.41, 0.52, 0.68, 0.92, 1.20, 1.60],
        [0.67, 0.55, 0.46, 0.37, 0.32, 0.29, 0.29, 0.30, 0.37],
        [0.70, 0.60, 0.51, 0.42, 0.34, 0.28, 0.26, 0.26, 0.29],
        [0.60, 0.52, 0.43, 0.33, 0.24, 0.17, 0.15, 0.17, 0.21]
    ])
    _zeta_c = np.array([
        [-0.01, -0.03, -0.01, 0.05, 0.13, 0.21, 0.29, 0.38, 0.46],
        [0.08, 0.00, -0.02, -0.01, 0.02, 0.08, 0.16, 0.24, 0.34],
        [-0.03, -0.06, -0.05, 0.00, 0.06, 0.12, 0.19, 0.27, 0.35],
        [0.04, -0.02, -0.04, -0.03, -0.01, 0.04, 0.12, 0.23, 0.37],
        [0.72, 0.48, 0.28, 0.13, 0.05, 0.04, 0.09, 0.18, 0.30],
        [-0.02, -0.04, -0.04, -0.01, 0.06, 0.13, 0.22, 0.30, 0.38],
        [0.10, 0.00, 0.01, -0.03, -0.01, 0.03, 0.10, 0.20, 0.30],
        [0.62, 0.38, 0.23, 0.13, 0.08, 0.05, 0.06, 0.10, 0.20]
    ])

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        Ab = duct_b.cross_section.area.to('m ** 2').magnitude
        Ac = duct_c.cross_section.area.to('m ** 2').magnitude
        As = duct_s.cross_section.area.to('m ** 2').magnitude
        Vb = duct_b.volume_flow_rate.to('m ** 3 / s').magnitude
        Vc = duct_c.volume_flow_rate.to('m ** 3 / s').magnitude
        self.Vb_on_Vc = Vb / Vc
        Ab_on_As = round(Ab / As, 2)
        Ab_on_Ac = round(Ab / Ac, 2)
        if (Ab_on_Ac, Ab_on_As) == (0.25, 0.25):
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[0])
            self.interp_zeta_c = interpolate.interp1d(self._Vb_on_Vc, self._zeta_c[0])
        elif (Ab_on_Ac, Ab_on_As) == (0.25, 0.33):
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[1])
            self.interp_zeta_c = interpolate.interp1d(self._Vb_on_Vc, self._zeta_c[1])
        elif (Ab_on_Ac, Ab_on_As) == (0.50, 0.50):
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[2])
            self.interp_zeta_c = interpolate.interp1d(self._Vb_on_Vc, self._zeta_c[2])
        elif (Ab_on_Ac, Ab_on_As) == (0.50, 0.67):
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[3])
            self.interp_zeta_c = interpolate.interp1d(self._Vb_on_Vc, self._zeta_c[3])
        elif (Ab_on_Ac, Ab_on_As) == (0.50, 1.00):
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[4])
            self.interp_zeta_c = interpolate.interp1d(self._Vb_on_Vc, self._zeta_c[4])
        elif (Ab_on_Ac, Ab_on_As) == (1.00, 1.00):
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[5])
            self.interp_zeta_c = interpolate.interp1d(self._Vb_on_Vc, self._zeta_c[5])
        elif (Ab_on_Ac, Ab_on_As) == (1.00, 1.33):
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[6])
            self.interp_zeta_c = interpolate.interp1d(self._Vb_on_Vc, self._zeta_c[6])
        elif (Ab_on_Ac, Ab_on_As) == (1.00, 2.00):
            self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b[7])
            self.interp_zeta_c = interpolate.interp1d(self._Vb_on_Vc, self._zeta_c[7])
        else:
            self.interp_zeta_b = None
            self.interp_zeta_c = None

    @property
    def zeta_b(self) -> float:
        if self.interp_zeta_b:
            zeta = self.interp_zeta_b(self.Vb_on_Vc)
            zeta = convert_zeta(zeta, self.duct_c, self.duct_b)
            return zeta
        else:
            return float('nan')

    @property
    def zeta_c(self) -> float:
        if self.interp_zeta_c:
            zeta = self.interp_zeta_c(self.Vb_on_Vc)
            return zeta
        else:
            return float('nan')

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class DivergingJunctionA11W(AbstractFitting):
    description = "symmetrical wye, dovetail, rectangular"

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        Ab = duct_b.cross_section.area.to('m ** 2').magnitude
        Ac = duct_c.cross_section.area.to('m ** 2').magnitude
        Ab_on_Ac = round(Ab / Ac, 1)
        if Ab_on_Ac == 0.5:
            self._zeta = 0.30
        elif Ab_on_Ac == 1.0:
            self._zeta = 0.25
        else:
            self._zeta = None

    @property
    def zeta_b(self) -> float:
        if self._zeta:
            zeta = convert_zeta(self._zeta, self.duct_c, self.duct_b)
            return zeta
        else:
            return float('nan')

    @property
    def zeta_c(self) -> float:
        if self._zeta:
            return self._zeta
        else:
            return float('nan')

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class DivergingJunctionA11X(AbstractFitting):
    description = "wye, rectangular and round"
    _vb_on_vc = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    _theta = np.array([15, 30, 45, 60, 90])
    _zeta = np.array([
        [0.81, 0.65, 0.51, 0.38, 0.28, 0.20, 0.11, 0.06, 0.14, 0.30, 0.51, 0.76, 1.00],
        [0.84, 0.69, 0.56, 0.44, 0.34, 0.26, 0.19, 0.15, 0.15, 0.30, 0.51, 0.76, 1.00],
        [0.87, 0.74, 0.63, 0.54, 0.45, 0.38, 0.29, 0.24, 0.23, 0.30, 0.51, 0.76, 1.00],
        [0.90, 0.82, 0.79, 0.66, 0.59, 0.53, 0.43, 0.36, 0.33, 0.39, 0.51, 0.76, 1.00],
        [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]
    ])

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, theta: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        self.theta = theta.to('deg').magnitude
        vb = duct_b.velocity.to('m / s').magnitude
        vc = duct_c.velocity.to('m / s').magnitude
        self.vb_on_vc = vb / vc
        self.interp_zeta = interpolate.interp2d(self._vb_on_vc, self._theta, self._zeta)

    @property
    def zeta_b(self) -> float:
        zeta = self.interp_zeta(self.vb_on_vc, self.theta)[0]
        zeta = convert_zeta(zeta, self.duct_c, self.duct_b)
        return zeta

    @property
    def zeta_c(self) -> float:
        zeta = self.interp_zeta(self.vb_on_vc, self.theta)[0]
        return zeta

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class DivergingJunctionA11Y(AbstractFitting):
    description = "tee, rectangular reducing, 45° branch"
    _Vb_on_Vc = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    _Vs_on_Vc = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    _zeta_b = np.array([1.16, 0.96, 0.82, 0.68, 0.56, 0.49, 0.47, 0.48, 0.50, 0.54])
    _zeta_m = np.array([0.21, 0.20, 0.20, 0.20, 0.20, 0.20, 0.22, 0.25, 0.35, 0.53])

    def __init__(self, duct_b: Duct, duct_c: Duct, duct_s: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct_b = duct_b
        self.duct_c = duct_c
        self.duct_s = duct_s
        Vb = duct_b.volume_flow_rate.to('m ** 3 / s').magnitude
        Vc = duct_c.volume_flow_rate.to('m ** 3 / s').magnitude
        Vs = duct_s.volume_flow_rate.to('m ** 3 / s').magnitude
        self.Vb_on_Vc = Vb / Vc
        self.Vs_on_Vc = Vs / Vc
        self.interp_zeta_b = interpolate.interp1d(self._Vb_on_Vc, self._zeta_b)
        self.interp_zeta_m = interpolate.interp1d(self._Vs_on_Vc, self._zeta_m)

    @property
    def zeta_b(self) -> float:
        zeta = self.interp_zeta_b(self.Vb_on_Vc)
        zeta = convert_zeta(zeta, self.duct_c, self.duct_b)
        return zeta

    @property
    def zeta_c(self) -> float:
        zeta = self.interp_zeta_m(self.Vs_on_Vc)
        return zeta

    @property
    def zeta_s(self) -> float:
        # refer `zeta_s` to velocity in straight duct
        zeta = convert_zeta(self.zeta_c, self.duct_c, self.duct_s)
        return zeta

    @property
    def pressure_drop_branch(self) -> Quantity:
        Pv = self.duct_b.velocity_pressure.to('Pa')
        return self.zeta_b * Pv

    @property
    def pressure_drop_main(self) -> Quantity:
        Pv = self.duct_c.velocity_pressure.to('Pa')
        return self.zeta_c * Pv


class EntryA12A(AbstractFitting):
    description = "duct mounted in wall, round and rectangular"
    _L_on_D = np.array([0.000, 0.002, 0.010, 0.050, 0.200, 0.500, 1.000])
    _zeta = np.array([
        [0.50, 0.57, 0.68, 0.80, 0.92, 1.00, 1.00],
        [0.50, 0.51, 0.52, 0.55, 0.66, 0.72, 0.72],
        [0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50]
    ])

    def __init__(self, duct: Duct, thickness: Quantity = Q_(0.0, 'mm'), ID: str = ''):
        super().__init__(ID)
        self.duct = duct
        D = duct.cross_section.hydraulic_diameter.to('mm').magnitude
        t = thickness.to('mm').magnitude
        L = duct.length.to('mm').magnitude
        t_on_D = round(t / D, 2)
        self.L_on_D = L / D
        if t_on_D == 0.00:
            self.interp_zeta = interpolate.interp1d(self._L_on_D, self._zeta[0])
        elif 0.00 < t_on_D < 0.05:
            self.interp_zeta = interpolate.interp1d(self._L_on_D, self._zeta[1])
        elif t_on_D >= 0.05:
            self.interp_zeta = interpolate.interp1d(self._L_on_D, self._zeta[2])
        else:
            self.interp_zeta = None

    @property
    def zeta(self) -> float:
        if self.interp_zeta:
            zeta = self.interp_zeta(self.L_on_D)
            return zeta
        else:
            return float('nan')

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct.velocity_pressure.to('Pa')
        return self.zeta * Pv


class EntryA12B(AbstractFitting):
    description = "smooth covering bell mouth, round, without end wall"
    _R_on_D = np.array([0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20])
    _zeta = np.array([1.0, 0.87, 0.74, 0.61, 0.51, 0.40, 0.32, 0.20, 0.15, 0.10, 0.06, 0.03])

    def __init__(self, duct: Duct, radius: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct = duct
        R = radius.to('mm').magnitude
        D = duct.cross_section.hydraulic_diameter.to('mm').magnitude
        self.R_on_D = R / D
        self.interp_zeta = interpolate.interp1d(self._R_on_D, self._zeta)

    @property
    def zeta(self) -> float:
        zeta = self.interp_zeta(self.R_on_D)
        return zeta

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct.velocity_pressure.to('Pa')
        return self.zeta * Pv


class EntryA12C(EntryA12B):
    description = "smooth covering bellmouth, round, with end wall"
    _zeta = np.array([0.50, 0.43, 0.36, 0.31, 0.26, 0.22, 0.20, 0.15, 0.12, 0.09, 0.06, 0.03])


class EntryA12D(AbstractFitting):
    description = "conical, converging bellmouth, round and rectangular without end wall"
    _L_on_D = np.array([0.025, 0.05, 0.10, 0.25, 0.60, 1.0])
    _theta = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 60.0, 100.0, 140.0, 180.0])
    _zeta = np.array([
        [1.0, 0.96, 0.93, 0.90, 0.86, 0.80, 0.69, 0.59, 0.50],
        [1.0, 0.93, 0.86, 0.80, 0.75, 0.67, 0.58, 0.53, 0.50],
        [1.0, 0.80, 0.67, 0.55, 0.48, 0.41, 0.41, 0.44, 0.50],
        [1.0, 0.68, 0.45, 0.30, 0.22, 0.17, 0.22, 0.34, 0.50],
        [1.0, 0.46, 0.27, 0.18, 0.14, 0.13, 0.21, 0.33, 0.50],
        [1.0, 0.32, 0.20, 0.14, 0.11, 0.10, 0.18, 0.30, 0.50]
    ])

    def __init__(self, duct: Duct, theta: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct = duct
        self.theta = theta.to('deg').magnitude
        L = duct.length.to('mm').magnitude
        D = duct.cross_section.hydraulic_diameter.to('mm').magnitude
        self.L_on_D = L / D
        self.interp_zeta = interpolate.interp2d(self._theta, self._L_on_D, self._zeta)

    @property
    def zeta(self) -> float:
        zeta = self.interp_zeta(self.theta, self.L_on_D)[0]
        return zeta

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct.velocity_pressure.to('Pa')
        return self.zeta * Pv


class EntryA12E(EntryA12D):
    description = "conical, converging bellmouth, round and rectangular with end wall"
    _L_on_D = np.array([0.025, 0.05, 0.075, 0.10, 0.15, 0.60])
    _zeta = np.array([
        [0.50, 0.47, 0.45, 0.43, 0.41, 0.40, 0.42, 0.45, 0.50],
        [0.50, 0.45, 0.41, 0.36, 0.33, 0.30, 0.35, 0.42, 0.50],
        [0.50, 0.42, 0.35, 0.30, 0.26, 0.23, 0.30, 0.40, 0.50],
        [0.50, 0.39, 0.32, 0.25, 0.22, 0.18, 0.27, 0.38, 0.50],
        [0.50, 0.37, 0.27, 0.20, 0.16, 0.15, 0.25, 0.37, 0.50],
        [0.50, 0.27, 0.18, 0.13, 0.11, 0.12, 0.23, 0.36, 0.50]
    ])


class EntryA12F(AbstractFitting):
    description = "intake hood"
    _L_on_D = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    _zeta = np.array([
        [2.6, 1.8, 1.5, 1.4, 1.3, 1.2, 1.2, 1.1, 1.1],
        [1.3, 0.77, 0.60, 0.48, 0.41, 0.30, 0.29, 0.28, 0.25]
    ])

    def __init__(self, duct: Duct, theta: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct = duct
        L = duct.length.to('mm').magnitude
        D = duct.cross_section.hydraulic_diameter.to('mm').magnitude
        self.L_on_D = L / D
        if theta == 0.0 * u.deg:
            self.interp_zeta = interpolate.interp1d(self._L_on_D, self._zeta[0])
        elif theta == 15.0 * u.deg:
            self.interp_zeta = interpolate.interp1d(self._L_on_D, self._zeta[1])
        else:
            self.interp_zeta = None

    @property
    def zeta(self) -> float:
        if self.interp_zeta:
            zeta = self.interp_zeta(self.L_on_D)
            return zeta
        else:
            return float('nan')

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct.velocity_pressure.to('Pa')
        return self.zeta * Pv


class EntryA12G(AbstractFitting):
    description = "hood, tapered, flanged or unflanged"
    _theta = np.array([0, 20, 40, 60, 80, 100, 120, 140, 160, 180])
    _zeta_round_hood = np.array([1.0, 0.11, 0.06, 0.09, 0.14, 0.18, 0.27, 0.32, 0.43, 0.50])
    _zeta_rect_hood = np.array([1.0, 0.19, 0.13, 0.16, 0.21, 0.27, 0.33, 0.43, 0.53, 0.62])

    def __init__(self, duct: Duct, theta: Quantity, shape: str, ID: str = ''):
        super().__init__(ID)
        self.duct = duct
        self.theta = theta.to('deg').magnitude
        if shape == 'round':
            self.interp_zeta = interpolate.interp1d(self._theta, self._zeta_round_hood)
        elif shape == 'rectangular':
            self.interp_zeta = interpolate.interp1d(self._theta, self._zeta_rect_hood)
        else:
            self.interp_zeta = None

    @property
    def zeta(self) -> float:
        if self.interp_zeta:
            zeta = self.interp_zeta(self.theta)
            return zeta
        else:
            return float('nan')

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct.velocity_pressure.to('Pa')
        return self.zeta * Pv


class ExitA13A(AbstractFitting):
    description = "exhaust hood"
    _L_on_D = np.array([0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.8, 1.0])
    _zeta = np.array([
        [4.0, 2.3, 1.9, 1.6, 1.4, 1.3, 1.2, 1.1, 1.0, 1.0],
        [2.6, 1.2, 1.0, 0.80, 0.70, 0.65, 0.60, 0.60, 0.60, 0.60]
    ])

    def __init__(self, duct: Duct, theta: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct = duct
        L = duct.length.to('mm').magnitude
        D = duct.cross_section.hydraulic_diameter.to('mm').magnitude
        self.L_on_D = L / D
        if theta == 0.0 * u.deg:
            self.interp_zeta = interpolate.interp1d(self._L_on_D, self._zeta[0])
        elif theta == 15.0 * u.deg:
            self.interp_zeta = interpolate.interp1d(self._L_on_D, self._zeta[1])
        else:
            self.interp_zeta = None

    @property
    def zeta(self) -> float:
        if self.interp_zeta:
            zeta = self.interp_zeta(self.L_on_D)
            return zeta
        else:
            return float('nan')

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct.velocity_pressure.to('Pa')
        return self.zeta * Pv


class ExitA13B(AbstractFitting):
    description = "exit, conical, round, with or without a wall"
    _theta = np.array([14, 16, 20, 30, 45, 60, 90])
    _As_on_A = np.array([2, 4, 6, 10, 16])
    _zeta = np.array([
        [0.33, 0.36, 0.44, 0.74, 0.97, 0.99, 1.0],
        [0.24, 0.28, 0.36, 0.54, 0.94, 1.0, 1.0],
        [0.22, 0.25, 0.32, 0.49, 0.94, 0.98, 1.0],
        [0.19, 0.23, 0.30, 0.50, 0.94, 0.72, 1.0],
        [0.17, 0.20, 0.27, 0.49, 0.94, 1.0, 1.0]
    ])

    def __init__(self, duct: Duct, theta: Quantity, D_exit: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct = duct
        self.theta = theta.to('deg').magnitude
        D_exit = D_exit.to('m').magnitude
        As = np.pi * D_exit ** 2 / 4
        A = duct.cross_section.area.to('m ** 2').magnitude
        self.As_on_A = As / A
        self.interp_zeta = interpolate.interp2d(self._theta, self._As_on_A, self._zeta)

    @property
    def zeta(self) -> float:
        zeta = self.interp_zeta(self.theta, self.As_on_A)[0]
        return zeta

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct.velocity_pressure.to('Pa')
        return self.zeta * Pv


class ExitA13C(AbstractFitting):
    description = "exit, plane diffuser, rectangular with or without wall"
    _theta = np.array([14, 20, 30, 45, 60, 90])
    _As_on_A = np.array([2, 4, 6])
    _zeta = np.array([
        [0.37, 0.38, 0.50, 0.75, 0.90, 1.1],
        [0.25, 0.37, 0.57, 0.82, 1.0, 1.1],
        [0.28, 0.47, 0.64, 0.87, 1.0, 1.1]
    ])

    def __init__(self, duct: Duct, theta: Quantity, width_exit: Quantity, height_exit: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct = duct
        self.theta = theta.to('deg').magnitude
        W = width_exit.to('m').magnitude
        H = height_exit.to('m').magnitude
        As = W * H
        A = duct.cross_section.area.to('m ** 2').magnitude
        self.As_on_A = As / A
        self.interp_zeta = interpolate.interp2d(self._theta, self._As_on_A, self._zeta)

    @property
    def zeta(self) -> float:
        zeta = self.interp_zeta(self.theta, self.As_on_A)[0]
        return zeta

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct.velocity_pressure.to('Pa')
        return self.zeta * Pv


class ExitA13D(ExitA13C):
    description = "exit, pyramidal diffuser, rectangular with or without wall"
    _theta = np.array([10, 14, 20, 30, 45, 60])
    _As_on_A = np.array([2, 4, 6, 10])
    _zeta = np.array([
        [0.44, 0.58, 0.70, 0.86, 1.0, 1.1],
        [0.31, 0.48, 0.61, 0.76, 0.94, 1.1],
        [0.29, 0.47, 0.62, 0.74, 0.94, 1.1],
        [0.26, 0.45, 0.60, 0.73, 0.89, 1.0]
    ])


class ExitA13J(AbstractFitting):
    description = "exit, abrupt, round and rectangular, with or without wall"

    def __init__(self, duct: Duct, ID: str = ''):
        super().__init__(ID)
        self.duct = duct

    @property
    def zeta(self) -> Quantity:
        return 1.0

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct.velocity_pressure.to('Pa')
        return self.zeta * Pv


class ObstructionA15A(AbstractFitting):
    description = "damper, butterfly, thin plate, round"
    _theta = np.array([0, 10, 20, 30, 40, 50, 60])
    _zeta = np.array([0.20, 0.52, 1.5, 4.5, 11, 29, 108])

    def __init__(self, duct: Duct, theta: Quantity, ID: str = ''):
        super().__init__(ID)
        self.duct = duct
        self.theta = theta.to('deg').magnitude
        self.interp_zeta = interpolate.interp1d(self._theta, self._zeta, bounds_error=False, fill_value='extrapolate')
        self.interp_theta = interpolate.interp1d(self._zeta, self._theta, bounds_error=False, fill_value='extrapolate')

    @property
    def zeta(self) -> float:
        zeta = self.interp_zeta(self.theta)
        return zeta

    @zeta.setter
    def zeta(self, value: float):
        self.theta = self.interp_theta(value)

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct.velocity_pressure.to('Pa')
        return self.zeta * Pv
