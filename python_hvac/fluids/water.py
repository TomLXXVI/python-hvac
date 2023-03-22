from typing import Dict
from iapws import IAPWS95
from .. import Quantity


Q_ = Quantity


class Water:
    units: Dict[str, str] = {
        'T': 'K',
        'P': 'MPa',
        'rho': 'kg / m ** 3',
        'v': 'm ** 3 / kg',
        'u': 'kJ / kg',
        'h': 'kJ / kg',
        's': 'kJ / kg / K',
        'cp': 'kJ / kg / K',
        'cv': 'kJ / kg / K',
        'x': 'frac',
        'k': 'W / m / K',
        'mu': 'Pa * s'
    }

    def __init__(self, **input_qties: Quantity):
        keys = list(input_qties.keys())
        qties = list(input_qties.values())
        values = [qty.to(self.units[key]).m for qty, key in zip(qties, keys)]
        inputs = {keys[0]: values[0], keys[1]: values[1]}
        self._state = IAPWS95(**inputs)

    @property
    def T(self) -> Quantity:
        return Q_(self._state.T, self.units['T'])

    @property
    def P(self) -> Quantity:
        return Q_(self._state.P, self.units['P'])

    @property
    def rho(self) -> Quantity:
        return Q_(self._state.rho, self.units['rho'])

    @property
    def v(self) -> Quantity:
        return Q_(self._state.v, self.units['v'])

    @property
    def u(self) -> Quantity:
        return Q_(self._state.u, self.units['u'])

    @property
    def h(self) -> Quantity:
        return Q_(self._state.h, self.units['h'])

    @property
    def s(self) -> Quantity:
        return Q_(self._state.s, self.units['s'])

    @property
    def cp(self) -> Quantity:
        return Q_(self._state.cp, self.units['cp'])

    @property
    def cv(self) -> Quantity:
        return Q_(self._state.cv, self.units['cv'])

    @property
    def x(self) -> Quantity:
        return Q_(self._state.x, self.units['x'])

    @property
    def k(self) -> Quantity:
        return Q_(self._state.k, self.units['k'])

    @property
    def mu(self) -> Quantity:
        return Q_(self._state.mu, self.units['mu'])
