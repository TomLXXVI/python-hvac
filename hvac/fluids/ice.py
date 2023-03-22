from typing import Dict

# noinspection PyProtectedMember
from iapws import _Ice

from .. import Quantity


Q_ = Quantity


class Ice:
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

    def __init__(self, P: Quantity, T: Quantity):
        keys = ['P', 'T']
        values = [P.to('MPa').m, T.to('K').m]
        inputs = {keys[0]: values[0], keys[1]: values[1]}
        self._state = _Ice(**inputs)

    @property
    def T(self) -> Quantity:
        return Q_(self._state['T'], self.units['T'])

    @property
    def P(self) -> Quantity:
        return Q_(self._state['P'], self.units['P'])

    @property
    def rho(self) -> Quantity:
        return Q_(self._state['rho'], self.units['rho'])

    @property
    def v(self) -> Quantity:
        return 1 / self.rho

    @property
    def u(self) -> Quantity:
        return Q_(self._state['u'], self.units['u'])

    @property
    def h(self) -> Quantity:
        return Q_(self._state['h'], self.units['h'])

    @property
    def s(self) -> Quantity:
        return Q_(self._state['s'], self.units['s'])

    @property
    def cp(self) -> Quantity:
        return Q_(self._state['cp'], self.units['cp'])
