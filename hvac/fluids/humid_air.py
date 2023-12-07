import math
from typing import Dict, Tuple
import warnings
from CoolProp.HumidAirProp import HAPropsSI
from .. import Quantity
from .constants import STANDARD_PRESSURE

Q_ = Quantity


class HumidAir:
    _P: Quantity = STANDARD_PRESSURE.to('Pa').m
    _coolprop_qties: Dict[str, Tuple[str, str]] = {
        'Tdb': ('Tdb', 'K'),
        'Twb': ('Twb', 'K'),
        'Tdp': ('Tdp', 'K'),
        'P': ('P', 'Pa'),
        'Pw': ('P_w', 'Pa'),
        'v': ('V', 'm ** 3 / kg'),
        'W': ('W', 'kg / kg'),
        'RH': ('RH', 'frac'),
        'h': ('H', 'J / kg'),
        'cp': ('cp', 'J / kg / K'),
        'mu': ('mu', 'Pa * s'),
        'k': ('k', 'W / m / K')
    }

    def __init__(self, **input_qties: Quantity):
        P = input_qties.pop('P', None)
        if P is not None: self.P = P
        keys = list(input_qties.keys())
        qties = list(input_qties.values())
        values = [
            qty.to(self._coolprop_qties[key][1]).m
            for qty, key in zip(qties, keys)
        ]
        # Only the first 2 input quantities are considered,
        # should there be more than 2 given:
        self._inputs = {keys[0]: values[0], keys[1]: values[1]}
        # Check the validity of input quantities:
        self._validate_inputs()

    def _validate_inputs(self):
        for k, v in self._inputs.items():
            if math.isnan(v) or v is None:
                raise ValueError(
                    f"Humid air state cannot be determined: "
                    f"parameter {k} is NaN or None."
                )
        RH = self._inputs.get('RH')
        if RH is not None and RH < 0.0:
            warnings.warn(
                message=(
                    "Negative value for RH detected. "
                    "RH has been reset to 0 %."
                ),
                category=RuntimeWarning
            )
            self._inputs['RH'] = 0.0
        W = self._inputs.get('W')
        if W is not None and W < 0.0:
            warnings.warn(
                message=(
                    "Negative value for W detected. "
                    "W has been reset to 0 kg/kg."
                ),
                category=RuntimeWarning
            )
            self._inputs['W'] = 0.0

    def __str__(self):
        return (
            f"{self.Tdb.to('degC'):~P.2f} DB, "
            f"{self.W.to('g/kg'):~P.2f} AH "
            f"({self.RH.to('pct'):~P.0f} RH)"
        )

    @property
    def P(self) -> Quantity:
        return Q_(self._P, 'Pa')

    @P.setter
    def P(self, v: Quantity):
        self._P = v.to('Pa').m

    def _get_quantity(self, key: str) -> Quantity:
        sym_in_key1, sym_in_key2 = tuple(self._inputs.keys())
        sym_in1 = self._coolprop_qties[sym_in_key1][0]
        val_in1 = self._inputs[sym_in_key1]
        sym_in2 = self._coolprop_qties[sym_in_key2][0]
        val_in2 = self._inputs[sym_in_key2]
        sym_out = self._coolprop_qties[key][0]
        try:
            val_out = HAPropsSI(sym_out, 'P', self._P, sym_in1, val_in1, sym_in2, val_in2)
        except ValueError:
            val_out = float('nan')
        unit_out = self._coolprop_qties[key][1]
        qty_out = Q_(val_out, unit_out)
        return qty_out

    @property
    def Tdb(self) -> Quantity:
        return self._get_quantity('Tdb')

    @property
    def W(self) -> Quantity:
        return self._get_quantity('W')

    @property
    def RH(self) -> Quantity:
        return self._get_quantity('RH')

    @property
    def h(self) -> Quantity:
        return self._get_quantity('h')

    @property
    def Twb(self) -> Quantity:
        return self._get_quantity('Twb')

    @property
    def Pw(self) -> Quantity:
        return self._get_quantity('Pw')

    @property
    def Tdp(self) -> Quantity:
        return self._get_quantity('Tdp')

    @property
    def v(self) -> Quantity:
        return self._get_quantity('v')

    @property
    def rho(self) -> Quantity:
        return 1 / self.v

    @property
    def cp(self) -> Quantity:
        return self._get_quantity('cp')

    @property
    def mu(self) -> Quantity:
        return self._get_quantity('mu')

    @property
    def k(self) -> Quantity:
        return self._get_quantity('k')
