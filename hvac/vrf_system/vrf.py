from typing import Callable
from enum import Enum
from inspect import signature
from hvac import Quantity
from hvac.logging import ModuleLogger
from .model import VRFModel
from .indoor_unit import IndoorUnit


Q_ = Quantity
logger = ModuleLogger.get_logger(__name__)


class WorkingMode(Enum):
    COOLING = 'cooling'
    HEATING = 'heating'


class VRFSystem:

    def __init__(
        self,
        Q_rated: Quantity,
        W_rated: Quantity,
        model_size: int = 0,
        CR_system: float | None = None,
        Leq_pipe: Quantity = Q_(10, 'm'),
        h_pipe: Quantity = Q_(0, 'm'),
        working_mode: WorkingMode = WorkingMode.COOLING,
        vrf_model_heating: VRFModel | None = None,
        vrf_model_cooling: VRFModel | None = None,
        indoor_units: list[tuple[IndoorUnit, Quantity]] | None = None
    ) -> None:
        """Creates a `VRFSystem` instance.

        Parameters
        ----------
        Q_rated: Quantity
            The available (full-load) cooling/heating capacity of the system at
            rated conditions.
            If `working_mode` is set to `WorkingMode.COOLING`, `Q_rated` should
            be the rated cooling capacity. If `working_mode` is set to
            `WorkingMode.HEATING`, `Q_rated` should be the rated heating capacity.
        W_rated: Quantity
            The input power taken by the VRF-system at rated cooling/heating
            conditions. This value must also be adapted according to the set
            working mode.
        model_size: int
            The model size of the outdoor unit.
        CR_system: float, optional, default None
            The installed combination ratio of the VRF-system. This is the
            ratio of the installed total indoor unit capacity (or the sum of
            the model sizes of the indoor units) to the rated capacity of the
            outdoor unit (or the model size of the outdoor unit).
            If indoor units are added to the `VRFSystem` instance using the
            method `add_indoor_unit`, this parameter won't be used and can
            be left to `None`.
        Leq_pipe: Quantity
            The maximum equivalent length of piping, i.e. the equivalent length
            (which must take the number of bends in the piping trajectory
            into account) between the outdoor unit and the furthest indoor unit.
        h_pipe: Quantity, default None
            The maximum of the heights between the outdoor unit and the indoor
            units or between indoor units.
        working_mode: WorkingMode, default WorkingMode.COOLING
            The working mode of the VRF system: either cooling, or heating.
        vrf_model_heating: VRFModel, default None
            The VRF model that applies to heating mode.
            The VRF model contains the functions that correct available capacity
            and input power depending on indoor air temperature, outdoor air
            temperature, equivalent piping length, combination ratio, part
            load ratio (PLR), and defrost cycli (only for heating). If working
            mode is set to heating, then the `vrf_model_heating` parameter must
            be set to a valid `VRFModel` object and cannot be left `None`.
        vrf_model_cooling: VRFModel, default None
            The VRF model that applies to cooling mode. If working
            mode is set to cooling, then the `vrf_model_cooling` parameter must
            be set to a valid `VRFModel` object and cannot be left `None`.
        indoor_units: list[tuple[IndoorUnit, Quantity]], default None
            List of the indoor units, together with the corresponding indoor air
            setpoint temperature of the room where they are installed. Indoor
            units can also be added after instantiation, using method
            `add_indoor_unit`.
        """
        self.Q_rated = Q_rated
        self.W_rated = W_rated
        self.model_size = model_size
        self.CR_system = CR_system
        self.Leq_pipe = Leq_pipe
        self.h_pipe = h_pipe
        self.working_mode = working_mode
        self.vrf_model_heating = vrf_model_heating
        self.vrf_model_cooling = vrf_model_cooling
        if (
            self.working_mode == WorkingMode.COOLING and
            self.vrf_model_cooling is not None
        ):
            self.vrf_model = self.vrf_model_cooling
        elif (
            self.working_mode == WorkingMode.HEATING and
            self.vrf_model_heating is not None
        ):
            self.vrf_model = self.vrf_model_heating
        else:
            raise TypeError('no valid VRF model attached')
        self.indoor_units: dict[str, tuple[IndoorUnit, Quantity]] = {}
        if indoor_units is not None:
            self.indoor_units = {
                iu.name: (iu, Tia)
                for iu, Tia in indoor_units
            }

    @staticmethod
    def _check_number_of_arguments(f: Callable) -> int:
        sig = signature(f)
        num_args = len(sig.parameters)
        return num_args

    @property
    def Tia_avg(self) -> Quantity:
        """Gets the load-weighted average indoor air temperature of the building.

        Raises
        ------
        This property implies that indoor units are added to the VRF-systems. If
        no indoor units are added, a `ValueError` exception will be raised.
        """
        if self.indoor_units is not None:
            Qiu_tot = sum(
                iu.get_capacity(Tia)
                for iu, Tia in self.indoor_units.values()
            )
            Tia_avg = sum(
                Tia.to('K') * iu.get_capacity(Tia)
                for iu, Tia in self.indoor_units.values()
            ) / Qiu_tot
            return Tia_avg
        else:
            raise ValueError('no indoor units were added to the VRF-system')

    def add_indoor_unit(self, iu: IndoorUnit, Tai: Quantity):
        """Adds a new or replaces an existing indoor unit of the VRF system.

        Parameters
        ----------
        iu:
            Indoor unit
        Tai:
            Air temperature set point in the room where the indoor unit is
            installed.
        """
        self.indoor_units[iu.name] = (iu, Tai)

    def _get_combination_ratio(self) -> float:
        if len(self.indoor_units):
            iu_size = sum(iu.model_size for iu, _ in self.indoor_units.values())
            CR = iu_size / self.model_size
            return CR
        else:
            return self.CR_system

    def get_full_load_outdoor_unit_capacity(
        self,
        Tia_avg: Quantity | None,
        Toa: Quantity
    ) -> Quantity:
        """Get the available cooling/heating capacity of the VRF outdoor unit at
        the given (load-weighted average) indoor air temperature `Tai_avg` and
        outdoor air temperature `Tao`.

        Notes
        -----
        1. If parameter `Tia_avg` is set to None, property `Tia_avg` of this
        class will be used to calculate the load-weighted average indoor air
        temperature based on the indoor air setpoint temperatures that were
        associated to the indoor units of the VRF-system.

        2. In cooling mode, indoor air temperature is wet-bulb, and outdoor air
        temperature is dry-bulb. In heating mode, indoor air temperature is
        dry-bulb, and outdoor air temperature wet-bulb.

        Raises
        ------
        If parameter `Tia_avg` is None and no indoor units were added to the
        VRF-system, a `ValueError` exception will be raised.
        """
        if Tia_avg is None: Tia_avg = self.Tia_avg
        Tia_avg = Tia_avg.to(self.vrf_model.units['temperature']).m
        Toa = Toa.to(self.vrf_model.units['temperature']).m
        
        # Temperature correction
        if self._check_number_of_arguments(self.vrf_model.CAPFT_fun) == 2:
            CAPFT = self.vrf_model.CAPFT_fun(Tia_avg, Toa)
        else:
            CAPFT = self.vrf_model.CAPFT_fun(Toa)

        # Correction if rated combination ratio (ratio of total rated indoor
        # unit capacity or sum of model sizes to outdoor unit capacity or model
        # size) is greater than 100 %.
        CR = self._get_combination_ratio()
        if CR > 1:
            CR_corr = self.vrf_model.CR_corr_fun(CR)
        else:
            CR_corr = 1.0

        # Correction for equivalent piping length and piping height
        Leq_pipe_corr = 1.0
        Leq_pipe = self.Leq_pipe.to(self.vrf_model.units['length']).m
        h_pipe = self.h_pipe.to(self.vrf_model.units['height']).m

        if self.working_mode == WorkingMode.COOLING:
            if self._check_number_of_arguments(self.vrf_model.Leq_pipe_corr_fun) == 3:
                Leq_pipe_corr = self.vrf_model.Leq_pipe_corr_fun(Leq_pipe, CR, h_pipe)
            else:
                Leq_pipe_corr = self.vrf_model.Leq_pipe_corr_fun(Leq_pipe, CR)
        if self.working_mode == WorkingMode.HEATING:
            if self._check_number_of_arguments(self.vrf_model.Leq_pipe_corr_fun) == 2:
                Leq_pipe_corr = self.vrf_model.Leq_pipe_corr_fun(Leq_pipe, h_pipe)
            else:
                Leq_pipe_corr = self.vrf_model.Leq_pipe_corr_fun(Leq_pipe)

        # Correction for defrost cycles if heating mode active
        if (
            self.working_mode == WorkingMode.HEATING and
            self.vrf_model.defrost_corr_fun is not None
        ):
            defrost_corr = self.vrf_model.defrost_corr_fun(Toa)
        else:
            defrost_corr = 1.0

        logger.debug(
            "calculate available capacity with 'Tia_avg' ="
            f" {Tia_avg} and 'Toa' = {Toa}"
        )
        logger.debug(f"CAPFT: {CAPFT}")
        logger.debug(f"CR_corr: {CR_corr}")
        logger.debug(f"Leq_pipe_corr: {Leq_pipe_corr}")

        # Available capacity of outdoor unit
        Q_ou = self.Q_rated * CAPFT * CR_corr * Leq_pipe_corr * defrost_corr
        return Q_ou

    def get_full_load_total_indoor_unit_capacity(self) -> Quantity:
        """Get the maximum deliverable cooling/heating power of all indoor units
        in the VRF system, depending on their associated indoor air setpoint
        temperature.

        Raises
        ------
        `ValueError` exception if no indoor units were added to the VRF-system.
        """
        if self.indoor_units is not None:
            Q_iu = sum(
                iu.get_capacity(Tia)
                for iu, Tia in self.indoor_units.values()
            )
            return Q_iu
        else:
            raise ValueError('no indoor units were added to the VRF-system')

    def get_available_capacity(
        self,
        Toa: Quantity
    ) -> Quantity:
        """Get the available cooling/heating power that the VRF-system can
        deliver at the given operating conditions under full-load. This is the
        smallest value returned from methods `get_full_load_outdoor_unit_capacity`
        and `get_full_load_total_indoor_unit_capacity`.

        Parameters
        ----------
        Toa:
            Outdoor air temperature (wet-bulb when heating, or dry-bulb when
            cooling).

        Raises
        ------
        This method implies that indoor units are added to the VRF-systems. If
        no indoor units are added, a `ValueError` exception will be raised.
        """
        Q_ou_fl = self.get_full_load_outdoor_unit_capacity(self.Tia_avg, Toa)
        Q_iu_fl = self.get_full_load_total_indoor_unit_capacity()
        return min(Q_ou_fl, Q_iu_fl)

    def get_input_power(
        self,
        Tia_avg: Quantity | None,
        Toa: Quantity,
        PLR: Quantity
    ) -> Quantity:
        """Get the input power taken up by VRF system at the given
        (load-weighted average) indoor air temperature `Tia_avg` and outdoor air
        temperature `Toa`, when the part-load ratio equals `PLR`.

        Notes
        -----
        1. If parameter `Tia_avg` is set to None, property `Tia_avg` of this
        class will be used to calculate the load-weighted average indoor air
        temperature based on the indoor air setpoint temperatures that were
        associated to the indoor units of the VRF-system.

        2. In cooling mode, indoor air temperature is wet-bulb, and outdoor air
        temperature is dry-bulb. In heating mode, indoor air temperature is
        dry-bulb, and outdoor air temperature wet-bulb.

        Raises
        ------
        If parameter `Tia_avg` is None and no indoor units were added to the
        VRF-system, a `ValueError` exception will be raised.
        """
        if Tia_avg is None: Tia_avg = self.Tia_avg
        Tia_avg = Tia_avg.to(self.vrf_model.units['temperature']).m
        Toa = Toa.to(self.vrf_model.units['temperature']).m
        PLR = PLR.to('frac').m
        if self._check_number_of_arguments(self.vrf_model.CAPFT_fun) == 2:
            CAPFT = self.vrf_model.CAPFT_fun(Tia_avg, Toa)
        else:
            CAPFT = self.vrf_model.CAPFT_fun(Toa)
        EIRFT = self.vrf_model.EIRFT_fun(Tia_avg, Toa)
        EIRFPLR = self.vrf_model.EIRFPLR_fun(PLR)
        HPRTF = self.vrf_model.HPRTF_fun(PLR)

        logger.debug(
            "calculate input power with 'Tia_avg' ="
            f" {Tia_avg}, 'Toa' = {Toa}, and 'PLR' = {PLR}"
        )
        logger.debug(f"CAPFT: {CAPFT}")
        logger.debug(f"EIRFT: {EIRFT}")
        logger.debug(f"EIRFPLR: {EIRFPLR}")
        logger.debug(f"HPRTF: {HPRTF}")

        W_input = self.W_rated * CAPFT * EIRFT * EIRFPLR * HPRTF
        return W_input
