from hvac import UNITS


class Units:
    unit_T = UNITS.Unit('K')
    unit_t = UNITS.Unit('s')
    unit_E = UNITS.Unit('J')
    unit_R = unit_T * unit_t / unit_E
    unit_C = unit_E / unit_T
    unit_Q = unit_E / unit_t
    _registry = []
    
    @classmethod
    def set_to_defaults(cls) -> None:
        cls.unit_E = UNITS.Unit('J')
        cls.unit_T = UNITS.Unit('K')
        cls.unit_t = UNITS.Unit('s')
        cls._update()
    
    @classmethod
    def set_time_unit(cls, t_unit: str) -> None:
        cls.unit_t = UNITS.Unit(t_unit)
        cls._update()
    
    @classmethod
    def set_temperature_unit(cls, T_unit: str) -> None:
        cls.unit_T = UNITS.Unit(T_unit)
        cls._update()
        
    @classmethod
    def set_energy_unit(cls, E_unit) -> None:
        cls.unit_E = UNITS.Unit(E_unit)
        cls._update()
        
    @classmethod
    def _update(cls) -> None:
        cls.unit_R = cls.unit_T * cls.unit_t / cls.unit_E
        cls.unit_C = cls.unit_E / cls.unit_T
        cls.unit_Q = cls.unit_E / cls.unit_t
        for class_obj in cls._registry:
            class_obj.update()
        
    @classmethod
    def register(cls, *class_objs) -> None:
        cls._registry.extend(class_objs)
        