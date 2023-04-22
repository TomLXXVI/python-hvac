from typing import Union, Dict, Optional, List
import math
import sympy
import pint
from hvac import Quantity


class Variable:

    def __init__(
        self,
        name: str,
        value: Union[float, int, None] = None,
        unit: Optional[str] = None
    ) -> None:
        self.name = name
        self.symbol = sympy.Symbol(name)
        self.base_value: Union[float, int, None] = None
        self.base_unit: Optional[str] = None
        self.default_unit: Optional[str] = None
        self._configure_variable(value, unit)

    def _configure_variable(
        self,
        value: Union[float, int, None],
        unit: Union[str, pint.Unit, None] = None
    ) -> None:
        """Configure the data attributes of `Variable` object depending on
        parameters `value` and `unit`.
        """
        if (isinstance(value, (int, float)) and not math.isnan(value)) and (unit is not None):
            # the variable contains a quantity
            value = Quantity(value, unit)
            self.default_unit = value.units
            # convert quantity to its base unit (all calculations will be done
            # in base units)
            value.ito_base_units()
            self.base_value = value.magnitude
            self.base_unit = value.units
        elif (value is None or math.isnan(value)) and (unit is not None):
            # the variable contains a quantity of which the magnitude is not known
            self.base_value = None
            self.default_unit = unit
            value = Quantity(float('nan'), unit)
            value.ito_base_units()
            self.base_unit = value.units
        else:
            # the variable contains a numeric value that may not be defined
            # yet (= None or NaN)
            if value is None or math.isnan(value):
                self.base_value = None
            else:
                self.base_value = value

    def __repr__(self):
        if self.default_unit:
            return f"{self.name} [{self.value:~P}]"
        else:
            return f"{self.name} [{self.base_value}]"
    
    @property
    def value(self) -> Union[float, int, Quantity, None]:
        if self.default_unit:
            # the variable contains a `Quantity` value
            if self.base_value is not None:
                # the `Variable` object has been assigned a real value
                # (not None): return the variable with its default unit (the
                # unit that was passed when creating the `Variable` object)
                value = Quantity(self.base_value, self.base_unit)
                return value.to(self.default_unit)
            else:
                # the `Variable` object hasn't been assigned a real value
                # (is None): return `Quantity` object with undefined magnitude
                # (= NaN).
                return Quantity(float('nan'), self.default_unit)
        else:
            # the `Variable` object contains only a numeric value (float or int)
            return self.base_value
    
    @value.setter
    def value(self, v: Union[float, int, Quantity, None]):
        if isinstance(v, Quantity):
            value = v.magnitude
            unit = v.units
        else:
            value = v
            unit = None
        self._configure_variable(value, unit)


class Equation:

    def __init__(self, variables: List[Variable], lhs: str):
        """
        Create equation.

        Parameters
        ----------
        variables: List[Variable]
            List of `Variable` objects that make up the equation.
        lhs: str
            Left-hand side of the equation. Right-hand side must always be zero.
            The variable names in the expression must correspond with the
            variable names in the list of `Variable` objects.
        """
        self.variables: Dict[str, Variable] = {v.name: v for v in variables}
        for variable in self.variables.values():
            exec(f"{variable.name} = self.variables['{variable.name}'].symbol")
        self.equation = sympy.Eq(eval(lhs), 0)

    def __setitem__(self, name: str, value: Union[float, int, Quantity, None]):
        """Assign `value` to the equation variable of which the name is `name`."""
        if name in self.variables.keys():
            self.variables[name].value = value

    def solve(self, unknown_variable_name: Optional[str] = None) -> Variable:
        """
        Solve equation for unknown variable with name `unknown_variable_name`.
        If no name is given, try to solve for any unknown variable in the
        equation. If there is more than 1 unknown variable in the equation, a
        `ValueError` is raised to inform that there are too many unknowns.
        If there are no unknown variables and no variable name was given, a
        `ValueError` is also raised to inform that the equation is already
        solved. However, in case a variable name was given, the corresponding
        variable will be returned.
        """
        # get symbols and values of all known variables
        known_vars = [
            (var.symbol, var.base_value)
            for var in self.variables.values()
            if var.base_value is not None
        ]

        # get all unknown variables
        unknown_vars = tuple(
            var for var in self.variables.values()
            if var.base_value is None
        )
        if len(unknown_vars) == 0:
            # equation already solved
            try:
                return self.variables[unknown_variable_name]
            except KeyError:
                raise ValueError('equation already solved') from None
        if len(unknown_vars) > 1:
            # too many unknown variables to solve with a single equation
            raise ValueError('too many unknowns')

        unknown_var = unknown_vars[0]

        # substitute known variables into the equation
        eq = self.equation.subs(known_vars)

        # solve equation for the unknown variable
        solution = sympy.solveset(eq, unknown_var.symbol)

        # assign the solution to the unknown variable
        try:
            # noinspection PyTypeChecker
            unknown_var.base_value = float(list(solution)[0])
        except (TypeError, IndexError):
            raise ValueError('no valid solution found') from None
        else:
            # return the solved unknown variable
            return unknown_var

    def __repr__(self):
        repr_ = (
            f"Equation: {self.equation.lhs} = 0\n"
            f"Variables: {[repr(var) for var in self.variables.values()]}"
        )
        return repr_
