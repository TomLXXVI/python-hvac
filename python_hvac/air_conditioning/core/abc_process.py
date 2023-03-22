from typing import Dict, List, Union, Tuple
from abc import ABC
from hvac import Quantity
from hvac.air_conditioning.core.equation import Equation, Variable


class Process(ABC):

    def __init__(self, equations: Dict[str, Equation]):
        self.equations = equations
        self._register: Dict[str, Tuple[List[Variable], List[Equation]]] = {}
        self._create_register()

    def _create_register(self):
        # Create register with references to all variables in the equations with the same name and the equations that
        # contain these variables.
        #
        # Equations can have the same variables. If a variable is solved in one equation, we want this variable also
        # to be updated in the other equations that have that same variable.
        #
        # In Python-code the register looks like:
        #
        #   register = {
        #       'variable_name': ([var1, var2, ...], [eq1, eq2, ...]),
        #       ...
        #   }
        #
        # Visualize the register as a table of which the row indices are the variable names and that has two columns:
        # the first column (with index 0) containing the list of variables that have the same name and the second column
        # (with index 1) the list of equations that contain each of these variables:
        #             +-----------------+------------------+-----------------+
        #             |                 |        0         |         1       |
        #             +-----------------+------------------+-----------------+
        #             |'variable_name1' | [var1, var2,...] | [eq1, eq2, ...] |
        #             +-----------------+------------------+-----------------+
        #             |'variable_name2' | [var1, var2,...] | [eq1, eq2, ...] |
        #             +-----------------+------------------+-----------------+
        #
        for equation in self.equations.values():
            for variable in equation.variables.values():
                if variable.name not in self._register.keys():
                    self._register[variable.name] = ([variable], [equation])
                else:
                    self._register[variable.name][0].append(variable)
                    self._register[variable.name][1].append(equation)

    def __setitem__(self, variable_name: str, variable_value: Union[float, int, Quantity, None]):
        """
        Assign `variable_value` to all variables in the register which have the name `variable_name`.
        Each time a variable is assigned a value, all the equations that contain that variable are tried to be solved
        """
        if variable_name in self._register.keys():
            # update the value of all variables that have `variable_name`
            variables = self._register[variable_name][0]
            for variable in variables:
                variable.value = variable_value

            # try to solve all equations that have this variable
            equations = self._register[variable_name][1]
            for equation in equations:
                try:
                    solved_variable = equation.solve()
                except ValueError:
                    # ignore when a solution cannot be found
                    pass
                else:
                    # if an equation can be solved, update the register again with this solution
                    self[solved_variable.name] = solved_variable.value

    def __getitem__(self, variable_name: str) -> Union[float, int, Quantity, None]:
        """
        Get the value from the variable of which the name is `variable_name`.
        """
        if variable_name in self._register.keys():
            variable = self._register[variable_name][0][0]
            return variable.value

    def solve(self, variable_name: str) -> Union[float, int, Quantity, None]:
        """Solve for `variable_name` with the equations that contain this variable."""
        if variable_name in self._register.keys():
            # get the equations that contain the variable to be solved for
            equations = self._register[variable_name][1]
            # try to find a solution with the available equations
            for equation in equations:
                try:
                    solved_variable = equation.solve(variable_name)
                except ValueError:
                    # skip this equation and continue with the next one (if any)
                    continue
                else:
                    # if there is a solution, update the variable's value in all equations that contain this variable
                    self[variable_name] = solved_variable.value
                    return solved_variable.value
            else:
                # if no solution was found, raise a `ValueError` exception
                raise ValueError(f'no solution found for {variable_name}')

    def show_equations(self):
        """Print the equations in which each process variable is used."""
        output = ""
        for variable_name in self._register.keys():
            title = f"Variable: {variable_name}\n"
            output += title
            output += "-" * len(title) + "\n"
            for equation in self._register[variable_name][1]:
                output += f"{equation}\n\n"
        print(output)
