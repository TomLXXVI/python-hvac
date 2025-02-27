"""BUILDING MATERIALS

Implements class `MaterialShelf` that encapsulates a Python shelf to store
`Material` objects on disk.

Notes
-----
Before adding materials on the shelf, the path to the shelf file must be
set with class variable `path` of class `MaterialShelf`.
"""
import shelve
import typing
import pandas as pd
from hvac import Quantity, UNITS
from ....cooling_load_calc import Material


MaterialRecord = ('MaterialTuple', ('ID', 'k', 'rho', 'c', 'R'))


class MaterialShelf:
    path: str
    units = {
        'k': 'W / (K * m)',
        'rho': 'kg / m**3',
        'c': 'J / (kg * K)',
        'R': 'K * m**2 / W'
    }

    @classmethod
    def add(cls, *records: MaterialRecord) -> None:
        """Adds a single or multiple new materials to the shelf.

        Parameters
        ----------
        records:
            One or more tuples with elements in the order: (ID, k, rho, c, R)
            where:
            - ID: name of the material (str)
            - k: thermal conductivity of the material (float | Quantity)
            - rho: mass density of material (float | Quantity)
            - c: specific heat of material (float | Quantity)
            - R: unit thermal resistance of material (float | Quantity | None)

        Notes
        -----
        It may happen that the thermal conductivity k of a material is not 
        available, but instead the unit thermal resistance R of the material is 
        specified. In that case, set k (index 1) to None, and enter R at index 4
        (last element of the tuple). In all cases where k is known, R (index 4) 
        can be omitted.

        k, rho, c and R can be assigned either a `Quantity` object with the
        appropriate units, or a float. If a float is given the corresponding
        units will be retrieved from the class attribute `units`, a dictionary
        with keys 'k', 'rho', 'c' and 'R' of which the corresponding values are
        standard SI-units: 
        - thermal conductivity 'W / (m * K)', 
        - mass density 'kg / m ** 3', 
        - specific heat 'J / (kg * K)'
        - unit thermal resistance 'm ** 2 * K / W'
        """
        with shelve.open(cls.path) as shelf:
            for record in records:
                ID = record[0]
                material = cls._create_material(*record[1:])
                shelf[ID] = material

    @classmethod
    def _create_material(
        cls,
        k: Quantity | float,
        rho: Quantity | float,
        c: Quantity | float,
        R: Quantity | float | None = None
    ) -> Material:
        """Creates and returns a `Material` object."""
        lst = [k, rho, c, R]
        for i, q, u in zip(range(len(lst)), lst, cls.units.values()):
            if isinstance(q, float):
                lst[i] = Quantity(q, u)
        material = Material(*(e for e in lst))
        return material

    @classmethod
    def load(cls, ID: str) -> Material:
        """Loads the `Material` object with the given ID from the shelf.
        Raises a `KeyError` exception if the given ID cannot be found on the
        shelf.
        """
        with shelve.open(cls.path) as shelf:
            try:
                material = shelf[ID]
            except KeyError:
                raise KeyError(f"Material '{ID}' could not be found.") from None
            else:
                return typing.cast(Material, material)

    @classmethod
    def overview(
        cls,
        detailed: bool = False,
        do_sort: bool = False
    ) -> list[str] | pd.DataFrame:
        """Returns a list with the IDs of all materials that are stored on the
        shelf if parameter `detailed` is `False`.
        Returns a Pandas Dataframe with the thermophysical properties of the
        materials, expressed in units taken from the class attribute `units`,
        if `detailed` is True.
        If `do_sort` is True, the materials are listed alphabetically.
        """
        with shelve.open(cls.path) as shelf:
            if not detailed:
                if not do_sort:
                    return list(shelf.keys())
                else:
                    return sorted(
                        list(shelf.keys()),
                        key=lambda item: item.lower()
                    )
            else:
                d = {
                    'ID': [],
                    'k': [],
                    'rho': [],
                    'c': [],
                    'R': []
                }
                max_ID_len = 0
                if not do_sort:
                    items = shelf.items()
                else:
                    items = sorted(
                        shelf.items(),
                        key=lambda item: item[0].lower()
                    )
                for ID, material in items:
                    if len(ID) > max_ID_len: max_ID_len = len(ID)
                    d['ID'].append(ID)
                    try:
                        d['k'].append(material.k.to(cls.units['k']).m)
                    except AttributeError:
                        d['k'].append(float('nan'))
                    d['rho'].append(material.rho.to(cls.units['rho']).m)
                    d['c'].append(material.c.to(cls.units['c']).m)
                    try:
                        d['R'].append(material.R.to(cls.units['R']).m)
                    except AttributeError:
                        d['R'].append(float('nan'))
                df = pd.DataFrame(d)
                df.columns = [
                    "ID",
                    f"k [{UNITS.Unit(cls.units['k']):~P}]",
                    f"rho [{UNITS.Unit(cls.units['rho']):~P}]",
                    f"c [{UNITS.Unit(cls.units['c']):~P}]",
                    f"R [{UNITS.Unit(cls.units['R']):~P}]"
                ]
                # Left align the string in column 'ID' (by default,
                # strings are also right-aligned in Pandas DataFrame):
                df['ID'] = df['ID'].str.ljust(max_ID_len)
                return df

    @classmethod
    def delete(cls, *IDs: str) -> tuple[str] | None:
        """Removes the material(s) with the given IDs from the shelf.
        If one or more IDs cannot be found, these IDs are returned, otherwise 
        None is returned.
        """
        IDs_not_found = []
        with shelve.open(cls.path) as shelf:
            for ID in IDs:
                try:
                    del shelf[ID]
                except KeyError:
                    IDs_not_found.append(ID)
        return tuple(IDs_not_found) or None

    @classmethod
    def search(cls, search_str: str) -> dict[str, Material] | None:
        """Searches for materials on the shelf that have `search_str` in their
        ID.
        Returns a dictionary with all the materials that have `search_str` in
        their ID. The dictionary keys are the material IDs that each map to
        their corresponding `Material` object. If no materials with `search_str`
        in their ID are found, `None` is returned.
        """
        with shelve.open(cls.path) as shelf:
            results = {k: m for k, m in shelf.items() if search_str in k}
        return results or None

    @classmethod
    def replace(cls, ID: str, record: MaterialRecord) -> None:
        """Replaces the material with the given ID by another material, which
        may also have a different ID. This method combines in a single call the
        removal of the material with the given ID and the addition of another
        material on the shelf.
        """
        ID_new = record[0]
        with shelve.open(cls.path) as shelf:
            try:
                material = typing.cast(Material, shelf[ID])
            except KeyError:
                raise KeyError(f"material '{ID}' not found") from None
            else:
                k_new = record[1] if record[1] is not None else material.k
                rho_new = record[2] if record[2] is not None else material.rho
                c_new = record[3] if record[3] is not None else material.c
                material_new = cls._create_material(k_new, rho_new, c_new)
                del shelf[ID]
                shelf[ID_new] = material_new

    @classmethod
    def export_to_excel(cls, file_path: str) -> None:
        """Exports the detailed overview of the shelf to an Excel file."""
        df = cls.overview(detailed=True)
        df.to_excel(file_path)
