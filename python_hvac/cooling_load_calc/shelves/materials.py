import typing
import shelve
import pandas as pd
from hvac import Quantity, UNITS
from ..core import Material


MaterialTuple = tuple[
    str,
    Quantity | float,
    Quantity | float,
    Quantity | float,
    Quantity | float | None
]


class MaterialsShelf:
    path: str
    units = {
        'k': 'W / (m * K)',
        'rho': 'kg / m ** 3',
        'c': 'J / (kg * K)',
        'R': 'm ** 2 * K / W'
    }

    @classmethod
    def add(cls, *records: MaterialTuple) -> None:
        """
        Adds a single or multiple new materials to the shelf.

        Parameters
        ----------
        records:
            One or more tuples with elements in the exact order (ID, k, rho, c, R)
            where:
            - ID (index 0): name of the material [str]
            - k (index 1): thermal conductivity of the material [float | Quantity]
            - rho (index 2): mass density of material [float | Quantity]
            - c (index 3): specific heat of material [float | Quantity]
            - R (index 4): unit thermal resistance of material [float | Quantity | None]

        Notes
        -----
        It may happen that k is not available, but instead the thermal unit
        resistance R of the material is specified. In such a case, set k (index 1)
        to None, and enter R at index 4 (last element of the tuple).
        In all cases where k is known, R (index 4) can be omitted.

        k, rho, c and R can be given either a `Quantity`-object with appropriate
        units, or a float. In case a float is given the corresponding units will
        be retrieved from the class attribute `units`, being a dictionary with
        keys 'k', 'rho', 'c' and 'R' of which the corresponding values are the
        standard SI-units for thermal conductivity ('W / (m * K)'), mass density
        ('kg / m ** 3'), specific heat ('J / (kg * K)'), and unit thermal
        resistance ('m ** 2 * K / W').

        Returns
        -------
        None
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
        """
        Creates and returns a `Material` object.
        """
        lst = [k, rho, c, R]
        for i, q, u in zip(range(len(lst)), lst, cls.units.values()):
            if isinstance(q, float):
                lst[i] = Quantity(q, u)
        material = Material(*(e for e in lst))
        return material

    @classmethod
    def load(cls, ID: str) -> Material:
        """
        Loads the material with the given ID from the shelf.

        Parameters
        ----------
        ID:
            The name of the material

        Returns
        -------
        Instance of class `Material` if the given ID is on the shelf.

        Raises
        ------
        `KeyError` if the given ID cannot be found on the shelf.
        """
        with shelve.open(cls.path) as shelf:
            try:
                material = shelf[ID]
            except KeyError:
                raise KeyError(f"material '{ID}' not found") from None
            else:
                return typing.cast(Material, material)

    @classmethod
    def overview(
        cls,
        detailed: bool = False,
        do_sort: bool = False
    ) -> list[str] | pd.DataFrame:
        """
        If parameter `detailed` is `False`, returns a list with the IDs of all
        materials that are stored on the shelf. Otherwise, returns a Pandas
        Dataframe with the thermophysical properties of the materials, expressed
        in the units taken from class attribute `units` (see also method `add`).
        If `do_sort` is True, the materials are sorted alphabetically.
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
                    items_view = shelf.items()
                else:
                    items_view = sorted(
                        shelf.items(),
                        key=lambda item: item[0].lower()
                    )
                for ID, material in items_view:
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
                    'ID',
                    f"k [{UNITS.Unit(cls.units['k']):~P}]",
                    f"rho [{UNITS.Unit(cls.units['rho']):~P}]",
                    f"c [{UNITS.Unit(cls.units['c']):~P}]",
                    f"R [{UNITS.Unit(cls.units['R']):~P}]"
                ]
                # left align the string in column 'material-ID' (by default,
                # strings are also right-aligned in Pandas DataFrame)
                df['ID'] = df['ID'].str.ljust(max_ID_len)
                return df

    @classmethod
    def delete(cls, *IDs: str) -> tuple[str] | None:
        """
        Removes the material(s) indicated by *IDs from the shelf.

        Returns
        -------
        If one or more material IDs have not been found, this ID or IDs are
        returned, otherwise None is returned.
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
        """
        Searches for materials on the shelf that have `search_str` in their ID.

        Returns
        -------
        A dictionary with all materials that have `search_str` in their ID. The
        dictionary keys are the material IDs that each map to their corresponding
        `Material` object. If no materials with `search_str` in their ID were
        found, `None` is returned.
        """
        with shelve.open(cls.path) as shelf:
            results = {k: m for k, m in shelf.items() if search_str in k}
        return results or None

    @classmethod
    def replace(cls, ID: str, record: MaterialTuple) -> None:
        """
        Replaces the material with given ID with another material, which may
        also have a different ID. This method combines the removal of the
        material with the given ID and adding another material on the shelf.
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
        df = cls.overview(detailed=True)
        df.to_excel(file_path)
