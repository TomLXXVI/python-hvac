import typing
import shelve
import pandas as pd
from hvac import UNITS
from ..core import ConstructionAssembly


class ConstructionAssembliesShelf:
    path: str
    units: dict[str, str] = {
        't': 'm',
        'R': 'm ** 2 * K / W',
        'U': 'W / (K * m ** 2)',
        'C': 'J / (K * m ** 2)'
    }

    @classmethod
    def add(cls, *constr_assemblies: ConstructionAssembly) -> None:
        """
        Adds a single or multiple new construction assemblies to the shelf.
        """
        with shelve.open(cls.path) as shelf:
            shelf.update({ca.ID: ca for ca in constr_assemblies})

    @classmethod
    def load(cls, ID: str) -> ConstructionAssembly:
        """
        Loads the construction assembly with the given ID from the shelf.

        Parameters
        ----------
        ID:
            The name of the construction assembly.

        Returns
        -------
        Instance of class `ConstructionAssembly` if the given ID is on the shelf.

        Raises
        ------
        `KeyError` if the given ID cannot be found on the shelf.
        """
        with shelve.open(cls.path) as shelf:
            try:
                ca = shelf[ID]
            except KeyError:
                raise KeyError(f"construction assembly '{ID}' not found") from None
            else:
                return typing.cast(ConstructionAssembly, ca)

    @classmethod
    def overview(
        cls,
        detailed: bool = False,
        do_sort: bool = False
    ) -> list[str] | pd.DataFrame:
        """
        If parameter `detailed` is `False`, returns a list with the IDs of all
        construction assemblies that are stored on the shelf. Otherwise, returns
        a Pandas Dataframe with the thickness, unit resistance, and unit
        conductance expressed in the units taken from class attribute `units`.
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
                    't': [],
                    'R': [],
                    'U': []
                }
                max_ID_len = 0
                if not do_sort:
                    items_view = shelf.items()
                else:
                    items_view = sorted(
                        shelf.items(),
                        key=lambda item: item[0].lower()
                    )
                for ID, ca in items_view:
                    if len(ID) > max_ID_len: max_ID_len = len(ID)
                    d['ID'].append(ID)
                    d['t'].append(ca.thickness.to(cls.units['t']).m)
                    d['R'].append(ca.R.to(cls.units['R']).m)
                    d['U'].append(ca.U.to(cls.units['U']).m)
                df = pd.DataFrame(d)
                df.columns = [
                    'ID',
                    f"t [{UNITS.Unit(cls.units['t']):~P}]",
                    f"R [{UNITS.Unit(cls.units['R']):~P}]",
                    f"U [{UNITS.Unit(cls.units['U']):~P}]"
                ]
                # left align the string in column 'ID' (by default,
                # strings are also right-aligned in Pandas DataFrame)
                df['ID'] = df['ID'].str.ljust(max_ID_len)
                return df

    @classmethod
    def delete(cls, *IDs: str) -> tuple[str] | None:
        """
        Removes the construction assemblies indicated by *IDs from the shelf.

        Returns
        -------
        If one or more IDs have not been found, this ID or IDs are returned,
        otherwise None is returned.
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
    def search(cls, search_str: str) -> dict[str, ConstructionAssembly] | None:
        """
        Searches for construction assemblies on the shelf that have `search_str`
        in their ID.

        Returns
        -------
        A dictionary with all construction assemblies that have `search_str` in
        their ID. The dictionary keys are the IDs that each map to their
        corresponding `ConstructionAssembly` object. If no construction assembles
        with `search_str` in their ID were found, `None` is returned.
        """
        with shelve.open(cls.path) as shelf:
            results = {k: m for k, m in shelf.items() if search_str in k}
        return results or None

    @classmethod
    def replace(cls, ID: str, ca_new: ConstructionAssembly) -> None:
        """
        Replaces the construction assembly with given ID with another construction
        assembly, which may also have a different ID. This method combines the
        removal of the construction assembly with the given ID and adding
        another construction assembly on the shelf.
        """
        with shelve.open(cls.path) as shelf:
            if ID in shelf.keys():
                del shelf[ID]
                shelf[ca_new.ID] = ca_new
            else:
                raise KeyError(f"construction assembly '{ID}' not found") from None

    @classmethod
    def export_to_excel(cls, file_path: str) -> None:
        multi_index_arr = []
        data = []
        with shelve.open(cls.path) as shelf:
            for ca in shelf.values():
                for layer in ca.layers.values():
                    multi_index_arr.append([ca.ID, layer.ID])
                    data.append([
                        layer.geometry.t.to(cls.units['t']).m,
                        layer.R.to(cls.units['R']).m,
                        layer.U.to(cls.units['U']).m,
                        layer.C.to(cls.units['C']).m
                    ])
                multi_index_arr.append([ca.ID, 'global'])
                data.append([
                    ca.thickness.to(cls.units['t']).m,
                    ca.R.to(cls.units['R']).m,
                    ca.U.to(cls.units['U']).m,
                    float('nan')
                ])
        multi_index_df = pd.DataFrame(multi_index_arr, columns=['assembly-ID', 'layer-ID'])
        multi_index = pd.MultiIndex.from_frame(multi_index_df)
        df = pd.DataFrame(
            data=data,
            index=multi_index,
            columns=[
                f"t [{UNITS.Unit(cls.units['t']):~P}]",
                f"R [{UNITS.Unit(cls.units['R']):~P}]",
                f"U [{UNITS.Unit(cls.units['U']):~P}]",
                f"C [{UNITS.Unit(cls.units['C']):~P}]"
            ]
        )
        df.to_excel(file_path)
