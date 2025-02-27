"""
CONSTRUCTION ASSEMBLY SHELF

Implements the class `ConstructionAssemblyShelf` that encapsulates a Python
shelf to store `ConstructionAssembly` objects on disk.

Notes
-----
Before adding construction assemblies on the shelf, the path to the shelf file
must be assigned to the class variable `path` of class `MaterialShelf`.
"""
import shelve
import typing
import pandas as pd
from hvac import UNITS
from ....cooling_load_calc import ConstructionAssembly


class ConstructionAssemblyShelf:
    path: str
    units: dict[str, str] = {
        't': 'm',
        'R': 'm ** 2 * K / W',
        'U': 'W / (K * m ** 2)',
        'C': 'J / (K * m ** 2)'
    }

    @classmethod
    def add(cls, *constr_assemblies: ConstructionAssembly) -> None:
        """Adds the given construction assemblies to the shelf."""
        with shelve.open(cls.path) as shelf:
            shelf.update({ca.ID: ca for ca in constr_assemblies})

    @classmethod
    def load(cls, ID: str) -> ConstructionAssembly:
        """Loads the construction assembly with the given ID from the shelf.
        Raises a `KeyError` exception if the given ID could not be found on the
        shelf.
        """
        with shelve.open(cls.path) as shelf:
            try:
                ca = shelf[ID]
            except KeyError:
                raise KeyError(
                    f"Construction assembly '{ID}' could not be found"
                ) from None
            else:
                return typing.cast(ConstructionAssembly, ca)

    @classmethod
    def overview(
        cls,
        detailed: bool = False,
        do_sort: bool = False
    ) -> list[str] | pd.DataFrame:
        """Returns a list with the IDs of all the construction assemblies
        stored on the shelf, if parameter `detailed` is `False`.
        Returns a Pandas Dataframe with the thickness, thermal unit resistance,
        and thermal unit conductance (transmittance) expressed in units taken
        from the class attribute `units`.
        If `do_sort` is True, the construction assemblies are listed
        alphabetically.
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
                # Left align the string in column 'ID' (by default,
                # strings are also right-aligned in Pandas DataFrame):
                df['ID'] = df['ID'].str.ljust(max_ID_len)
                return df

    @classmethod
    def delete(cls, *IDs: str) -> tuple[str] | None:
        """Removes the construction assemblies with the given IDs from the shelf.
        If one or more IDs could not be found, these IDs are returned, otherwise
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
    def search(cls, search_str: str) -> dict[str, ConstructionAssembly] | None:
        """Searches for construction assemblies on the shelf that have
        `search_str` in their ID.
        Returns a dictionary with all construction assemblies that have
        `search_str` in their ID. The dictionary keys are the IDs that each map
        to their corresponding `ConstructionAssembly` object. If no construction
        assemblies with `search_str` in their ID were found, `None` is returned.
        """
        with shelve.open(cls.path) as shelf:
            results = {k: m for k, m in shelf.items() if search_str in k}
        return results or None

    @classmethod
    def replace(cls, ID: str, ca_new: ConstructionAssembly) -> None:
        """Replaces the construction assembly with the given ID with another
        construction assembly, which may also have a different ID. This method
        combines in a single call the removal of the construction assembly with
        the given ID and the addition of another construction assembly on the
        shelf.
        """
        with shelve.open(cls.path) as shelf:
            if ID in shelf.keys():
                del shelf[ID]
                shelf[ca_new.ID] = ca_new
            else:
                raise KeyError(f"construction assembly '{ID}' not found") from None

    @classmethod
    def export_to_excel(cls, file_path: str) -> None:
        """Exports the detailed overview of the shelf to an Excel file."""
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
        multi_index_df = pd.DataFrame(
            multi_index_arr,
            columns=['assembly-ID', 'layer-ID']
        )
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
