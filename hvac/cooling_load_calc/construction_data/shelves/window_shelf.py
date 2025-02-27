"""
WINDOW THERMAL PROPERTIES SHELF

Implements the class `WindowPropertiesShelf` that encapsulates a Python
shelf to store `WindowThermalProperties` objects on disk.

Notes
-----
Before adding `WindowThermalProperties` objects on the shelf, the path to the
shelf file must be assigned to the class variable `path` of class
`MaterialShelf`.
"""
import shelve
import typing
import pandas as pd
from ....cooling_load_calc import WindowThermalProperties


class WindowPropertiesShelf:
    path: str
    units: dict[str, str] = {
        'U': 'W / (m ** 2 * K)'
    }

    @classmethod
    def add(cls, *wnd_props: WindowThermalProperties) -> None:
        """Adds a single or multiple `WindowThermalProperties` objects to the
        shelf.
        """
        with shelve.open(cls.path) as shelf:
            shelf.update({wp.ID: wp for wp in wnd_props})

    @classmethod
    def load(cls, ID: str) -> WindowThermalProperties:
        """Loads the `WindowThermalProperties` object with the given ID from the
        shelf. Raises a `KeyError` exception if the given ID cannot be found on 
        the shelf.
        """
        with shelve.open(cls.path) as shelf:
            try:
                wnd_props = shelf[ID]
            except KeyError:
                raise KeyError(f"window-type '{ID}' not found") from None
            else:
                return typing.cast(WindowThermalProperties, wnd_props)

    @classmethod
    def overview(
        cls,
        detailed: bool = False,
        do_sort: bool = False
    ) -> list[str] | pd.DataFrame:
        """If parameter `detailed` is `False`, returns a list with the IDs of all
        `WindowThermalProperties` objects that are stored on the shelf.
        Otherwise, returns a Pandas Dataframe with: the U-value of the entire
        window, the center-of-glass SHGC for direct radiation at normal
        incidence, the center-of-glass SHGC for diffuse radiation, and the SHGC
        of the entire window at normal incidence.
        If `do_sort` is True, the IDs are sorted alphabetically.
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
                    'U': [],
                    'SHGC_dir(0°)': [],
                    'SHGC_dir(40°)': [],
                    'SHGC_dir(50°)': [],
                    'SHGC_dir(60°)': [],
                    'SHGC_dir(70°)': [],
                    'SHGC_dir(80°)': [],
                    'SHGC_dif': [],
                    'SHGC_wnd(0°)': []
                }
                max_ID_len = 0
                if not do_sort:
                    items_view = shelf.items()
                else:
                    items_view = sorted(
                        shelf.items(),
                        key=lambda item: item[0].lower()
                    )
                for ID, wnd_props in items_view:
                    if len(ID) > max_ID_len: max_ID_len = len(ID)
                    d['ID'].append(ID)
                    d['U'].append(wnd_props.U.to(cls.units['U']).m)
                    d['SHGC_dir(0°)'].append(wnd_props.SHGC.cog_dir(0))
                    d['SHGC_dir(40°)'].append(wnd_props.SHGC.cog_dir(40))
                    d['SHGC_dir(50°)'].append(wnd_props.SHGC.cog_dir(50))
                    d['SHGC_dir(60°)'].append(wnd_props.SHGC.cog_dir(60))
                    d['SHGC_dir(70°)'].append(wnd_props.SHGC.cog_dir(70))
                    d['SHGC_dir(80°)'].append(wnd_props.SHGC.cog_dir(80))
                    d['SHGC_dif'].append(wnd_props.SHGC.cog_dif)
                    d['SHGC_wnd(0°)'].append(wnd_props.SHGC.wnd)
                df = pd.DataFrame(d)
                # Left align the string in column 'ID' (by default,
                # strings are also right-aligned in Pandas DataFrame)
                df['ID'] = df['ID'].str.ljust(max_ID_len)
                return df

    @classmethod
    def delete(cls, *IDs: str) -> tuple[str] | None:
        """Removes the `WindowThermalProperties` objects with the given IDs from
        the shelf.
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
    def search(cls, search_str: str) -> dict[str, WindowThermalProperties] | None:
        """Searches for `WindowThermalProperties` objects on the shelf that have
        `search_str` in their ID.
        Returns a dictionary with all `WindowThermalProperties` objects that
        have `search_str` in their ID. The dictionary keys are the IDs that each
        map to their corresponding `WindowThermalProperties` object. If no
        `WindowThermalProperties` objects with `search_str` in their ID are
        found, `None` is returned.
        """
        with shelve.open(cls.path) as shelf:
            results = {k: m for k, m in shelf.items() if search_str in k}
        return results or None

    @classmethod
    def replace(cls, ID: str, wnd_props_new: WindowThermalProperties) -> None:
        """Replaces the `WindowThermalProperties` object with the given ID with
        another `WindowThermalProperties` object, which may also have a
        different ID. This method combines in a single call the removal of the
        `WindowThermalProperties` object with the given ID and the addition of
        another `WindowThermalProperties` object on the shelf.
        """
        with shelve.open(cls.path) as shelf:
            if ID in shelf.keys():
                del shelf[ID]
                shelf[wnd_props_new.ID] = wnd_props_new
            else:
                raise KeyError(f"window-type '{ID}' not found") from None

    @classmethod
    def export_to_excel(cls, file_path: str) -> None:
        """Exports the detailed overview of the shelf to an Excel file."""
        df = cls.overview(detailed=True)
        df.to_excel(file_path)
