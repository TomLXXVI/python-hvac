import typing
import shelve
import pandas as pd
from ..core import WindowThermalProperties


class WindowPropertiesShelf:
    path: str
    units: dict[str, str] = {
        'U': 'W / (m ** 2 * K)'
    }

    @classmethod
    def add(cls, *wnd_props: WindowThermalProperties) -> None:
        """
        Adds a single or multiple new window-types to the shelf.
        """
        with shelve.open(cls.path) as shelf:
            shelf.update({wp.ID: wp for wp in wnd_props})

    @classmethod
    def load(cls, ID: str) -> WindowThermalProperties:
        """
        Loads the window-type with the given ID from the shelf.

        Parameters
        ----------
        ID:
            The name of the window-type.

        Returns
        -------
        Instance of class `WindowThermalProperties` if the given ID is on the
        shelf.

        Raises
        ------
        `KeyError` if the given ID cannot be found on the shelf.
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
        """
        If parameter `detailed` is `False`, returns a list with the IDs of all
        window-types that are stored on the shelf. Otherwise, returns
        a Pandas Dataframe with the U-value of the entire window, the
        center-of-glass SHGC for direct radiation at normal incidence, the
        center-of-glass SHGC for diffuse radiation, and the SHGC of the entire
        window at normal incidence.
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
                # left align the string in column 'ID' (by default,
                # strings are also right-aligned in Pandas DataFrame)
                df['ID'] = df['ID'].str.ljust(max_ID_len)
                return df

    @classmethod
    def delete(cls, *IDs: str) -> tuple[str] | None:
        """
        Removes the window-types indicated by *IDs from the shelf.

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
    def search(cls, search_str: str) -> dict[str, WindowThermalProperties] | None:
        """
        Searches for window-types on the shelf that have `search_str`
        in their ID.

        Returns
        -------
        A dictionary with all window-types that have `search_str` in
        their ID. The dictionary keys are the IDs that each map to their
        corresponding `WindowThermalProperties` object. If no window-types
        with `search_str` in their ID were found, `None` is returned.
        """
        with shelve.open(cls.path) as shelf:
            results = {k: m for k, m in shelf.items() if search_str in k}
        return results or None

    @classmethod
    def replace(cls, ID: str, wnd_props_new: WindowThermalProperties) -> None:
        """
        Replaces the window-type with given ID with another window-type, which
        may also have a different ID. This method combines the removal of the
        window-type with the given ID and adding another window-type on the shelf.
        """
        with shelve.open(cls.path) as shelf:
            if ID in shelf.keys():
                del shelf[ID]
                shelf[wnd_props_new.ID] = wnd_props_new
            else:
                raise KeyError(f"window-type '{ID}' not found") from None

    @classmethod
    def export_to_excel(cls, file_path: str) -> None:
        df = cls.overview(detailed=True)
        df.to_excel(file_path)
