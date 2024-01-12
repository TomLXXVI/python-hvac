"""
03. WINDOW PROPERTIES
Creates `WindowProperties` objects and stores them on the window properties
shelf.

Window properties can be found in ASHRAE Fundamentals 2017, Ch. 15, table 4 and
table 10.
"""
import pandas as pd
from hvac import Quantity
from hvac.cooling_load_calc.core import WindowThermalProperties
from hvac.cooling_load_calc.wtcb.setup import WindowPropertiesShelf, db_path

Q_ = Quantity


def create_window_5a() -> WindowThermalProperties:
    # Operable window, wood/vinyl with uncoated double glazing, CLR-CLR
    wnd_props = WindowThermalProperties(
        ID='window-5a-operable-wood/vinyl',
        U=Q_(2.86, 'W / (m ** 2 * K)'),
        SHGC_cog_dir={
            0.00: 0.76,
            40.0: 0.74,
            50.0: 0.71,
            60.0: 0.64,
            70.0: 0.50,
            80.0: 0.26
        },
        SHGC_cog_dif=0.66,
        SHGC_wnd=0.62
    )
    return wnd_props


def main():
    wnd_5a = create_window_5a()

    WindowPropertiesShelf.add(wnd_5a)

    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.max_colwidth', None,
        'display.colheader_justify', 'center',
        'display.width', None
    ):
        print(WindowPropertiesShelf.overview(detailed=True))

    WindowPropertiesShelf.export_to_excel(str(db_path / 'windows.ods'))


if __name__ == '__main__':
    main()
