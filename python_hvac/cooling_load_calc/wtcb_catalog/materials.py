"""
01. BUILDING MATERIALS

Add the building materials to the materials shelf that will be
needed to create the construction assemblies (see construction_assemblies.py).
"""
import pandas as pd
from hvac.cooling_load_calc.wtcb_catalog import MaterialsShelf, db_path


def main():
    MaterialsShelf.add(
        ('gipspleister', 0.56, 1300.0, 840.0),
        ('blokken gebakken aarde, 1200 kg/m3', 0.51, 1200.0, 840.0),
        ('bakstenen gebakken aarde, 1500 kg/m3', 1.59, 1500.0, 840.0),
        ('minerale wol, plaat + glasvlies, 22 kg/m3', 0.04, 22.0, 1030.0),
        ('blokken cellenbeton, gelijmd, 600 kg/m3', 0.22, 600.0, 840.0),
        ('volle betonblokken, geëxpandeerde klei', 0.57, 1400.0, 840.0),
        ('volle blokken halfzwaar beton, 1700 kg/m3', 1.19, 1700.0, 1000.0),
        ('holle betonblokken, geëxpandeerde klei, 1200 kg/m3, t=14 cm', None, 1200.0, 1000.0, 0.30),
        ('polyurethaan (PUR), plaat, 150 kg/m3', 0.04, 150.0, 1470.0),
        ('holle betonblokken, geëxpandeerde klei, 1200 kg/m3, t=19 cm', None, 1200.0, 1000.0, 0.35),
        ('holle blokken zwaar beton, 1400 kg/m3, t=14 cm', None, 1200.0, 1000.0, 0.11),
        ('volle blokken halfzwaar beton, 1800 kg/m3', 1.69, 1800.0, 1000.0),
        ('volle blokken halfzwaar beton, 1700 kg/m3', 1.19, 1800.0, 1000.0),
        ('gipsplaat, tussen 2 lagen karton', 0.21, 700.0, 1000.0),
        ('OSB-plaat', 0.13, 650.0, 1700.0),
        ('cementpleister', 1.55, 1900.0, 840.0),
        ('blokken gebakken aarde, 1500 kg/m3', 0.51, 1500.0, 840.0),
        ('geëxtrudeerd polystyreen (XPS), plaat', 0.04, 30.0, 1500.0),
        ('kleidakpannen', 1.0, 2000.0, 800.0),
        ('multiplexplaat', 0.17, 1000.0, 1680.0),
        ('keramische dakpannen', 1.3, 2300.0, 840.0),
        ('minerale wol, onbeklede plaat, 135 kg/m3', 0.04, 135.0, 840.0),
        ('naaldhout', 0.16, 550.0, 1630.0),
        ('bitumen', 0.17, 1050.0, 1000.0),
        ('EPDM', 0.25, 1150.0, 1000.0),
        ('geëxpandeerd polystyreen (EPS), plaat', 0.04, 30.0, 1500.0),
        ('licht beton, 1600 kg/m3', 1.20, 1600.0, 840.0),
        ('licht beton, 1900 kg/m3', 1.40, 1900.0, 840.0),
        ('grindbeton, verdicht, gewapend', 2.0, 2500.0, 840.0),
        ('grindbeton, verdicht, ongewapend', 1.7, 2400.0, 840.0),
        ('beton, gewapend (1 % staal)', 2.3, 2300.0, 1000.0),
        ('beton, gewapend (2 % staal)', 2.5, 2400.0, 1000.0),
        ('welfsel, zwaar beton, t=12 cm', None, 2300.0, 840.0, 0.11),
        ('beton, gewapend, 2400 kg/m3', 2.50, 2400.0, 840.0)
    )
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.max_colwidth', None,
        'display.colheader_justify', 'center',
        'display.width', None
    ):
        print(MaterialsShelf.overview(detailed=True))
    MaterialsShelf.export_to_excel(str(db_path / 'materials.ods'))


if __name__ == '__main__':
    main()
