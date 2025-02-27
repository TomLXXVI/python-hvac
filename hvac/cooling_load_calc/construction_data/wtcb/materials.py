"""
01. BUILDING MATERIALS

In function `main()` building materials are created and stored on the material
shelf (see '../shelves/material_shelf.py').
"""
import pandas as pd
from hvac.cooling_load_calc.construction_data.shelves import MaterialShelf
from hvac.cooling_load_calc.construction_data.wtcb.setup import db_path


def main(show: bool = True):
    MaterialShelf.add(
        ('gypsum-plaster', 0.56, 1300.0, 840.0),
        ('terracotta-block-1200kg/m3', 0.51, 1200.0, 840.0),
        ('terracotta-brick-1500kg/m3', 1.59, 1500.0, 840.0),
        ('mineral-wool-glass-fibre-cover-sheet-22kg/m3', 0.04, 22.0, 1030.0),
        ('concrete-aerated-block-glued-600kg/m3', 0.22, 600.0, 840.0),
        ('expanded-clay-solid-block', 0.57, 1400.0, 840.0),
        ('concrete-medium-solid-block-1700kg/m3', 1.19, 1700.0, 1000.0),
        ('expanded-clay-hollow-block-1200kg/m3-t=14cm', None, 1200.0, 1000.0, 0.30),
        ('polyurethane-sheet-150kg/m3', 0.04, 150.0, 1470.0),
        ('expanded-clay-hollow-block-1200kg/m3-t=19cm', None, 1200.0, 1000.0, 0.35),
        ('concrete-heavy-hollow-block-1400kg/m3-t=14cm', None, 1200.0, 1000.0, 0.11),
        ('concrete-medium-solid-block-1800kg/m3', 1.69, 1800.0, 1000.0),
        ('gypsum-cardboard', 0.21, 700.0, 1000.0),
        ('osb-board', 0.13, 650.0, 1700.0),
        ('cement-plaster', 1.55, 1900.0, 840.0),
        ('terracotta-block-1500kg/m3', 0.51, 1500.0, 840.0),
        ('polystyrene-extruded-sheet', 0.04, 30.0, 1500.0),
        ('roof-tile-clay', 1.0, 2000.0, 800.0),
        ('plywood', 0.17, 1000.0, 1680.0),
        ('roof-tile-ceramic', 1.3, 2300.0, 840.0),
        ('mineral-wool-sheet-uncovered-135kg/m3', 0.04, 135.0, 840.0),
        ('wood-pine', 0.16, 550.0, 1630.0),
        ('bitumen', 0.17, 1050.0, 1000.0),
        ('epdm', 0.25, 1150.0, 1000.0),
        ('polystyrene-expanded-sheet', 0.04, 30.0, 1500.0),
        ('concrete-light-1600kg/m3', 1.20, 1600.0, 840.0),
        ('concrete-light-1900kg/m3', 1.40, 1900.0, 840.0),
        ('gravel-concrete-compacted-reinforced', 2.0, 2500.0, 840.0),
        ('gravel-concrete-compacted', 1.7, 2400.0, 840.0),
        ('concrete-reinforced-1%-steel', 2.3, 2300.0, 1000.0),
        ('concrete-reinforced-2%-steel', 2.5, 2400.0, 1000.0),
        ('precast-slab-heavy-concrete-t=12cm', None, 2300.0, 840.0, 0.11),
        ('concrete-reinforced-2400kg/m3', 2.50, 2400.0, 840.0)
    )
    if show:
        with pd.option_context(
            'display.max_rows', None,
            'display.max_columns', None,
            'display.max_colwidth', None,
            'display.colheader_justify', 'center',
            'display.width', None
        ):
            print(MaterialShelf.overview(detailed=True))
    MaterialShelf.export_to_excel(str(db_path / 'materials.ods'))


if __name__ == '__main__':
    main()
