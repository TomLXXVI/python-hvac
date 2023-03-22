"""
00. WTCB CATALOG INITIALISATION

Set the paths to the materials and construction assemblies shelves.
"""
from pathlib import Path
from hvac.cooling_load_calc.shelves import (
    MaterialsShelf,
    ConstructionAssembliesShelf,
    WindowPropertiesShelf
)

# the shelves are stored in the user's home directory in the directory
# wtcb_database
db_path = Path.home() / 'wtcb_database'
db_path.mkdir(exist_ok=True)

MaterialsShelf.path = str(db_path / 'wtcb_materials.db')
ConstructionAssembliesShelf.path = str(db_path / 'wtcb_construction_assemblies.db')
WindowPropertiesShelf.path = str(db_path / 'window_properties.db')
