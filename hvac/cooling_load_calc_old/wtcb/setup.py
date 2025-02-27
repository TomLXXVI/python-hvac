"""
00. WTCB CATALOG INITIALISATION
Sets the file paths to the shelves with materials, construction assemblies and
window properties.
"""
from pathlib import Path

from hvac.cooling_load_calc_old.shelves import (
    MaterialShelf,
    ConstructionAssemblyShelf,
    WindowPropertiesShelf
)

# The shelves are stored in the user's home directory in a directory
# 'wtcb-database'.
db_path = Path.home() / 'wtcb-database'
db_path.mkdir(exist_ok=True)

# Assign the path to the shelf file to their respective class.
MaterialShelf.path = str(db_path / 'wtcb-materials.db')
ConstructionAssemblyShelf.path = str(db_path / 'wtcb-construction-assemblies.db')
WindowPropertiesShelf.path = str(db_path / 'window-properties.db')
