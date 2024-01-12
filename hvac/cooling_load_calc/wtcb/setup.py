"""
00. WTCB CATALOG INITIALISATION
Set the file paths to the material, construction assembly and window properties
shelves.
"""
from pathlib import Path
from hvac.cooling_load_calc.shelves import (
    MaterialShelf,
    ConstructionAssemblyShelf,
    WindowPropertiesShelf
)

# The shelves are stored in the user's home directory in a directory called
# "wtcb-database":
db_path = Path.home() / 'wtcb-database'
db_path.mkdir(exist_ok=True)

MaterialShelf.path = str(db_path / 'wtcb-materials.db')
ConstructionAssemblyShelf.path = str(db_path / 'wtcb-construction-assemblies.db')
WindowPropertiesShelf.path = str(db_path / 'window-properties.db')
