"""00. WTCB CATALOG INITIALISATION

Sets the file paths to the shelves with materials, construction assemblies and
window properties. The shelves are stored in the user's home directory in a 
directory 'wtcb-database'.
"""
from pathlib import Path

from ..shelves import (
    MaterialShelf,
    ConstructionAssemblyShelf,
    WindowPropertiesShelf
)

# The shelves are stored in the user's home directory in a directory
# 'wtcb-database'.
db_path = Path.home() / 'wtcb-database'
db_path.mkdir(exist_ok=True)

# Assign the path to the shelf file to their respective class.
MaterialShelf.path = str(db_path / 'materials')
ConstructionAssemblyShelf.path = str(db_path / 'construction-assemblies')
WindowPropertiesShelf.path = str(db_path / 'window-properties')
