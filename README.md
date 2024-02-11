# xyz2atoms / 2024 python homework
Generates sequential numbering to atoms within a nucleotide and creates a plot using Plotly to visualize names/bond length/indexes.
Built upon xyz2graph.
![](img/example.png)

## Installation
```
python -m pip install git+https://github.com/MelnychenkoM/xyz2atoms.git
```
## Requirements
- numpy
- pandas
- plotly
## Usage
```
# Import module
from xyz2atoms import Molecule
from xyz2atoms import to_plotly_figure, torsion_angle

# Create molecule object and load xyz file.
mol = Molecule()
mol.read_xyz('example/20_rf_xyz_dftV3#r15_084_GbbbB.xyz')

# Create plotly plot.
to_plotly_figure(mol)

# Load pdb file with different models
pdb_mol = Molecule()
pdb_mol.read_pdb('example/example.pdb', model=1)
to_plotly_figure(pdb_mol)

# Access by index
mol[1]

# Access by atom name
mol["C4'"]

# Calculate torsion angle
torsion_angle(mol["C1'"], mol["C2'"], mol["C3'"], mol["C4'"])
```