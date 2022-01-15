from rdkit.Chem import Draw
from rdkit.Chem import AllChem as Chem
import sys
import os
molfile = sys.argv[1]
#filepath = sys.argv[2]
mol = Chem.MolFromMolFile(molfile+".mol")
mol.Compute2DCoords()
img = Draw.MolsToGridImage([mol],molsPerRow=1,useSVG=True)
with open(molfile+".svg", "w") as f:  
    f.write(img)
