"""
This is used to extract the carbon skeleton of terpenoids
"""
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import re

def get_skeleton(smi, cat=None):
    mol = Chem.RWMol(Chem.MolFromSmiles(smi))
    # delete aromatic flag on atoms
    Chem.Kekulize(mol, clearAromaticFlags=True)

    # search c-o, c-n bond
    bond = '[#6]~[!#6]'
    bond_mol = Chem.MolFromSmarts(bond)
    bs = mol.GetSubstructMatches(bond_mol)

    # match and break the bond
    for match in bs:
        mol.RemoveBond(match[0], match[1])
    
    frags = Chem.GetMolFrags(mol, asMols=True)

    # select scaffold
    core_smi = ""
    for frag in frags:
        scaffold = MurckoScaffold.MakeScaffoldGeneric(frag)
        if not cat:
            if abs(Chem.MolToSmiles(scaffold, isomericSmiles=False).count('C')) > abs(core_smi.count('C')):
                core_smi = Chem.MolToSmiles(scaffold)
        if cat == 'Sesquiterpenoids':
            if abs(Chem.MolToSmiles(scaffold, isomericSmiles=False).count('C')-15) < abs(core_smi.count('C')-15):
                core_smi = Chem.MolToSmiles(scaffold)
        if cat == 'Diterpenoids':
            if abs(Chem.MolToSmiles(scaffold, isomericSmiles=False).count('C')-20) < abs(core_smi.count('C')-20):
                core_smi = Chem.MolToSmiles(scaffold)

        if cat == 'Triterpenoids':
            if abs(Chem.MolToSmiles(scaffold, isomericSmiles=False).count('C')-30) < abs(core_smi.count('C')-30):
                core_smi = Chem.MolToSmiles(scaffold)
    core_smi_new = re.sub('\[C.{0,2}]', 'C', core_smi)

    return core_smi_new

if __name__ == "__main__":
    smi = "CC(C)=CCC[C@@H](C)C1=CCC(C)=CC1"
    #category = "Sesquiterpenoids"
    print(get_skeleton(smi))
