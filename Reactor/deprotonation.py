"""
This is used to exhaustively deprotonate the carbacation
"""
from rdkit import Chem

if __name__=="__main__":
    cations = []
    for line in open("network.txt"):
        smi = line.split("\t")[3]
        if line.split("\t")[3] in cations:
            continue
        cations.append(smi)
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        flag = False
        for atom in mol.GetAtoms():
            if atom.GetFormalCharge() == 1:
                cation_idx = atom.GetIdx()
                flag = True
        if not flag:
            continue
        for atom in mol.GetAtomWithIdx(cation_idx).GetNeighbors():
            if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                for atom1 in atom.GetNeighbors():
                    if atom1.GetSymbol() == "H":
                        mol_new = Chem.RWMol(mol)
                        mol_new.RemoveAtom(atom1.GetIdx())
                        mol_new.GetAtomWithIdx(cation_idx).SetFormalCharge(0)
                        mol_new.GetBondBetweenAtoms(cation_idx, atom.GetIdx()).SetBondType(Chem.rdchem.BondType.DOUBLE)
                        mol_new = Chem.RemoveHs(mol_new)
                        with open("deprotonate.txt", "a") as f:
                            f.write(smi + "\t" + Chem.MolToSmiles(mol_new) + "\n")
                        break