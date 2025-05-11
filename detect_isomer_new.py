from rdkit import Chem

def detect_potential_stereochemistry(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Invalid SMILES input.")
        return

    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    potential_ez = []

    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and not bond.IsInRing():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()

            begin_neighbors = []
            end_neighbors = []

            # Collect begin_atom neighbors + implicit H
            for nbr in begin_atom.GetNeighbors():
                if nbr.GetIdx() != end_atom.GetIdx():
                    begin_neighbors.append(nbr.GetSymbol())
            begin_neighbors += ["H"] * begin_atom.GetTotalNumHs(includeNeighbors=False)

            # Collect end_atom neighbors + implicit H
            for nbr in end_atom.GetNeighbors():
                if nbr.GetIdx() != begin_atom.GetIdx():
                    end_neighbors.append(nbr.GetSymbol())
            end_neighbors += ["H"] * end_atom.GetTotalNumHs(includeNeighbors=False)

            # Decision logic
            if len(begin_neighbors) >= 1 and len(end_neighbors) >= 1:
                if len(set(begin_neighbors)) > 1 or len(set(end_neighbors)) > 1:
                    potential_ez.append((bond.GetIdx(), begin_atom.GetIdx(), end_atom.GetIdx()))

    # Output
    if chiral_centers:
        print(f"Potential Chiral Centers (atom_idx, R/S if known): {chiral_centers}")
    else:
        print("No potential chiral centers detected.")

    if potential_ez:
        print(f"Potential E/Z double bonds (bond_idx, atom1_idx, atom2_idx): {potential_ez}")
    else:
        print("No potential E/Z double bonds detected.")

# Example usage
if __name__ == "__main__":
    smiles_list = [
        "C1CCC=CC1",     
        "CC(C)=C(C)C",    
        "C/C=C/C",       
        "CC=CC",         
        "CC(C)(F)Cl",    
        "CC(C)Br",       
        "CC(C(=O)O)O",   
        "CC(C(C)Br)Cl"  
    ]
    for smiles in smiles_list:
        print(f"\nAnalyzing SMILES: {smiles}")
        detect_potential_stereochemistry(smiles)
