import unittest
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen


class LigandLipinskiTest(unittest.TestCase):
    """Sanity tests for Lipinski's Rule of Five compliance"""

    def setUp(self):
        """Set up test ligands with known Lipinski compliance"""
        # Known drug-like molecules (passing Lipinski)
        self.drug_like_smiles = {
            'aspirin': 'CC(=O)Oc1ccccc1C(=O)O',  # MW: 180, logP: 1.19, HBD: 1, HBA: 4
            'ibuprofen': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',  # MW: 206, logP: 3.97, HBD: 1, HBA: 2
            'caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # MW: 194, logP: 0.16, HBD: 0, HBA: 3
            'paracetamol': 'CC(=O)Nc1ccc(O)cc1',  # MW: 151, logP: 0.46, HBD: 2, HBA: 2
        }
        
        # Known non-drug-like molecules (violating Lipinski)
        self.non_drug_like_smiles = {
            'large_molecule': 'C' * 150,  # Very large MW > 500
            'polar_molecule': 'OCCO' * 20,  # Many H-bond donors/acceptors
        }

    def test_lipinski_molecular_weight(self):
        """Test that drug-like molecules have MW <= 500"""
        for name, smiles in self.drug_like_smiles.items():
            mol = Chem.MolFromSmiles(smiles)
            self.assertIsNotNone(mol, f"Failed to parse SMILES for {name}")
            mw = Descriptors.MolWt(mol)
            self.assertLessEqual(mw, 500, f"{name}: MW {mw} exceeds 500")

    def test_lipinski_logp(self):
        """Test that drug-like molecules have logP <= 5"""
        for name, smiles in self.drug_like_smiles.items():
            mol = Chem.MolFromSmiles(smiles)
            self.assertIsNotNone(mol, f"Failed to parse SMILES for {name}")
            logp = Crippen.MolLogP(mol)
            self.assertLessEqual(logp, 5, f"{name}: logP {logp} exceeds 5")

    def test_lipinski_hbd(self):
        """Test that drug-like molecules have HBD <= 5"""
        for name, smiles in self.drug_like_smiles.items():
            mol = Chem.MolFromSmiles(smiles)
            self.assertIsNotNone(mol, f"Failed to parse SMILES for {name}")
            hbd = Descriptors.NumHDonors(mol)
            self.assertLessEqual(hbd, 5, f"{name}: HBD {hbd} exceeds 5")

    def test_lipinski_hba(self):
        """Test that drug-like molecules have HBA <= 10"""
        for name, smiles in self.drug_like_smiles.items():
            mol = Chem.MolFromSmiles(smiles)
            self.assertIsNotNone(mol, f"Failed to parse SMILES for {name}")
            hba = Descriptors.NumHAcceptors(mol)
            self.assertLessEqual(hba, 10, f"{name}: HBA {hba} exceeds 10")

    def test_non_drug_like_molecules(self):
        """Test that non-drug-like molecules violate at least one Lipinski rule"""
        for name, smiles in self.non_drug_like_smiles.items():
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            
            # At least one rule should be violated
            violation = mw > 500 or logp > 5 or hbd > 5 or hba > 10
            self.assertTrue(violation, f"{name}: Expected to violate Lipinski's rule")

    def test_lipinski_complete_rule(self):
        """Test complete Lipinski's Rule of Five for drug-like molecules"""
        for name, smiles in self.drug_like_smiles.items():
            mol = Chem.MolFromSmiles(smiles)
            self.assertIsNotNone(mol, f"Failed to parse SMILES for {name}")
            
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            
            # All rules should pass
            self.assertLessEqual(mw, 500, f"{name} MW violation")
            self.assertLessEqual(logp, 5, f"{name} logP violation")
            self.assertLessEqual(hbd, 5, f"{name} HBD violation")
            self.assertLessEqual(hba, 10, f"{name} HBA violation")


if __name__ == '__main__':
    unittest.main()
