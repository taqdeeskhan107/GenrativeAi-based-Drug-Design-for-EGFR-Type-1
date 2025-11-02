#!/usr/bin/env python3
"""
Co-crystal Ligand Extractor
Extracts ligand structures from protein-ligand co-crystal structures.
Author: TAQDEES
"""

import os
import sys
import json
import numpy as np
from pathlib import Path
import argparse
import logging
from datetime import datetime

# BioPython imports for structure handling
try:
    from Bio.PDB import PDBParser, PDBIO, Select
    print("✅ BioPython imported successfully")
except ImportError:
    print("❌ BioPython not found. Please install: pip install biopython")
    exit(1)

# Optional RDKit import for additional ligand processing
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
    print("✅ RDKit available for ligand processing")
except ImportError:
    print("⚠️ RDKit not available. Limited ligand processing capabilities.")
    RDKIT_AVAILABLE = False

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class LigandSelector(Select):
    """Custom selector for extracting ligand residues."""
    
    def __init__(self, ligand_id):
        self.ligand_id = ligand_id
    
    def accept_residue(self, residue):
        """Accept only the specified ligand residue."""
        return residue.get_resname() == self.ligand_id

def extract_ligand(pdb_file, ligand_id, output_dir=None, analyze=True):
    """
    Extract ligand structure from a protein-ligand complex.
    
    Args:
        pdb_file (str): Path to input PDB file
        ligand_id (str): Three-letter code of the ligand (e.g., 'AQ4' for erlotinib)
        output_dir (str): Directory for output files (default: same as input)
        analyze (bool): Whether to analyze ligand properties
    
    Returns:
        dict: Ligand information and analysis results
    """
    
    # Setup output directory
    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(pdb_file))
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Initialize results dictionary
    results = {
        'input_file': pdb_file,
        'ligand_id': ligand_id,
        'extraction_time': datetime.now().isoformat(),
        'success': False,
        'error': None,
        'ligand_info': {}
    }
    
    try:
        # Parse PDB structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('complex', pdb_file)
        
        # Find ligand
        ligand_atoms = []
        ligand_residue = None
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() == ligand_id:
                        ligand_residue = residue
                        ligand_atoms = [atom for atom in residue.get_atoms()]
                        logger.info(f"✅ Found ligand {ligand_id} in chain {chain.id}")
                        break
                if ligand_atoms:
                    break
            if ligand_atoms:
                break
        
        if not ligand_atoms:
            error_msg = f"Ligand {ligand_id} not found in structure!"
            logger.error(f"❌ {error_msg}")
            results['error'] = error_msg
            return results
        
        # Basic ligand information
        results['ligand_info'] = {
            'n_atoms': len(ligand_atoms),
            'chain': ligand_residue.get_parent().id,
            'residue_number': ligand_residue.get_id()[1],
            'center_of_mass': {
                'x': float(np.mean([atom.get_coord()[0] for atom in ligand_atoms])),
                'y': float(np.mean([atom.get_coord()[1] for atom in ligand_atoms])),
                'z': float(np.mean([atom.get_coord()[2] for atom in ligand_atoms]))
            }
        }
        
        # Extract ligand to separate PDB file
        io = PDBIO()
        io.set_structure(structure)
        
        ligand_file = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(pdb_file))[0]}_{ligand_id}.pdb")
        io.save(ligand_file, LigandSelector(ligand_id))
        logger.info(f"💾 Saved ligand structure to: {ligand_file}")
        
        results['output_file'] = ligand_file
        results['success'] = True
        
        # Additional analysis with RDKit if available
        if analyze and RDKIT_AVAILABLE:
            try:
                # Convert PDB to mol using RDKit
                mol = Chem.MolFromPDBFile(ligand_file)
                if mol is not None:
                    # Calculate basic properties
                    results['ligand_info']['properties'] = {
                        'molecular_weight': float(Chem.Descriptors.ExactMolWt(mol)),
                        'rotatable_bonds': int(Chem.Descriptors.NumRotatableBonds(mol)),
                        'h_bond_donors': int(Chem.Descriptors.NumHDonors(mol)),
                        'h_bond_acceptors': int(Chem.Descriptors.NumHAcceptors(mol)),
                        'logp': float(Chem.Descriptors.MolLogP(mol)),
                        'tpsa': float(Chem.Descriptors.TPSA(mol)),
                        'aromatic_rings': int(Chem.Descriptors.NumAromaticRings(mol))
                    }
                    
                    # Generate 2D depiction
                    img_file = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(pdb_file))[0]}_{ligand_id}.png")
                    AllChem.Compute2DCoords(mol)
                    Chem.Draw.MolToFile(mol, img_file)
                    results['ligand_info']['2d_image'] = img_file
                    logger.info(f"🎨 Generated 2D structure image: {img_file}")
            
            except Exception as e:
                logger.warning(f"⚠️ RDKit analysis failed: {str(e)}")
                results['ligand_info']['rdkit_analysis_error'] = str(e)
        
        # Save results summary
        summary_file = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(pdb_file))[0]}_{ligand_id}_info.json")
        with open(summary_file, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"📊 Saved ligand information to: {summary_file}")
        
    except Exception as e:
        error_msg = str(e)
        logger.error(f"❌ Error during ligand extraction: {error_msg}")
        results['error'] = error_msg
    
    return results

def main():
    """Main function for command-line usage."""
    
    parser = argparse.ArgumentParser(
        description="Extract and analyze co-crystallized ligands from PDB structures."
    )
    parser.add_argument("pdb_file", help="Input PDB file path")
    parser.add_argument("ligand_id", help="Three-letter code of the ligand (e.g., AQ4)")
    parser.add_argument("--output-dir", help="Output directory (default: same as input)")
    parser.add_argument("--no-analysis", action="store_true", help="Skip ligand analysis")
    
    args = parser.parse_args()
    
    # Print banner
    print("\n🧬 Co-crystal Ligand Extractor")
    print("="*40)
    
    # Extract ligand
    results = extract_ligand(
        args.pdb_file,
        args.ligand_id,
        args.output_dir,
        not args.no_analysis
    )
    
    # Print results summary
    if results['success']:
        print("\n✅ Ligand extraction successful!")
        print(f"📊 Extracted {results['ligand_info']['n_atoms']} atoms")
        print(f"💾 Output file: {results['output_file']}")
        
        if 'properties' in results['ligand_info']:
            props = results['ligand_info']['properties']
            print("\n📝 Ligand Properties:")
            print(f"   Molecular Weight: {props['molecular_weight']:.1f}")
            print(f"   Rotatable Bonds: {props['rotatable_bonds']}")
            print(f"   H-bond Donors: {props['h_bond_donors']}")
            print(f"   H-bond Acceptors: {props['h_bond_acceptors']}")
            print(f"   LogP: {props['logp']:.2f}")
            print(f"   TPSA: {props['tpsa']:.1f}")
            print(f"   Aromatic Rings: {props['aromatic_rings']}")
    else:
        print(f"\n❌ Ligand extraction failed: {results['error']}")
    
    print("\n" + "="*40)

if __name__ == "__main__":
    main()
