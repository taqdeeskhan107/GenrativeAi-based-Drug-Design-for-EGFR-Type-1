#!/usr/bin/env python3
"""
Receptor Preparation for Cross-Docking
Prepares protein receptor structures for cross-docking validation.
Author: TAQDEES
"""

import os
import sys
import json
import logging
from pathlib import Path
import subprocess
import argparse
from datetime import datetime

# BioPython imports for structure processing
try:
    from Bio.PDB import PDBParser, PDBIO, Select, Structure
    from Bio.PDB.Polypeptide import three_to_one
    print("✅ BioPython imported successfully")
except ImportError:
    print("❌ BioPython not found. Please install: pip install biopython")
    exit(1)

# Optional import for structure optimization
try:
    from pymol import cmd
    PYMOL_AVAILABLE = True
    print("✅ PyMOL available for structure optimization")
except ImportError:
    print("⚠️ PyMOL not available. Limited structure optimization capabilities.")
    PYMOL_AVAILABLE = False

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class ProteinSelect(Select):
    """Custom selector for protein residues."""
    
    def accept_residue(self, residue):
        """Accept only standard amino acid residues."""
        return residue.get_id()[0] == ' '  # Standard amino acids have blank insertion code

def download_pdb(pdb_id, output_dir):
    """Download PDB structure if not already present."""
    
    pdb_file = os.path.join(output_dir, f"{pdb_id}.pdb")
    if os.path.exists(pdb_file):
        logger.info(f"📂 Using existing PDB file: {pdb_file}")
        return pdb_file
    
    logger.info(f"📥 Downloading {pdb_id}...")
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        subprocess.run(['wget', '-O', pdb_file, url], check=True)
        logger.info(f"✅ Downloaded {pdb_id} successfully")
        return pdb_file
    except subprocess.CalledProcessError as e:
        logger.error(f"❌ Failed to download {pdb_id}: {e}")
        return None

def clean_structure(pdb_file, chain_id=None):
    """Clean protein structure by removing waters, ligands, and other non-protein components."""
    
    logger.info(f"🧹 Cleaning structure: {pdb_file}")
    parser = PDBParser(QUIET=True)
    
    try:
        structure = parser.get_structure('protein', pdb_file)
        clean_struct = Structure.Structure('clean')
        model = structure[0]  # First model
        
        clean_model = model.copy()
        clean_model.detach_parent()
        
        # Remove non-protein residues and specified chains
        for chain in list(clean_model.get_chains()):
            if chain_id and chain.id != chain_id:
                clean_model.detach_child(chain.id)
                continue
                
            for residue in list(chain.get_residues()):
                if residue.get_id()[0] != ' ':  # Non-protein residue
                    chain.detach_child(residue.get_id())
        
        clean_struct.add(clean_model)
        
        # Save clean structure
        output_file = os.path.splitext(pdb_file)[0] + "_clean.pdb"
        io = PDBIO()
        io.set_structure(clean_struct)
        io.save(output_file)
        
        logger.info(f"✅ Saved clean structure to: {output_file}")
        return output_file
        
    except Exception as e:
        logger.error(f"❌ Error cleaning structure: {e}")
        return None

def optimize_structure(pdb_file):
    """Optimize protein structure using PyMOL (if available)."""
    
    if not PYMOL_AVAILABLE:
        logger.warning("⚠️ PyMOL not available, skipping structure optimization")
        return pdb_file
    
    logger.info(f"🔧 Optimizing structure: {pdb_file}")
    try:
        output_file = os.path.splitext(pdb_file)[0] + "_optimized.pdb"
        
        cmd.load(pdb_file, "protein")
        
        # Basic structure optimization
        cmd.h_add("protein")  # Add hydrogens
        cmd.remove("solvent")  # Remove water molecules
        cmd.sort()  # Sort atoms
        
        # Energy minimization (if supported)
        try:
            cmd.minimize("protein", method="conjugate_gradients", nsteps=500)
        except:
            logger.warning("⚠️ Energy minimization not supported in this PyMOL version")
        
        cmd.save(output_file, "protein")
        cmd.delete("protein")
        
        logger.info(f"✅ Saved optimized structure to: {output_file}")
        return output_file
        
    except Exception as e:
        logger.error(f"❌ Error optimizing structure: {e}")
        return pdb_file

def prepare_receptor(pdb_file, working_dir, box_center=None, box_size=None):
    """
    Prepare receptor for docking using AutoDockTools.
    
    Args:
        pdb_file (str): Input PDB file path
        working_dir (str): Working directory for output files
        box_center (tuple): Optional (x, y, z) coordinates for grid box center
        box_size (tuple): Optional (x, y, z) dimensions for grid box
    """
    
    logger.info(f"🔄 Preparing receptor: {pdb_file}")
    
    try:
        # Convert to PDBQT format
        receptor_name = os.path.splitext(os.path.basename(pdb_file))[0]
        output_file = os.path.join(working_dir, f"{receptor_name}.pdbqt")
        
        subprocess.run([
            'prepare_receptor4.py',
            '-r', pdb_file,
            '-o', output_file,
            '-A', 'hydrogens',  # Add hydrogens
            '-U', 'nphs_lps_waters'  # Remove non-polar hydrogens and waters
        ], check=True)
        
        logger.info(f"✅ Generated PDBQT file: {output_file}")
        
        # Create configuration file if box parameters provided
        if box_center and box_size:
            config_file = os.path.join(working_dir, f"{receptor_name}_config.txt")
            with open(config_file, 'w') as f:
                f.write(f"center_x = {box_center[0]}\n")
                f.write(f"center_y = {box_center[1]}\n")
                f.write(f"center_z = {box_center[2]}\n")
                f.write(f"size_x = {box_size[0]}\n")
                f.write(f"size_y = {box_size[1]}\n")
                f.write(f"size_z = {box_size[2]}\n")
                f.write("exhaustiveness = 16\n")
            
            logger.info(f"📝 Created configuration file: {config_file}")
            
        return output_file
        
    except subprocess.CalledProcessError as e:
        logger.error(f"❌ Error preparing receptor: {e}")
        return None

def prepare_egfr_receptors(pdb_ids, output_dir, optimize=True):
    """
    Prepare multiple EGFR structures for cross-docking.
    
    Args:
        pdb_ids (list): List of PDB IDs to prepare
        output_dir (str): Output directory
        optimize (bool): Whether to optimize structures with PyMOL
    """
    
    # Create output directories
    working_dir = os.path.join(output_dir, "receptors")
    Path(working_dir).mkdir(parents=True, exist_ok=True)
    
    results = []
    
    # Default grid box parameters for EGFR ATP binding site
    default_box = {
        'center': (70.0, 50.0, 55.0),
        'size': (25.0, 25.0, 25.0)
    }
    
    for pdb_id in pdb_ids:
        logger.info(f"\n🔄 Processing {pdb_id}...")
        
        result = {
            'pdb_id': pdb_id,
            'success': False,
            'files': {},
            'timestamp': datetime.now().isoformat()
        }
        
        try:
            # Step 1: Download structure
            pdb_file = download_pdb(pdb_id, working_dir)
            if not pdb_file:
                raise Exception(f"Failed to download {pdb_id}")
            result['files']['original'] = pdb_file
            
            # Step 2: Clean structure
            clean_file = clean_structure(pdb_file)
            if not clean_file:
                raise Exception(f"Failed to clean {pdb_id}")
            result['files']['clean'] = clean_file
            
            # Step 3: Optimize structure (if enabled)
            if optimize:
                opt_file = optimize_structure(clean_file)
                if opt_file:
                    result['files']['optimized'] = opt_file
                    prep_input = opt_file
            else:
                prep_input = clean_file
            
            # Step 4: Prepare for docking
            pdbqt_file = prepare_receptor(
                prep_input,
                working_dir,
                box_center=default_box['center'],
                box_size=default_box['size']
            )
            
            if not pdbqt_file:
                raise Exception(f"Failed to prepare {pdb_id} for docking")
            result['files']['pdbqt'] = pdbqt_file
            
            result['success'] = True
            logger.info(f"✅ Successfully prepared {pdb_id}")
            
        except Exception as e:
            logger.error(f"❌ Error processing {pdb_id}: {e}")
            result['error'] = str(e)
        
        results.append(result)
    
    # Save summary
    summary_file = os.path.join(output_dir, "receptor_preparation_summary.json")
    with open(summary_file, 'w') as f:
        json.dump({
            'receptors': results,
            'timestamp': datetime.now().isoformat(),
            'total_success': sum(1 for r in results if r['success']),
            'total_failed': sum(1 for r in results if not r['success'])
        }, f, indent=2)
    
    logger.info(f"\n📊 Summary saved to: {summary_file}")
    return results

def main():
    """Main function for command-line usage."""
    
    parser = argparse.ArgumentParser(
        description="Prepare protein structures for cross-docking validation."
    )
    parser.add_argument(
        "--pdb-ids",
        nargs="+",
        default=["1M17", "2ITY", "2J6M", "4WKQ", "4ZAU"],
        help="PDB IDs to prepare (default: EGFR Type 1 inhibitor structures)"
    )
    parser.add_argument(
        "--output-dir",
        default="cross_docking_workspace",
        help="Output directory for prepared structures"
    )
    parser.add_argument(
        "--no-optimize",
        action="store_true",
        help="Skip structure optimization step"
    )
    
    args = parser.parse_args()
    
    # Print banner
    print("\n🧬 Cross-Docking Receptor Preparation")
    print("="*40)
    
    # Prepare receptors
    results = prepare_egfr_receptors(
        args.pdb_ids,
        args.output_dir,
        not args.no_optimize
    )
    
    # Print final summary
    print("\n📊 Preparation Summary")
    print("-"*20)
    successful = sum(1 for r in results if r['success'])
    print(f"✅ Successfully prepared: {successful}/{len(results)}")
    if successful < len(results):
        print("\n❌ Failed preparations:")
        for r in results:
            if not r['success']:
                print(f"   {r['pdb_id']}: {r.get('error', 'Unknown error')}")
    
    print("\n" + "="*40)

if __name__ == "__main__":
    main()
