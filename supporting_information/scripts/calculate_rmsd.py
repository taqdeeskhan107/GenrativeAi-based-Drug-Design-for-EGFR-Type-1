#!/usr/bin/env python3
"""
RMSD Calculator for EGFR Type 1 Inhibitors
Calculates Root Mean Square Deviation between predicted and crystal poses.
Author: TAQDEES
"""

import os
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')
import json
import pandas as pd

def calculate_atomic_rmsd(coords1, coords2):
    """
    Calculate RMSD between two sets of coordinates.
    
    Args:
        coords1 (np.array): First set of 3D coordinates
        coords2 (np.array): Second set of 3D coordinates
        
    Returns:
        float: RMSD value
    """
    if len(coords1) != len(coords2):
        raise ValueError("Coordinate sets must have same number of atoms")
    
    squared_diff = np.sum((coords1 - coords2) ** 2, axis=1)
    rmsd = np.sqrt(np.mean(squared_diff))
    return rmsd

def load_coordinates(pdb_file, ligand_only=True):
    """
    Load atomic coordinates from PDB file.
    
    Args:
        pdb_file (str): Path to PDB file
        ligand_only (bool): If True, only load ligand coordinates
        
    Returns:
        np.array: Array of coordinates
    """
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or (not ligand_only and line.startswith('HETATM')):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    return np.array(coords)

def calculate_pose_rmsd(crystal_file, docked_file, output_file=None):
    """
    Calculate RMSD between crystal and docked poses.
    
    Args:
        crystal_file (str): Path to crystal structure PDB
        docked_file (str): Path to docked pose PDB
        output_file (str): Optional path to save results
        
    Returns:
        dict: RMSD calculation results
    """
    try:
        # Load coordinates
        crystal_coords = load_coordinates(crystal_file)
        docked_coords = load_coordinates(docked_file)
        
        # Calculate RMSD
        rmsd = calculate_atomic_rmsd(crystal_coords, docked_coords)
        
        results = {
            'rmsd': float(rmsd),
            'crystal_atoms': len(crystal_coords),
            'docked_atoms': len(docked_coords),
            'success': True
        }
        
        # Save results if output file specified
        if output_file:
            with open(output_file, 'w') as f:
                json.dump(results, f, indent=4)
        
        return results
    
    except Exception as e:
        print(f"Error calculating RMSD: {e}")
        return {'success': False, 'error': str(e)}

def calculate_multiple_rmsd(crystal_dir, docked_dir, output_file=None):
    """
    Calculate RMSD for multiple poses.
    
    Args:
        crystal_dir (str): Directory with crystal structures
        docked_dir (str): Directory with docked poses
        output_file (str): Optional path to save results
        
    Returns:
        pd.DataFrame: RMSD results for all poses
    """
    results = []
    
    crystal_files = Path(crystal_dir).glob('*.pdb')
    for crystal_file in crystal_files:
        ligand_name = crystal_file.stem
        docked_file = Path(docked_dir) / f"{ligand_name}_docked.pdb"
        
        if docked_file.exists():
            rmsd_result = calculate_pose_rmsd(str(crystal_file), str(docked_file))
            if rmsd_result['success']:
                results.append({
                    'ligand': ligand_name,
                    'rmsd': rmsd_result['rmsd'],
                    'crystal_atoms': rmsd_result['crystal_atoms'],
                    'docked_atoms': rmsd_result['docked_atoms']
                })
    
    df = pd.DataFrame(results)
    
    if output_file:
        df.to_csv(output_file, index=False)
    
    return df

def main():
    """Main function for RMSD calculation"""
    print("🎯 EGFR Type 1 Inhibitor RMSD Calculator")
    print("="*50)
    
    # Default directories
    crystal_dir = "crystal_structures"
    docked_dir = "docked_poses"
    output_file = "rmsd_results.csv"
    
    if not os.path.exists(crystal_dir):
        print(f"❌ Crystal structure directory not found: {crystal_dir}")
        return
    
    if not os.path.exists(docked_dir):
        print(f"❌ Docked poses directory not found: {docked_dir}")
        return
    
    print(f"📊 Calculating RMSD values...")
    results_df = calculate_multiple_rmsd(crystal_dir, docked_dir, output_file)
    
    # Print summary
    print("\n📋 RMSD Calculation Summary")
    print("-"*30)
    print(f"Total poses analyzed: {len(results_df)}")
    print(f"Average RMSD: {results_df['rmsd'].mean():.2f} Å")
    print(f"Min RMSD: {results_df['rmsd'].min():.2f} Å")
    print(f"Max RMSD: {results_df['rmsd'].max():.2f} Å")
    
    # Success criteria
    good_poses = (results_df['rmsd'] < 2.0).sum()
    print(f"\n✅ Poses with RMSD < 2.0 Å: {good_poses}/{len(results_df)}")
    print(f"Success rate: {(good_poses/len(results_df))*100:.1f}%")
    
    print(f"\n📁 Results saved to: {output_file}")

if __name__ == "__main__":
    main()
