#!/usr/bin/env python3
"""
Cross-Docking Validation Pipeline
Author: TAQDEES
Description: Performs cross-docking validation using multiple EGFR crystal structures
"""

import os
import pandas as pd
import numpy as np
import subprocess
import logging
from datetime import datetime
from pathlib import Path
import json
import shutil

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
    print("✅ RDKit available for molecule processing")
except ImportError:
    print("❌ RDKit not available. Please install: conda install -c conda-forge rdkit")
    RDKIT_AVAILABLE = False
    exit(1)

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class CrossDockingValidator:
    """Cross-docking validation for EGFR Type 1 inhibitors"""
    
    def __init__(self, output_dir="cross_docking_validation"):
        """Initialize cross-docking validator"""
        self.output_dir = output_dir
        self.structures_dir = os.path.join(output_dir, "structures")
        self.results_dir = os.path.join(output_dir, "results")
        self.ligands_dir = os.path.join(output_dir, "ligands")
        
        # Create directories
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.structures_dir, exist_ok=True)
        os.makedirs(self.results_dir, exist_ok=True)
        os.makedirs(self.ligands_dir, exist_ok=True)
        
        # PDB IDs for cross-docking (EGFR structures with Type 1 inhibitors)
        self.pdb_structures = [
            "1M17",  # Erlotinib
            "2ITY",  # Gefitinib
            "2J6M",  # Lapatinib
            "4WKQ",  # Afatinib
            "4ZAU"   # Icotinib
        ]
        
        logger.info(f"🔄 Initialized cross-docking validation in: {output_dir}")

    def prepare_structures(self):
        """Prepare receptor structures for cross-docking"""
        
        logger.info("🧬 Preparing receptor structures...")
        results = []
        
        for pdb_id in self.pdb_structures:
            try:
                # Download PDB if not exists
                pdb_file = os.path.join(self.structures_dir, f"{pdb_id}.pdb")
                if not os.path.exists(pdb_file):
                    logger.info(f"📥 Downloading {pdb_id}...")
                    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                    subprocess.run(['wget', '-O', pdb_file, url], check=True)
                
                # Clean structure
                clean_file = os.path.join(self.structures_dir, f"{pdb_id}_clean.pdb")
                logger.info(f"🧹 Cleaning {pdb_id}...")
                
                with open(pdb_file, 'r') as f_in, open(clean_file, 'w') as f_out:
                    for line in f_in:
                        if line.startswith('ATOM') and 'A' in line[21]:
                            f_out.write(line)
                
                # Convert to PDBQT
                receptor_file = os.path.join(self.structures_dir, f"{pdb_id}_receptor.pdbqt")
                logger.info(f"🛠️ Converting {pdb_id} to PDBQT...")
                subprocess.run([
                    'prepare_receptor4.py',
                    '-r', clean_file,
                    '-o', receptor_file
                ], check=True)
                
                results.append({
                    'pdb_id': pdb_id,
                    'success': True,
                    'files': {
                        'original': pdb_file,
                        'clean': clean_file,
                        'pdbqt': receptor_file
                    }
                })
                
                logger.info(f"✅ {pdb_id} prepared successfully")
                
            except Exception as e:
                logger.error(f"❌ Error preparing {pdb_id}: {e}")
                results.append({
                    'pdb_id': pdb_id,
                    'success': False,
                    'error': str(e)
                })
        
        return results

    def extract_native_ligands(self):
        """Extract co-crystallized ligands from structures"""
        
        logger.info("🔍 Extracting native ligands...")
        results = []
        
        for pdb_id in self.pdb_structures:
            try:
                pdb_file = os.path.join(self.structures_dir, f"{pdb_id}.pdb")
                ligand_file = os.path.join(self.ligands_dir, f"{pdb_id}_native.pdb")
                
                # Extract HETATM lines for ligand
                with open(pdb_file, 'r') as f_in, open(ligand_file, 'w') as f_out:
                    for line in f_in:
                        if line.startswith('HETATM') and 'A' in line[21]:
                            f_out.write(line)
                
                # Convert to PDBQT
                ligand_pdbqt = os.path.join(self.ligands_dir, f"{pdb_id}_native.pdbqt")
                subprocess.run([
                    'prepare_ligand4.py',
                    '-l', ligand_file,
                    '-o', ligand_pdbqt
                ], check=True)
                
                results.append({
                    'pdb_id': pdb_id,
                    'success': True,
                    'files': {
                        'pdb': ligand_file,
                        'pdbqt': ligand_pdbqt
                    }
                })
                
                logger.info(f"✅ Extracted ligand from {pdb_id}")
                
            except Exception as e:
                logger.error(f"❌ Error extracting ligand from {pdb_id}: {e}")
                results.append({
                    'pdb_id': pdb_id,
                    'success': False,
                    'error': str(e)
                })
        
        return results

    def run_cross_docking(self):
        """Perform cross-docking with all ligands against all receptors"""
        
        logger.info("🔄 Starting cross-docking validation...")
        results = []
        
        for receptor_id in self.pdb_structures:
            receptor_file = os.path.join(self.structures_dir, f"{receptor_id}_receptor.pdbqt")
            
            for ligand_id in self.pdb_structures:
                ligand_file = os.path.join(self.ligands_dir, f"{ligand_id}_native.pdbqt")
                output_file = os.path.join(self.results_dir, f"{ligand_id}_vs_{receptor_id}.pdbqt")
                
                try:
                    # Run AutoDock Vina
                    logger.info(f"🎯 Docking {ligand_id} ligand against {receptor_id} receptor...")
                    
                    # Use standard grid parameters for EGFR ATP binding site
                    subprocess.run([
                        'vina',
                        '--receptor', receptor_file,
                        '--ligand', ligand_file,
                        '--center_x', '70.0',
                        '--center_y', '50.0',
                        '--center_z', '55.0',
                        '--size_x', '25.0',
                        '--size_y', '25.0',
                        '--size_z', '25.0',
                        '--exhaustiveness', '16',
                        '--out', output_file
                    ], check=True)
                    
                    # Parse docking score
                    score = float('inf')
                    with open(output_file, 'r') as f:
                        for line in f:
                            if "REMARK VINA RESULT" in line:
                                score = float(line.split()[3])
                                break
                    
                    results.append({
                        'ligand_id': ligand_id,
                        'receptor_id': receptor_id,
                        'docking_score': score,
                        'success': True,
                        'output_file': output_file
                    })
                    
                    logger.info(f"✅ {ligand_id} vs {receptor_id}: {score:.2f} kcal/mol")
                    
                except Exception as e:
                    logger.error(f"❌ Error in {ligand_id} vs {receptor_id}: {e}")
                    results.append({
                        'ligand_id': ligand_id,
                        'receptor_id': receptor_id,
                        'success': False,
                        'error': str(e)
                    })
        
        return pd.DataFrame(results)

    def analyze_results(self, df):
        """Analyze cross-docking results"""
        
        logger.info("📊 Analyzing cross-docking results...")
        
        # Filter successful docking runs
        successful = df[df['success'] == True].copy()
        
        if len(successful) == 0:
            logger.error("No successful docking results to analyze")
            return None
        
        # Create cross-docking matrix
        matrix = successful.pivot(
            index='ligand_id',
            columns='receptor_id',
            values='docking_score'
        )
        
        # Calculate statistics
        stats = {
            'total_runs': len(df),
            'successful_runs': len(successful),
            'success_rate': len(successful) / len(df) * 100,
            'average_score': successful['docking_score'].mean(),
            'best_score': successful['docking_score'].min(),
            'cross_docking_matrix': matrix.to_dict(),
            'self_docking_scores': []
        }
        
        # Analyze self-docking results
        for pdb_id in self.pdb_structures:
            self_dock = successful[
                (successful['ligand_id'] == pdb_id) &
                (successful['receptor_id'] == pdb_id)
            ]
            if len(self_dock) > 0:
                stats['self_docking_scores'].append({
                    'pdb_id': pdb_id,
                    'score': float(self_dock['docking_score'].iloc[0])
                })
        
        return stats

    def generate_report(self, stats):
        """Generate validation report"""
        
        report_file = os.path.join(self.output_dir, 'cross_docking_report.md')
        
        report = f"""# Cross-Docking Validation Report
Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## Summary
- Total docking runs: {stats['total_runs']}
- Successful runs: {stats['successful_runs']} ({stats['success_rate']:.1f}%)
- Average docking score: {stats['average_score']:.2f} kcal/mol
- Best docking score: {stats['best_score']:.2f} kcal/mol

## Self-Docking Results
"""
        
        for result in stats['self_docking_scores']:
            report += f"- {result['pdb_id']}: {result['score']:.2f} kcal/mol\n"
        
        with open(report_file, 'w') as f:
            f.write(report)
        
        # Save detailed statistics
        stats_file = os.path.join(self.output_dir, 'cross_docking_statistics.json')
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=4)
        
        logger.info(f"📋 Report saved to: {report_file}")
        logger.info(f"📊 Statistics saved to: {stats_file}")

def main():
    """Main function for cross-docking validation"""
    
    print("🔄 EGFR Type 1 Inhibitor Cross-Docking Validation")
    print("="*60)
    
    validator = CrossDockingValidator()
    
    # Prepare structures
    print("\n1️⃣ Preparing receptor structures...")
    prep_results = validator.prepare_structures()
    
    if not any(r['success'] for r in prep_results):
        print("❌ Failed to prepare any structures")
        return
    
    # Extract ligands
    print("\n2️⃣ Extracting native ligands...")
    ligand_results = validator.extract_native_ligands()
    
    if not any(r['success'] for r in ligand_results):
        print("❌ Failed to extract any ligands")
        return
    
    # Run cross-docking
    print("\n3️⃣ Performing cross-docking...")
    docking_results = validator.run_cross_docking()
    
    # Analyze results
    print("\n4️⃣ Analyzing results...")
    stats = validator.analyze_results(docking_results)
    
    if stats:
        # Generate report
        print("\n5️⃣ Generating validation report...")
        validator.generate_report(stats)
        
        print("\n✅ Cross-docking validation completed successfully!")
        print(f"📁 Results available in: {validator.output_dir}")
    else:
        print("\n❌ Failed to analyze docking results")

if __name__ == "__main__":
    main()
