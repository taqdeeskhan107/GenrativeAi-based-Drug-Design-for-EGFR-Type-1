#!/usr/bin/env python3
"""
Main Validation Pipeline
Runs comprehensive validation of EGFR Type 1 inhibitor dataset and models.
Author: TAQDEES
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import logging
import argparse

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class ValidationPipeline:
    """Comprehensive validation pipeline for EGFR Type 1 inhibitor dataset."""
    
    def __init__(self, workspace_dir):
        """
        Initialize validation pipeline.
        
        Args:
            workspace_dir (str): Base workspace directory containing all data
        """
        self.workspace_dir = Path(workspace_dir)
        self.results_dir = self.workspace_dir / "validation_results"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # Configuration
        self.config = {
            'data_checks': True,
            'chemical_validation': True,
            'activity_validation': True,
            'structural_validation': True,
            'albumin_validation': True,
            'model_validation': True
        }
        
        # Results storage
        self.validation_results = {
            'timestamp': datetime.now().isoformat(),
            'overall_status': 'pending',
            'steps_completed': [],
            'steps_failed': [],
            'results': {}
        }
    
    def validate_data_structure(self):
        """Validate presence and structure of all required data files."""
        
        logger.info("🔍 Validating data structure...")
        
        required_files = {
            'dataset': self.workspace_dir / "processed_data/egfr_type1_filtered.csv",
            'binding_site': self.workspace_dir / "binding_site_data/binding_site_features.csv",
            'albumin_data': self.workspace_dir / "albumin_binding_analysis/egfr_type1_albumin_binding.csv",
            'model_predictions': self.workspace_dir / "model_predictions/bindingforge_ic50_predictions.csv"
        }
        
        results = {
            'valid': True,
            'missing_files': [],
            'file_stats': {}
        }
        
        for file_type, file_path in required_files.items():
            if file_path.exists():
                try:
                    stats = file_path.stat()
                    results['file_stats'][file_type] = {
                        'size': stats.st_size,
                        'modified': datetime.fromtimestamp(stats.st_mtime).isoformat()
                    }
                    logger.info(f"✅ Found {file_type}: {file_path.name}")
                except Exception as e:
                    results['valid'] = False
                    results['missing_files'].append(f"{file_type} (error: {e})")
                    logger.error(f"❌ Error checking {file_type}: {e}")
            else:
                results['valid'] = False
                results['missing_files'].append(file_type)
                logger.error(f"❌ Missing {file_type}: {file_path}")
        
        self.validation_results['results']['data_structure'] = results
        if results['valid']:
            self.validation_results['steps_completed'].append('data_structure')
        else:
            self.validation_results['steps_failed'].append('data_structure')
        
        return results['valid']
    
    def validate_chemical_data(self):
        """Validate chemical structures and properties."""
        
        if not self.config['chemical_validation']:
            logger.info("⏩ Skipping chemical validation")
            return True
        
        logger.info("🧪 Validating chemical structures...")
        
        try:
            # Import RDKit for chemical validation
            from rdkit import Chem
            from rdkit.Chem import Descriptors, QED
            
            # Load dataset
            dataset_file = self.workspace_dir / "processed_data/egfr_type1_filtered.csv"
            df = pd.read_csv(dataset_file)
            
            results = {
                'valid': True,
                'total_compounds': len(df),
                'valid_structures': 0,
                'property_stats': {},
                'alerts': []
            }
            
            # Validate SMILES
            valid_mols = []
            for smiles in df['smiles']:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    results['valid_structures'] += 1
                    valid_mols.append(mol)
            
            # Calculate validity percentage
            validity_ratio = results['valid_structures'] / results['total_compounds']
            results['structure_validity'] = validity_ratio
            
            if validity_ratio < 0.95:  # Less than 95% valid structures
                results['valid'] = False
                results['alerts'].append(f"Low structure validity: {validity_ratio*100:.1f}%")
            
            # Analyze properties
            if valid_mols:
                mw_values = [Descriptors.ExactMolWt(mol) for mol in valid_mols]
                logp_values = [Descriptors.MolLogP(mol) for mol in valid_mols]
                qed_values = [QED.default(mol) for mol in valid_mols]
                
                results['property_stats'] = {
                    'molecular_weight': {
                        'min': float(np.min(mw_values)),
                        'max': float(np.max(mw_values)),
                        'mean': float(np.mean(mw_values)),
                        'std': float(np.std(mw_values))
                    },
                    'logp': {
                        'min': float(np.min(logp_values)),
                        'max': float(np.max(logp_values)),
                        'mean': float(np.mean(logp_values)),
                        'std': float(np.std(logp_values))
                    },
                    'qed': {
                        'min': float(np.min(qed_values)),
                        'max': float(np.max(qed_values)),
                        'mean': float(np.mean(qed_values)),
                        'std': float(np.std(qed_values))
                    }
                }
                
                # Check property ranges
                if np.mean(mw_values) > 500:
                    results['alerts'].append("High average molecular weight")
                if np.mean(logp_values) > 3.6:
                    results['alerts'].append("High average LogP")
                if np.mean(qed_values) < 0.3:
                    results['alerts'].append("Low average drug-likeness")
            
            logger.info(f"✅ Validated {results['valid_structures']} structures ({validity_ratio*100:.1f}%)")
            if results['alerts']:
                for alert in results['alerts']:
                    logger.warning(f"⚠️ {alert}")
            
            self.validation_results['results']['chemical_validation'] = results
            if results['valid']:
                self.validation_results['steps_completed'].append('chemical_validation')
            else:
                self.validation_results['steps_failed'].append('chemical_validation')
            
            return results['valid']
            
        except ImportError:
            logger.error("❌ RDKit not available for chemical validation")
            self.validation_results['steps_failed'].append('chemical_validation')
            return False
        except Exception as e:
            logger.error(f"❌ Error in chemical validation: {e}")
            self.validation_results['steps_failed'].append('chemical_validation')
            return False
    
    def validate_activity_data(self):
        """Validate biological activity data."""
        
        if not self.config['activity_validation']:
            logger.info("⏩ Skipping activity validation")
            return True
        
        logger.info("🎯 Validating activity data...")
        
        try:
            # Load dataset
            dataset_file = self.workspace_dir / "processed_data/egfr_type1_filtered.csv"
            df = pd.read_csv(dataset_file)
            
            results = {
                'valid': True,
                'total_compounds': len(df),
                'activity_stats': {},
                'alerts': []
            }
            
            # Check activity values
            if 'standard_value' in df.columns:
                valid_activities = df['standard_value'].notna()
                activity_values = df.loc[valid_activities, 'standard_value']
                
                results['activity_stats'] = {
                    'total_activities': len(activity_values),
                    'valid_ratio': len(activity_values) / len(df),
                    'min_activity': float(activity_values.min()),
                    'max_activity': float(activity_values.max()),
                    'mean_activity': float(activity_values.mean()),
                    'median_activity': float(activity_values.median())
                }
                
                # Activity distribution analysis
                potent = (activity_values <= 100).sum()  # IC50 ≤ 100 nM
                moderate = ((activity_values > 100) & (activity_values <= 1000)).sum()
                weak = (activity_values > 1000).sum()
                
                results['activity_distribution'] = {
                    'potent_inhibitors': int(potent),
                    'moderate_inhibitors': int(moderate),
                    'weak_inhibitors': int(weak),
                    'potent_ratio': float(potent / len(activity_values))
                }
                
                # Validation checks
                if results['activity_stats']['valid_ratio'] < 0.9:
                    results['valid'] = False
                    results['alerts'].append("Low activity data completeness")
                
                if results['activity_distribution']['potent_ratio'] < 0.3:
                    results['alerts'].append("Low proportion of potent inhibitors")
                
                logger.info(f"✅ Validated {len(activity_values)} activity values")
                logger.info(f"📊 Potent inhibitors: {potent} ({results['activity_distribution']['potent_ratio']*100:.1f}%)")
            else:
                results['valid'] = False
                results['alerts'].append("Missing activity data column")
            
            self.validation_results['results']['activity_validation'] = results
            if results['valid']:
                self.validation_results['steps_completed'].append('activity_validation')
            else:
                self.validation_results['steps_failed'].append('activity_validation')
            
            return results['valid']
            
        except Exception as e:
            logger.error(f"❌ Error in activity validation: {e}")
            self.validation_results['steps_failed'].append('activity_validation')
            return False
    
    def validate_structural_data(self):
        """Validate structural data and binding site analysis."""
        
        if not self.config['structural_validation']:
            logger.info("⏩ Skipping structural validation")
            return True
        
        logger.info("🧬 Validating structural data...")
        
        try:
            binding_site_file = self.workspace_dir / "binding_site_data/binding_site_features.csv"
            
            results = {
                'valid': True,
                'binding_site_stats': {},
                'alerts': []
            }
            
            if binding_site_file.exists():
                df_binding = pd.read_csv(binding_site_file)
                
                results['binding_site_stats'] = {
                    'total_residues': len(df_binding),
                    'unique_chains': df_binding['chain'].nunique(),
                    'residue_types': df_binding['residue_name'].nunique()
                }
                
                # Validate expected number of binding site residues (typically 15-25)
                if not (15 <= results['binding_site_stats']['total_residues'] <= 25):
                    results['alerts'].append("Unusual number of binding site residues")
                
                # Analyze residue composition if available
                if 'category' in df_binding.columns:
                    composition = df_binding['category'].value_counts()
                    results['binding_site_stats']['composition'] = composition.to_dict()
                    
                    # Check for essential residue types
                    required_types = {'Charged', 'Polar', 'Hydrophobic'}
                    missing_types = required_types - set(composition.index)
                    if missing_types:
                        results['alerts'].append(f"Missing residue types: {missing_types}")
                
                logger.info(f"✅ Validated binding site: {results['binding_site_stats']['total_residues']} residues")
            else:
                results['valid'] = False
                results['alerts'].append("Missing binding site data file")
            
            self.validation_results['results']['structural_validation'] = results
            if results['valid']:
                self.validation_results['steps_completed'].append('structural_validation')
            else:
                self.validation_results['steps_failed'].append('structural_validation')
            
            return results['valid']
            
        except Exception as e:
            logger.error(f"❌ Error in structural validation: {e}")
            self.validation_results['steps_failed'].append('structural_validation')
            return False
    
    def run_validation(self):
        """Run complete validation pipeline."""
        
        logger.info("\n🔄 Starting validation pipeline...")
        start_time = datetime.now()
        
        # Run validation steps
        data_valid = self.validate_data_structure()
        if not data_valid and not self.config.get('continue_on_error', False):
            logger.error("❌ Data structure validation failed. Stopping pipeline.")
            return False
        
        chemical_valid = self.validate_chemical_data()
        activity_valid = self.validate_activity_data()
        structural_valid = self.validate_structural_data()
        
        # Generate summary
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        self.validation_results.update({
            'duration_seconds': duration,
            'completed_at': end_time.isoformat(),
            'overall_status': 'passed' if all([data_valid, chemical_valid, activity_valid, structural_valid]) else 'failed'
        })
        
        # Save results
        self.save_validation_report()
        
        # Print summary
        self.print_validation_summary()
        
        return self.validation_results['overall_status'] == 'passed'
    
    def save_validation_report(self):
        """Save validation report to file."""
        
        report_file = self.results_dir / f"validation_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        
        with open(report_file, 'w') as f:
            json.dump(self.validation_results, f, indent=2)
        
        logger.info(f"📝 Validation report saved to: {report_file}")
    
    def print_validation_summary(self):
        """Print validation summary to console."""
        
        print("\n" + "="*60)
        print("📊 Validation Pipeline Summary")
        print("="*60)
        
        print(f"\n🏁 Overall Status: {self.validation_results['overall_status'].upper()}")
        print(f"⏱️ Duration: {self.validation_results['duration_seconds']:.1f} seconds")
        
        print("\n✅ Completed Steps:")
        for step in self.validation_results['steps_completed']:
            print(f"   • {step}")
        
        if self.validation_results['steps_failed']:
            print("\n❌ Failed Steps:")
            for step in self.validation_results['steps_failed']:
                print(f"   • {step}")
        
        if self.validation_results['results'].get('chemical_validation', {}).get('alerts'):
            print("\n⚠️ Chemical Alerts:")
            for alert in self.validation_results['results']['chemical_validation']['alerts']:
                print(f"   • {alert}")
        
        if self.validation_results['results'].get('activity_validation', {}).get('alerts'):
            print("\n⚠️ Activity Alerts:")
            for alert in self.validation_results['results']['activity_validation']['alerts']:
                print(f"   • {alert}")
        
        print("\n" + "="*60)

def main():
    """Main function for command-line usage."""
    
    parser = argparse.ArgumentParser(
        description="Run comprehensive validation pipeline for EGFR Type 1 inhibitor dataset."
    )
    parser.add_argument(
        "--workspace",
        default=".",
        help="Base workspace directory containing all data (default: current directory)"
    )
    parser.add_argument(
        "--skip-chemical",
        action="store_true",
        help="Skip chemical structure validation"
    )
    parser.add_argument(
        "--skip-activity",
        action="store_true",
        help="Skip activity data validation"
    )
    parser.add_argument(
        "--skip-structural",
        action="store_true",
        help="Skip structural data validation"
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue validation even if some steps fail"
    )
    
    args = parser.parse_args()
    
    # Initialize and configure pipeline
    pipeline = ValidationPipeline(args.workspace)
    
    pipeline.config.update({
        'chemical_validation': not args.skip_chemical,
        'activity_validation': not args.skip_activity,
        'structural_validation': not args.skip_structural,
        'continue_on_error': args.continue_on_error
    })
    
    # Run validation
    success = pipeline.run_validation()
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
