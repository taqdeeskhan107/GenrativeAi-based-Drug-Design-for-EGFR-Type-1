#!/usr/bin/env python3
"""
Simplified Validation Script
Quick validation checks for EGFR Type 1 inhibitor dataset.
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

# Optional RDKit import for chemical validation
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
    print("✅ RDKit available for chemical validation")
except ImportError:
    print("⚠️ RDKit not available. Limited validation capabilities.")
    RDKIT_AVAILABLE = False

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class SimpleValidator:
    """Quick validation checks for dataset quality."""
    
    def __init__(self, dataset_path):
        """Initialize validator with dataset path."""
        self.dataset_path = Path(dataset_path)
        self.results = {
            'timestamp': datetime.now().isoformat(),
            'dataset': str(dataset_path),
            'validation_summary': {},
            'warnings': [],
            'status': 'pending'
        }
    
    def check_data_completeness(self, df):
        """Check for missing values and data completeness."""
        
        total_rows = len(df)
        results = {
            'total_entries': total_rows,
            'missing_data': {}
        }
        
        # Check each column for missing values
        for column in df.columns:
            missing = df[column].isna().sum()
            if missing > 0:
                missing_pct = (missing / total_rows) * 100
                results['missing_data'][column] = {
                    'count': int(missing),
                    'percentage': float(missing_pct)
                }
                if missing_pct > 10:  # More than 10% missing
                    self.results['warnings'].append(
                        f"High missing data in {column}: {missing_pct:.1f}%"
                    )
        
        return results
    
    def validate_chemical_structures(self, df):
        """Basic chemical structure validation."""
        
        if not RDKIT_AVAILABLE or 'smiles' not in df.columns:
            return None
        
        results = {
            'total_structures': len(df),
            'valid_structures': 0,
            'invalid_smiles': []
        }
        
        # Validate SMILES strings
        for idx, smiles in enumerate(df['smiles']):
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    results['valid_structures'] += 1
                else:
                    results['invalid_smiles'].append(idx)
            except:
                results['invalid_smiles'].append(idx)
        
        validity_ratio = results['valid_structures'] / results['total_structures']
        if validity_ratio < 0.95:
            self.results['warnings'].append(
                f"Low structure validity: {validity_ratio*100:.1f}%"
            )
        
        return results
    
    def check_activity_data(self, df):
        """Basic validation of activity data."""
        
        if 'standard_value' not in df.columns:
            return None
        
        results = {
            'total_activities': len(df),
            'valid_activities': 0,
            'activity_range': {},
            'potency_distribution': {}
        }
        
        # Filter valid activities
        valid_mask = df['standard_value'].notna() & (df['standard_value'] > 0)
        valid_activities = df.loc[valid_mask, 'standard_value']
        results['valid_activities'] = len(valid_activities)
        
        if len(valid_activities) > 0:
            # Activity statistics
            results['activity_range'] = {
                'min': float(valid_activities.min()),
                'max': float(valid_activities.max()),
                'mean': float(valid_activities.mean()),
                'median': float(valid_activities.median())
            }
            
            # Potency distribution
            potent = (valid_activities <= 100).sum()
            moderate = ((valid_activities > 100) & (valid_activities <= 1000)).sum()
            weak = (valid_activities > 1000).sum()
            
            results['potency_distribution'] = {
                'potent_le_100nM': int(potent),
                'moderate_100_1000nM': int(moderate),
                'weak_gt_1000nM': int(weak)
            }
            
            # Check for concerning distributions
            if potent / len(valid_activities) < 0.2:  # Less than 20% potent compounds
                self.results['warnings'].append(
                    "Low proportion of potent compounds (<100 nM)"
                )
        
        return results
    
    def check_property_ranges(self, df):
        """Check for unusual property ranges."""
        
        if not RDKIT_AVAILABLE:
            return None
        
        results = {
            'molecular_weight': {},
            'logp': {},
            'out_of_range': []
        }
        
        # Property filters
        mw_range = (200, 600)  # Typical drug-like range
        logp_range = (-1, 5)   # Typical drug-like range
        
        for idx, row in df.iterrows():
            try:
                mol = Chem.MolFromSmiles(row['smiles'])
                if mol is not None:
                    mw = Descriptors.ExactMolWt(mol)
                    logp = Descriptors.MolLogP(mol)
                    
                    # Check ranges
                    if not (mw_range[0] <= mw <= mw_range[1]):
                        results['out_of_range'].append({
                            'index': idx,
                            'property': 'MW',
                            'value': float(mw)
                        })
                    
                    if not (logp_range[0] <= logp <= logp_range[1]):
                        results['out_of_range'].append({
                            'index': idx,
                            'property': 'LogP',
                            'value': float(logp)
                        })
            except:
                continue
        
        if results['out_of_range']:
            self.results['warnings'].append(
                f"Found {len(results['out_of_range'])} compounds with unusual properties"
            )
        
        return results
    
    def run_validation(self):
        """Run all validation checks."""
        
        logger.info(f"🔍 Validating dataset: {self.dataset_path}")
        
        try:
            # Load dataset
            df = pd.read_csv(self.dataset_path)
            logger.info(f"📊 Loaded {len(df)} entries")
            
            # Run checks
            validation_results = {
                'completeness': self.check_data_completeness(df),
                'chemistry': self.validate_chemical_structures(df),
                'activity': self.check_activity_data(df),
                'properties': self.check_property_ranges(df)
            }
            
            self.results['validation_summary'] = validation_results
            self.results['status'] = 'completed'
            
            # Print summary
            self.print_summary()
            
            return len(self.results['warnings']) == 0
            
        except Exception as e:
            logger.error(f"❌ Validation failed: {e}")
            self.results['status'] = 'failed'
            self.results['error'] = str(e)
            return False
    
    def print_summary(self):
        """Print validation summary to console."""
        
        print("\n" + "="*60)
        print("📊 Quick Validation Summary")
        print("="*60)
        
        if self.results['status'] == 'completed':
            summary = self.results['validation_summary']
            
            # Data completeness
            if 'completeness' in summary:
                print(f"\n📝 Data Completeness:")
                print(f"   Total entries: {summary['completeness']['total_entries']}")
                if summary['completeness']['missing_data']:
                    print("   Missing data in columns:")
                    for col, stats in summary['completeness']['missing_data'].items():
                        print(f"   - {col}: {stats['count']} ({stats['percentage']:.1f}%)")
            
            # Chemical structures
            if summary.get('chemistry'):
                chem = summary['chemistry']
                valid_pct = (chem['valid_structures'] / chem['total_structures']) * 100
                print(f"\n🧪 Chemical Structures:")
                print(f"   Valid structures: {chem['valid_structures']}/{chem['total_structures']} ({valid_pct:.1f}%)")
            
            # Activity data
            if summary.get('activity'):
                act = summary['activity']
                if 'potency_distribution' in act:
                    dist = act['potency_distribution']
                    print(f"\n🎯 Activity Distribution:")
                    print(f"   Potent (≤100 nM): {dist['potent_le_100nM']}")
                    print(f"   Moderate (100-1000 nM): {dist['moderate_100_1000nM']}")
                    print(f"   Weak (>1000 nM): {dist['weak_gt_1000nM']}")
            
            # Warnings
            if self.results['warnings']:
                print("\n⚠️ Validation Warnings:")
                for warning in self.results['warnings']:
                    print(f"   • {warning}")
            else:
                print("\n✅ No validation warnings")
        else:
            print(f"\n❌ Validation failed: {self.results.get('error', 'Unknown error')}")
        
        print("\n" + "="*60)
    
    def save_results(self, output_file):
        """Save validation results to JSON file."""
        
        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        logger.info(f"💾 Results saved to: {output_file}")

def main():
    """Main function for command-line usage."""
    
    parser = argparse.ArgumentParser(
        description="Quick validation of EGFR Type 1 inhibitor dataset."
    )
    parser.add_argument(
        "dataset",
        help="Path to dataset CSV file"
    )
    parser.add_argument(
        "--output",
        help="Output file for validation results (JSON)",
        default="validation_results.json"
    )
    
    args = parser.parse_args()
    
    # Run validation
    validator = SimpleValidator(args.dataset)
    success = validator.run_validation()
    
    # Save results if specified
    if args.output:
        validator.save_results(args.output)
    
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
