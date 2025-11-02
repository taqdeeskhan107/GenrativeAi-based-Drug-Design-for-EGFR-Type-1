#!/usr/bin/env python3
"""
Validation Report Compiler
Generates comprehensive validation reports for EGFR Type 1 inhibitor dataset.
Author: TAQDEES
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')

def load_validation_results(base_dir):
    """
    Load validation results from various components.
    
    Args:
        base_dir (str): Base project directory
        
    Returns:
        dict: Combined validation results
    """
    results = {}
    
    try:
        # Load chemical validation
        chem_file = Path(base_dir) / "validation" / "chemical_validation.json"
        if chem_file.exists():
            with open(chem_file, 'r') as f:
                results['chemistry'] = json.load(f)
        
        # Load activity validation
        activity_file = Path(base_dir) / "validation" / "activity_validation.json"
        if activity_file.exists():
            with open(activity_file, 'r') as f:
                results['activity'] = json.load(f)
        
        # Load binding site validation
        binding_file = Path(base_dir) / "validation" / "binding_site_validation.json"
        if binding_file.exists():
            with open(binding_file, 'r') as f:
                results['binding_site'] = json.load(f)
        
        # Load albumin validation
        albumin_file = Path(base_dir) / "validation" / "albumin_validation.json"
        if albumin_file.exists():
            with open(albumin_file, 'r') as f:
                results['albumin'] = json.load(f)
        
    except Exception as e:
        print(f"Error loading validation results: {e}")
    
    return results

def calculate_overall_score(validation_results):
    """
    Calculate overall dataset quality score.
    
    Args:
        validation_results (dict): Combined validation results
        
    Returns:
        float: Overall score (0-100)
    """
    scores = {
        'chemistry': 0,
        'activity': 0,
        'binding_site': 0,
        'albumin': 0
    }
    
    # Chemical score (40%)
    if 'chemistry' in validation_results:
        chem = validation_results['chemistry']
        if 'valid_smiles' in chem and 'total_molecules' in chem:
            validity_score = (chem['valid_smiles'] / chem['total_molecules']) * 40
            scores['chemistry'] = validity_score
    
    # Activity score (25%)
    if 'activity' in validation_results:
        act = validation_results['activity']
        if 'potency_analysis' in act:
            potency = act['potency_analysis']
            activity_score = (
                (potency.get('highly_active_percent', 0) * 0.6) +
                (potency.get('active_percent', 0) * 0.4)
            ) * 0.25
            scores['activity'] = min(25, activity_score)
    
    # Binding site score (20%)
    if 'binding_site' in validation_results:
        site = validation_results['binding_site']
        if site.get('binding_site_valid', False):
            scores['binding_site'] = 20
    
    # Albumin optimization score (15%)
    if 'albumin' in validation_results:
        alb = validation_results['albumin']
        if alb.get('optimization_success', False):
            scores['albumin'] = 15
    
    # Calculate total score
    total_score = sum(scores.values())
    return min(100, total_score)

def generate_validation_report(base_dir, validation_results):
    """
    Generate comprehensive validation report.
    
    Args:
        base_dir (str): Base project directory
        validation_results (dict): Validation results
        
    Returns:
        str: Markdown formatted report
    """
    
    overall_score = calculate_overall_score(validation_results)
    
    report = f"""# BindingForge Dataset Validation Report

## Overall Quality Score: {overall_score:.1f}/100

### Summary
{'✅ PASSED' if overall_score >= 80 else '⚠️ NEEDS REVIEW' if overall_score >= 60 else '❌ FAILED'}

## 1. Chemical Structure Validation
"""
    
    if 'chemistry' in validation_results:
        chem = validation_results['chemistry']
        report += f"""
- **Valid Structures:** {chem.get('valid_smiles', 0)}/{chem.get('total_molecules', 0)} ({chem.get('valid_smiles', 0)/chem.get('total_molecules', 1)*100:.1f}%)
"""
        
        if 'property_validation' in chem:
            prop = chem['property_validation']
            if 'MW' in prop:
                report += f"- **MW Range:** {prop['MW']['min']:.1f} - {prop['MW']['max']:.1f} Da\n"
            if 'LogP' in prop:
                report += f"- **LogP Optimization:** {prop['LogP']['logp_le_36_percent']:.1f}% ≤ 3.6\n"
            if 'QED' in prop:
                report += f"- **Drug-like (QED ≥ 0.3):** {prop['QED']['drug_like_percent']:.1f}%\n"
        
        if 'chemical_diversity' in chem:
            div = chem['chemical_diversity']
            report += f"- **Scaffold Diversity:** {div.get('unique_scaffolds', 0)} unique scaffolds ({div.get('diversity_ratio', 0)*100:.1f}%)\n"
    
    report += "\n## 2. Activity Validation\n"
    
    if 'activity' in validation_results:
        act = validation_results['activity']
        if 'potency_analysis' in act:
            pot = act['potency_analysis']
            report += f"""
- **Highly Active (≤100 nM):** {pot.get('highly_active_percent', 0):.1f}%
- **Moderate Activity:** {pot.get('active_percent', 0):.1f}%
- **Activity Range:** {pot.get('min_activity', 0):.1f} - {pot.get('max_activity', 0):.1f} nM
"""
    
    report += "\n## 3. Binding Site Analysis\n"
    
    if 'binding_site' in validation_results:
        site = validation_results['binding_site']
        report += f"""
- **Binding Site Characterization:** {'✅ Complete' if site.get('binding_site_valid', False) else '❌ Incomplete'}
- **Key Interactions:** {site.get('key_interactions', 0)} validated
- **Pocket Volume:** {site.get('pocket_volume', 0):.1f} Å³
"""
    
    report += "\n## 4. Albumin Binding Optimization\n"
    
    if 'albumin' in validation_results:
        alb = validation_results['albumin']
        report += f"""
- **Optimization Status:** {'✅ Successful' if alb.get('optimization_success', False) else '❌ Needs improvement'}
- **Binding Reduction:** {alb.get('binding_reduction', 0):.1f}%
- **Retained Activity:** {alb.get('activity_retention', 0):.1f}%
"""
    
    report += f"""
## Recommendations

1. {'✅ Dataset ready for use' if overall_score >= 80 else '🔧 Address validation issues'}
2. {'✅ Chemical structures validated' if validation_results.get('chemistry', {}).get('valid_smiles', 0) >= validation_results.get('chemistry', {}).get('total_molecules', 1) * 0.95 else '⚠️ Review chemical structures'}
3. {'✅ Activity distribution appropriate' if validation_results.get('activity', {}).get('potency_analysis', {}).get('highly_active_percent', 0) >= 30 else '⚠️ Check activity distribution'}
4. {'✅ Binding site complete' if validation_results.get('binding_site', {}).get('binding_site_valid', False) else '❌ Complete binding site analysis'}
5. {'✅ Albumin optimization successful' if validation_results.get('albumin', {}).get('optimization_success', False) else '⚠️ Review albumin binding'}

---
Generated by BindingForge Validation Pipeline v1.0
"""
    
    return report

def save_validation_report(report, base_dir):
    """
    Save validation report to file.
    
    Args:
        report (str): Generated report content
        base_dir (str): Base directory for output
    """
    output_dir = Path(base_dir) / "documentation"
    output_dir.mkdir(exist_ok=True)
    
    report_file = output_dir / "dataset_validation_report.md"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(f"📋 Validation report saved: {report_file}")

def main():
    """Main function for report compilation"""
    print("📊 BindingForge Validation Report Compiler")
    print("="*50)
    
    base_dir = Path.cwd()
    
    # Load validation results
    print("📥 Loading validation results...")
    validation_results = load_validation_results(base_dir)
    
    if not validation_results:
        print("❌ No validation results found!")
        return
    
    # Generate report
    print("📝 Generating validation report...")
    report = generate_validation_report(base_dir, validation_results)
    
    # Save report
    print("💾 Saving validation report...")
    save_validation_report(report, base_dir)
    
    # Calculate overall score
    overall_score = calculate_overall_score(validation_results)
    
    print(f"\n📊 Overall Quality Score: {overall_score:.1f}/100")
    if overall_score >= 80:
        print("✅ Validation PASSED - Dataset ready for use!")
    elif overall_score >= 60:
        print("⚠️ Validation PARTIAL - Review recommendations")
    else:
        print("❌ Validation FAILED - Major improvements needed")

if __name__ == "__main__":
    main()
