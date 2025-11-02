#!/bin/bash
# EGFR Receptor Preparation Script
# Download and prepare 1M17 structure for docking

echo "🧬 Preparing EGFR receptor (1M17)..."

# Download PDB structure
cd Combined_EGFR_Docking_Analysis\03_Receptor
wget -O 1M17.pdb "https://files.rcsb.org/download/1M17.pdb"

# Clean structure (remove waters, ligands)
grep "^ATOM" 1M17.pdb | grep " A " > 1M17_clean.pdb

# Convert to PDBQT (requires AutoDockTools)
echo "Converting to PDBQT format..."
prepare_receptor4.py -r 1M17_clean.pdb -o 1M17_receptor.pdbqt

echo "✅ Receptor preparation complete!"
echo "Files created:"
echo "  - 1M17.pdb (original structure)"
echo "  - 1M17_clean.pdb (cleaned structure)"  
echo "  - 1M17_receptor.pdbqt (ready for docking)"
