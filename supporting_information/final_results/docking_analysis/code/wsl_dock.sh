#!/bin/bash
# Simple WSL AutoDock Vina Docking Script
# Save as: wsl_dock.sh in your WSL home directory

echo "🎯 WSL AutoDock Vina EGFR Docking"
echo "================================="

# Set paths (Windows paths in WSL format)
BASE_DIR="/mnt/c/Users/admin/BF-final-version/Combined_EGFR_Docking_Analysis"
RECEPTOR="$BASE_DIR/03_Receptor/egfr_receptor.pdbqt"
LIGANDS_DIR="$BASE_DIR/02_Structures/PDBQT_Files"
OUTPUT_DIR="$BASE_DIR/04_Docking/Vina_Outputs"
LOGS_DIR="$BASE_DIR/04_Docking/Vina_Logs"

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOGS_DIR"

# EGFR binding site coordinates (from your successful previous docking)
CENTER_X=25.0
CENTER_Y=4.0
CENTER_Z=44.0
SIZE_X=20.0
SIZE_Y=20.0
SIZE_Z=20.0

echo "📋 Docking Parameters:"
echo "   Receptor: $RECEPTOR"
echo "   Binding site: ($CENTER_X, $CENTER_Y, $CENTER_Z)"
echo "   Search box: ${SIZE_X}x${SIZE_Y}x${SIZE_Z} Å"
echo ""

# Check if receptor exists
if [ ! -f "$RECEPTOR" ]; then
    echo "❌ Receptor file not found: $RECEPTOR"
    exit 1
fi

echo "✅ Receptor found: $(basename $RECEPTOR)"

# Count ligand files
LIGAND_COUNT=$(ls "$LIGANDS_DIR"/*.pdbqt 2>/dev/null | wc -l)
echo "✅ Found $LIGAND_COUNT ligand files"

if [ $LIGAND_COUNT -eq 0 ]; then
    echo "❌ No ligand files found in: $LIGANDS_DIR"
    exit 1
fi

echo ""
echo "🚀 Starting docking..."

# Track results
SUCCESS_COUNT=0
TOTAL_COUNT=0

# Loop through all PDBQT files
for ligand_file in "$LIGANDS_DIR"/*.pdbqt; do
    if [ -f "$ligand_file" ]; then
        TOTAL_COUNT=$((TOTAL_COUNT + 1))
        mol_name=$(basename "$ligand_file" .pdbqt)
        
        echo "🧬 [$TOTAL_COUNT/$LIGAND_COUNT] Docking $mol_name..."
        
        output_file="$OUTPUT_DIR/${mol_name}_docked.pdbqt"
        log_file="$LOGS_DIR/${mol_name}_log.txt"
        
        # Run AutoDock Vina
        vina \
            --ligand "$ligand_file" \
            --receptor "$RECEPTOR" \
            --center_x $CENTER_X \
            --center_y $CENTER_Y \
            --center_z $CENTER_Z \
            --size_x $SIZE_X \
            --size_y $SIZE_Y \
            --size_z $SIZE_Z \
            --out "$output_file" \
            --log "$log_file" \
            --exhaustiveness 8 \
            --num_modes 9 \
            2>/dev/null
        
        # Check if docking was successful
        if [ $? -eq 0 ] && [ -f "$log_file" ]; then
            # Extract best score from log file
            score=$(grep -A 3 "mode |   affinity | dist from best mode" "$log_file" | tail -1 | awk '{print $2}')
            if [ ! -z "$score" ]; then
                echo "   ✅ Score: $score kcal/mol"
                SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
                echo "$mol_name,$score" >> "$BASE_DIR/05_Results/wsl_scores.csv"
            else
                echo "   ⚠️  Could not parse score"
            fi
        else
            echo "   ❌ Docking failed"
        fi
    fi
done

echo ""
echo "🎉 Docking completed!"
echo "📊 Results: $SUCCESS_COUNT/$TOTAL_COUNT successful"

# Create results summary
RESULTS_DIR="$BASE_DIR/05_Results"
mkdir -p "$RESULTS_DIR"

if [ $SUCCESS_COUNT -gt 0 ]; then
    echo ""
    echo "🏆 TOP 10 RESULTS:"
    echo "=================="
    echo "Rank | Score (kcal/mol) | Molecule"
    echo "=================================="
    
    # Sort results and show top 10
    if [ -f "$RESULTS_DIR/wsl_scores.csv" ]; then
        sort -t, -k2 -n "$RESULTS_DIR/wsl_scores.csv" | head -10 | nl -w3 -s". " | while read line; do
            rank=$(echo "$line" | awk '{print $1}')
            mol=$(echo "$line" | awk -F, '{print $1}' | awk '{print $2}')
            score=$(echo "$line" | awk -F, '{print $2}')
            printf "%3s | %13s | %s\n" "$rank" "$score" "$mol"
        done
    fi
    
    echo ""
    echo "📁 Results saved to:"
    echo "   - Docked poses: $OUTPUT_DIR"
    echo "   - Log files: $LOGS_DIR"
    echo "   - Scores: $RESULTS_DIR/wsl_scores.csv"
fi

echo ""
echo "✅ WSL AutoDock Vina docking completed!"