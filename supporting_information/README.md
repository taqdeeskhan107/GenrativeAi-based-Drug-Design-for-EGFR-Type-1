# Supporting Information for BindingForge Project

This directory contains the supporting information and supplementary data for the project, focusing on EGFR Type 1 inhibitor analysis and validation.

## Important Note on Project Naming
Please note that "BindingForge" is a working title used during the development and coding phase of this project. The published paper uses its formal academic title. The name "BindingForge" appears in code, file names, and development documentation but should not be confused with the paper's official title.

## Important Note on Molecule Numbering
For clarity and presentation in the paper, we have adopted the following molecule numbering convention:
- Original molecule #02 → Referred to as Molecule 01 in the paper
- Original molecule #03 → Referred to as Molecule 02 in the paper
- Original molecule #17 → Referred to as Molecule 03 in the paper
- Original molecule #18 → Referred to as Molecule 04 in the paper
- Original molecule #16 → Referred to as Molecule 05 in the paper

This standardized numbering (01-05) is used consistently throughout the manuscript and figures for better readability and presentation.

## Directory Structure

```
BF-final-version/
├── Data Processing/
│   ├── 01_database_setup.py           # Initial database setup
│   ├── 02_chembl_data_collection.py   # Data collection from ChEMBL
│   └── 03_process_type1_inhibitors.py # Type 1 inhibitor processing
├── Structural Analysis/
│   ├── 04_binding_site_extraction.py  # Binding site analysis
│   ├── 05_albumin_binding_analysis.py # Albumin binding studies
│   └── binding_site_data/            # Binding site extraction results
├── Docking and Validation/
│   ├── 03_Receptor/                   # Receptor preparation files
│   ├── 04_Docking/                    # Core docking results
│   ├── Combined_EGFR_Docking_Analysis/# Comprehensive docking analysis
│   ├── docking_workspace/             # Active docking workspace
│   └── validation/                    # Validation results and metrics
├── Results and Visualization/
│   ├── 06_data_visualization.py       # Visualization scripts
│   ├── figures/                       # Generated figures
│   │   ├── figure_1_vina_correlation.png
│   │   ├── figure_2_cbdock_correlation.png
│   │   └── [additional figure files]
│   └── results/                       # Analysis results
├── Model Data/
│   ├── checkpoints/                   # Model checkpoints
│   ├── training_logs/                 # Training progress logs
│   └── model_data/                    # Core model files
└── Documentation/
    ├── documentation/                 # Project documentation
    ├── PyMOL_Analysis_Instructions.md # PyMOL usage guide
    └── WSL_SETUP_GUIDE.md            # Setup instructions
```

## Figure Mapping

### Main Analysis Figures
- **Figure 1**: AutoDock Vina Correlation Analysis
  - Location: `figure_1_vina_correlation.png`
  
- **Figure 2**: CB-Dock Correlation Analysis
  - Location: `figure_2_cbdock_correlation.png`
  
- **Figure 3**: Composite Correlation Analysis
  - Location: `figure_3_composite_correlation.png`
  
- **Figure 4**: Model Performance Metrics
  - Location: `figure_4_model_performance.png`

### Prediction Analysis Figures
- **Figure 5**: Prediction Analysis Series
  - Vina Predictions: `figure_5_vina_predictions.png`
  - CB-Dock Predictions: `figure_5_cbdock_predictions.png`
  - Composite Predictions: `figure_5_composite_predictions.png`
  - Ensemble Predictions: `figure_5_ensemble_predictions.png`

### Advanced Analysis Figures
- **Figure 9**: Comprehensive Comparison
  - Location: `figure_9_comprehensive_comparison.png`
  
- **Figure 10**: Ranking Correlation Analysis
  - Location: `figure_10_ranking_correlation.png`
  
- **Figure 11**: Activity Heatmap
  - Location: `figure_11_activity_heatmap.png`
  
- **Figure 12**: Confidence Intervals
  - Location: `figure_12_confidence_intervals.png`

### Supplementary Visualizations
- **Molecular Structure**: Example visualization
  - Location: `molecule_02_structure.png`
  
- **SAR Analysis**: Structure-Activity Relationship
  - Location: `sar_heatmap.png`
  
- **Literature Analysis**: IC50 Validation
  - Location: `bindingforge_literature_ic50_analysis.png`

## Data Files and Results

### Core Dataset Files
- **IC50 Predictions**: 
  - Location: `bindingforge_ic50_predictions_jcim.csv`
  - Contents: Comprehensive IC50 predictions and analysis

### Validation Results
- **Focused Validation**: 
  - Location: `focused_validation_results.csv`
  - Contents: Detailed validation metrics for selected compounds

- **Simple Analysis**: 
  - Location: `simple_analysis_results.txt`
  - Contents: Quick analysis output and preliminary results

- **Validation Summary**: 
  - Location: `bindingforge_validation_summary.txt`
  - Contents: Overall validation metrics and conclusions

### Analysis Outputs
- **Docking Pipeline Log**: 
  - Location: `docking_pipeline.log`
  - Contents: Complete docking execution records

- **Test Logs**: 
  - Location: `test.log`, `test_log.txt`
  - Contents: Testing and validation process logs

## Model Implementation and Training

### Core Components
- **Generator Scripts**
  - `egfr_molecule_generator.py`: Main molecule generation implementation
  - `comprehensive_egfr_docking.py`: Complete docking pipeline

### Data Processing Pipeline
- **Setup and Collection**
  - `01_database_setup.py`: Database initialization
  - `02_chembl_data_collection.py`: ChEMBL data extraction
  - `03_process_type1_inhibitors.py`: Type 1 inhibitor processing

### Analysis Tools
- **Binding Analysis**
  - `04_binding_site_extraction.py`: Binding site analysis
  - `05_albumin_binding_analysis.py`: Albumin binding studies
  - `create_pymol_analysis.py`: PyMOL visualization generation

### Validation Suite
- **Core Validation**
  - `07_dataset_validation.py`: Comprehensive validation
  - `quick_validation_analysis.py`: Rapid testing module

### Model Resources
- **Configuration**
  - `tokenizer.json`: Molecule tokenization configuration
  - Located in project root

### Training Data
- **Logs and Checkpoints**
  - Location: `training_logs/`: Training process logs
  - Location: `checkpoints/`: Model checkpoints

### Results
- **Analysis Outputs**
  - Location: `results/`: Generated results
  - Location: `visualizations/`: Visual analysis outputs

## Key Validation Files

### Docking Analysis
- **Erlotinib Validation**
  - Location: `final_results/docking_analysis/05_Results/erlotinib_validation.txt`
  - Contains: RMSD validation (4.41 Å) and docking protocol details
  
- **Docking Results**
  - Location: `final_results/docking_analysis/05_Results/final_docking_results.csv`
  - Contains: Binding affinity scores (-7.7 to -8.4 kcal/mol)

### Machine Learning Metrics
- **Model Validation**
  - Location: `final_results/docking_analysis/05_Results/model_validation_metrics.csv`
  - Contains:
    - Training loss: 0.142
    - Validation loss: 0.163
    - Test loss: 0.171
    - Generation metrics:
      - Validity: 82.6%
      - Uniqueness: 28.3%
      - Novelty: 79.5%

### Structural Analysis
- **PLIP Analysis**
  - Location: `final_results/plip_analysis/`
  - Files: PLIPMOL1.png through PLIPMOL5.png
  - Contains: Protein-ligand interaction profiles
  
- **Cavity Analysis**
  - Location: `final_results/cavity_analysis/`
  - Contains: Binding site characterization for all lead molecules

## Validation Scripts
The validation functionality is distributed across several core files:

- `07_dataset_validation.py`: Comprehensive dataset validation including RMSD calculations and validation reports
- `combined_docking_pipeline.py`: Complete docking pipeline with receptor preparation and cross-docking validation
- `04_binding_site_extraction.py`: Binding site and co-crystal ligand extraction
- `quick_validation_analysis.py`: Simplified validation for rapid testing

## Dataset Files
- **Raw Data**
  - Location: `dataset/raw_data/egfr_raw_chembl.csv`
  - Contains: Original ChEMBL extractions
  
- **Processed Data**
  - Location: `dataset/processed_data/egfr_type1_filtered.csv`
  - Contains: Filtered and processed EGFR inhibitors

