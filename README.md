# GWPred

**Code repository for manuscript**  
_Predicting trace element concentrations and mixtures in U.S. private drinking water using machine learning methods_
---

## 📁 Project Structure

```text
GWPred/
├── LICENSE               # License information
├── GWPred.Rproj          # RStudio project file
├── README.md             # Project documentation
├── renv.lock             # Records the exact versions of R and R packages to ensure reproducibility
├── renv/                 # Directory containing renv infrastructure for the project
├── Prediction Maps/      # Directory containing prediction maps produced by the final model
└── R/                    # Main analysis scripts
    ├── 0_helper_fct.R/            # Script to download or load raw data
    ├── 1_All_Data_Prep.R/         # Data cleaning and preprocessing
    ├── 2_Predictors_Extract.R/    # Extract predictor variables from various sources
    ├── 3_Predictor_Selection.R/   # Select important predictor variables
    ├── 4_MICE.R/                  # Multiple Imputation by Chained Equations for missing data
    ├── 5_Regression_Model.R/      # Develop predictive models
    ├── 6_Model_Eval_Reg.R/        # Evaluate model performance
    ├── 7_Shapley_Analysis.R/      # Conduct Shapley value analysis for model interpretability
    ├── 8_geoSOM_Grid.R/           # Generate geoSOM grid for mixture analysis
    ├── 9_Map_Predict.R/           # Generate prediction maps
    └── 10_Bootstrap_Predict.R/    # Estimate prediction uncertainty using bootstrap methods
```
