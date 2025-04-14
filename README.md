#  Matrix completion with Missing Structurally and Sporadically (Macomss)

This repository supports the paper **"Integrated Analysis for Electronic Health Records with Structured and Sporadic Missingness"** by Jianbin Tan, Yan Zhang, Chuan Hong, T. Tony Cai, Tianxi Cai, and Anru R. Zhang. 


## Abstract
This repository contains code and additional resources supporting our study.

## 1. Data
- **Electronic Health Records (EHRs)** from **Duke University Health System (DUHS)**.
- Data accessed through **Duke Clinical Research Datamart (CRDM)** (Hurst et al., 2021).
- Includes data from three hospitals, each treated as an **isolated site**:
  - **Duke Raleigh Hospital (DRAH)** – Site 1
  - **Duke Regional Hospital (DRH)** – Site 2
  - **Duke University Hospital (DUH)** – Site 3


## 2. Code Structure
### Overview
Our study introduces a imputation algorithm to impute matrix with both sporadic and structure missing.

### Reproducibility
- The folder **`Results for Recovery Accuracy`**: Contains code for simulation experiments for Figure 2(A) and 2(B).
- The folder **`Row_increased_case`**: Contains code for simulation experiments for Figure 3 .
- The folder **`Feature_increased_case`**: Contains code for simulation experiments for Figure 4.
- The folder **`Real_Data_Experiment`**: Contains code for real data experiment in Section 4.3 Figure 5.
- The folder **`Heterogeneity_Experiment`**: Contains code for  real data experiment in Section 4.3 Figure 6.

Refer to the `readme.txt` files in each folder for detailed usage instructions.
