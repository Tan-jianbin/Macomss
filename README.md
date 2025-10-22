# Matrix Completion with Structured and Sporadic Missingness (Macomss)

This codebase is developed by **Yan Zhang**, **Anru Zhang**, and **Jianbin Tan**. The repository accompanies the paper **“Integrated Analysis for Electronic Health Records with Structured and Sporadic Missingness”** by **Jianbin Tan**, **Yan Zhang**, **Chuan Hong**, **T.Tony Cai**, **Tianxi Cai**, and **Anru Zhang** ([ArXiv:2506.09208](https://arxiv.org/abs/2506.09208)). The paper is published in [Journal of Biomedical Informatics](https://www.sciencedirect.com/science/article/pii/S1532046425001625?via%3Dihub).

## Citation
If you find our code, data, or methodology useful, you may cite us as:

    @article{TAN2025104933,
      title   = {Integrated analysis for electronic health records with structured and sporadic missingness},
      author  = {Jianbin Tan and Yan Zhang and Chuan Hong and T. Tony Cai and Tianxi Cai and Anru R. Zhang},
      journal = {Journal of Biomedical Informatics},
      volume  = {171},
      pages   = {104933},
      year    = {2025},
      doi     = {https://doi.org/10.1016/j.jbi.2025.104933}
    }

---

## Overview

This repository provides the implementation of a novel matrix completion algorithm tailored for data with **both structured and sporadic missingness**, as encountered in multi-site electronic health record (EHR) systems. The method supports robust imputation across heterogeneous healthcare systems by leveraging structural information and cross-site similarities.

## 1. Data Description

- The analysis is based on **Electronic Health Records (EHRs)** from the **Duke University Health System (DUHS)**.
- Data were obtained through the **Duke Clinical Research Datamart (CRDM)**.
- The dataset includes records from three hospitals, each treated as an **independent site**:
  - **Duke Raleigh Hospital (DRAH)** – Site 1  
  - **Duke Regional Hospital (DRH)** – Site 2  
  - **Duke University Hospital (DUH)** – Site 3  

## 2. Code Structure

### Algorithm Implementation
The repository implements our proposed imputation framework for matrix recovery under mixed missingness patterns.

### Reproducibility
To facilitate replication of results, the code is organized as follows:

- **`Results_for_Recovery_Accuracy/`**  
  Contains simulation scripts for generating Figure 2

- **`Row_increased_case/`**  
  Scripts for row-wise scalability experiments in Figure 3.

- **`Feature_increased_case/`**  
  Scripts for feature-wise scalability experiments in Figure 4.
  
  For detailed usage and instructions, please refer to the `readme` files included in each subfolder.