# Matrix Completion with Structured and Sporadic Missingness (Macomss)

This codebase is developed by Anru R. Zhang, Yan Zhang, and Jianbin Tan. The repository accompanies the paper **"Integrated Analysis for Electronic Health Records with Structured and Sporadic Missingness"** by **Jianbin Tan**, **Yan Zhang**, **Chuan Hong**, **T. Tony Cai**, **Tianxi Cai**, and **Anru R. Zhang**.

## Overview

This repository provides the implementation of a novel matrix completion algorithm tailored for data with **both structured and sporadic missingness**, as encountered in multi-site electronic health record (EHR) systems. The method supports robust imputation across heterogeneous healthcare systems by leveraging structural information and cross-site similarities.

## 1. Data Description

- The analysis is based on **Electronic Health Records (EHRs)** from the **Duke University Health System (DUHS)**.
- Data were obtained through the **Duke Clinical Research Datamart (CRDM)** (Hurst et al., 2021).
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
  
  For detailed usage and instructions, please refer to the individual `readme.txt` files included in each subfolder.