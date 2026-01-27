# synthetic_dynasty_analysis_replicate

## Class Mobility in the Era of Rising Inequality

### A Synthetic Dynasty Analysis

**Authors:**
Geoffrey T. Wodtke, Weiqi Wang, Kristina Butaeva, and Steven Durlauf

This repository contains the replication code for the paper:

> **Class Mobility in the Era of Rising Inequality: A Synthetic Dynasty Analysis**

It provides all scripts necessary to reproduce the **baseline empirical results and figures locally**, as well as the codes for replicating the **full set of results using high-performance computing (HPC)** resources.

> **Disclaimer:** All remaining faults are the responsibility of **Weiqi Wang**.

---

## Repository Structure

The repository is organized as follows:

* `Functions/`
  Core function definitions (estimation routines, plotting utilities, data processing helpers)

* `Data/`
  Contains supporting data files, including a local copy of the Morgan (2017) occupational crosswalk

* **Main directory scripts** (run in order to replicate baseline results):

  * `00_generate_main_results_data.R`
  * `01_generate_main_results_estimation.R`
  * `02_create_fig1.R`
  * …
  * `10_create_fig9.R`

  These scripts are ordered sequentially to allow full replication of the **main results and figures**.

* `Robust_Codes/`
  Scripts for **robustness analyses**, organized by topic (e.g. age, gender, race, alternative measures, smoothing choices, class typologies, etc.)

* `sample.sh`
  A sample bash script with suggested computational parameters for running large jobs on HPC systems (e.g. SLURM environments)

---

## Data Sources

This project relies on two external data sources:

* **GSS 2024 data**
  The data are fetched programmatically within the scripts using the R package `gssr`.
  Users do **not** need to manually download the GSS data.

* **Morgan (2017) occupational crosswalk**
  Available at: [https://osf.io/9nkrw/overview](https://osf.io/9nkrw/overview)

  For convenience and reproducibility, a local copy of the crosswalk is also included in the `Data/` folder.

---

## Reproducing Baseline Results Locally

The files included in this repository are sufficient to reproduce:

* All **baseline estimation results**
* On a standard desktop or laptop environment

To replicate the main results, run the master scripts in order:

* `master_main_results.R`
* `master_robust_ABCD.R`
* `master_robust_EFGH.R`


Typical runtime is approximately **2–4 hours per major estimation block**, depending on hardware.

---

## Full Replication (Recommended: HPC)

Running the full set of robustness analyses can be computationally intensive.
We therefore recommend using **high-performance computing (HPC)** resources.

A sample SLURM-compatible submission script is provided:

* `sample.sh`

This file includes suggested settings for memory, CPU cores, and runtime.

---

## Acknowledgements

* Some plotting functions were developed with reference to outputs generated using **ChatGPT-4o**.
* Functions related to **Hellinger’s dependence measure** are adapted from code provided by **Prof. Thomas Coleman**, whom the authors thank for introducing valuable perspectives on dependence measures.

---

##  Contact

For questions regarding the code or replication:

**Weiqi Wang**
📧 [weiqiw@uchicago.edu](mailto:weiqiw@uchicago.edu)

