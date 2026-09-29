# synthetic_dynasty_analysis_replicate

## Class Mobility in the Era of Rising Inequality

### A Synthetic Dynasty Analysis

**Authors:**
Geoffrey T. Wodtke, Weiqi Wang, Kristina Butaeva, and Steven Durlauf

This repository contains the replication code for the paper:

> **Class Mobility in the Era of Rising Inequality: A Synthetic Dynasty Analysis**

It provides all scripts necessary to reproduce the **baseline empirical results and figures locally**, as well as the codes for replicating the **full set of results using high-performance computing (HPC)** resources.


---

## Repository Structure

The repository is organized as follows:

* `Functions/`
  Core function definitions (estimation routines, plotting utilities, data processing helpers)

* `NLSY_fun/`
  Function definitions, configuration, and the figure-data verification manifest for the NLSY validation pipeline

* `perturbation_fun/`
  Function definitions for the sparsity perturbation analysis (inputs, mechanisms, measures, bootstrap trends, figures, tables)

* `Data/`
  Contains supporting data files, including a local copy of the Morgan (2017) occupational crosswalk.
  The NLSY input data are hosted externally (see **Data Sources** below) and should be placed under `Data/NLSY/` before running the NLSY pipeline.

* **Main directory scripts** (run in order to replicate baseline results):

  * `00_generate_main_results_data.R`
  * `01_generate_main_results_estimation.R`
  * `02_create_fig1.R`
  * …
  * `10_create_fig9.R`
  * `11_create_figD1_D2.R`

  These scripts are ordered sequentially to allow full replication of the **main results and figures** (Figures 1–9 and Appendix Figures D.1–D.2).

* `Robust_Codes/`
  Scripts for **robustness analyses**, organized by topic. Script prefixes and figure names follow the appendix lettering of the current manuscript:

  | Folder | Prefix | Manuscript appendix |
  |---|---|---|
  | `Age/` | `RA` | A. Alternative Age Restrictions |
  | `Class_Typology/` | `RB` | B. Alternative Class Typologies |
  | `Weighted/` | `RC` | C. Weighted Estimates |
  | `Alternative_Measure/` | `RE` | E. Alternative Measures of Mobility |
  | `Gender/` | `RF1`–`RF14` | F. Trends by Gender (father/mother–child, son, daughter, higher-status parent) |
  | `Education/` | `RF15`–`RF24` | F. Trends by Gender (educational attainment) |
  | `Race/` | `RG` | G. Trends by Race |
  | `Alternative_Smoothing/` | `RH` | H. Alternative Approaches to Smoothing |
  | `NLSY/` | `00`–`08` | I. Income and Occupational Mobility in the NLSY |
  | `Extended_Cohorts/` | `RJ` | J. An Extended Series of Cohorts |
  | `Perturbation/` | `01`–`03` | K. Sparsity and Estimation Uncertainty |

  (Appendix D, Confidence Intervals for Class-Specific Measures, is produced by the main-directory script `11_create_figD1_D2.R`.)

* `sample.sh`
  A sample bash script with suggested computational parameters for running large jobs on HPC systems (e.g. SLURM environments)

---

## Data Sources

This project relies on the following external data sources:

* **GSS 1972–2024 Cross-Sectional Cumulative Data, Release 2**
  All GSS-based analyses use **Release 2** of the 1972–2024 cumulative file
  (NORC, October 2025; Stata file `gss7224_r2.dta`).
  When run locally, the data are fetched programmatically via the R package
  `gssr`, which bundles this release as of its October 2025 update; users do
  **not** need to manually download the GSS data. On HPC systems, place
  `gss7224_r2.dta` under `Data/`. Note that NORC has since issued later
  releases (e.g., Release 3, March 2026); exact replication requires
  Release 2.

* **Morgan (2017) occupational crosswalk**
  Available at: [https://osf.io/9nkrw/overview](https://osf.io/9nkrw/overview)

  For convenience and reproducibility, a local copy of the crosswalk is also included in the `Data/` folder.

* **NLSY data for the income and occupational mobility analyses (Appendix I)**
  The NLSY79 and NLSY97 input files (public-use extracts, custom sampling weights, and the PCE price deflator) are hosted at:

  [https://www.dropbox.com/scl/fo/oj62jugo5lq7y7q16m9he/APjJF6_kYE8UQr44zT4vlXM?rlkey=mq003fcktit76500l69jup9rk&st=f870ynvr&dl=0](https://www.dropbox.com/scl/fo/oj62jugo5lq7y7q16m9he/APjJF6_kYE8UQr44zT4vlXM?rlkey=mq003fcktit76500l69jup9rk&st=f870ynvr&dl=0)

  Download the folder and place its contents under `Data/NLSY/` before running the NLSY pipeline.
  The underlying extracts originate from the [NLS Investigator](https://www.nlsinfo.org/investigator/) (public-use data); custom weights were produced by the NLS custom weighting service.

* **Intermediate main-results estimates (for exact reproducibility)**
  Although the full pipeline is seed-controlled, small numerical differences across computing environments (e.g., BLAS/OS floating-point variation on HPC systems) can prevent byte-exact reproduction of downstream tables and figures. To guarantee exact reproducibility, the authors also provide the intermediate estimation results used in the paper — `main_rst_baseline.rds`, `main_rst_boot.rds`, `main_rst_bc.rds`, and the analysis sample `gss_fc_occ10_5class.rds` — at:

  [https://www.dropbox.com/scl/fo/1pwccchpxg1pdtaf988dr/AOtzjzcbigbAVAYUb_xevEs?rlkey=banzdnevqw63t8bbfen4ne1hf&st=y7xaedxo&dl=0](https://www.dropbox.com/scl/fo/1pwccchpxg1pdtaf988dr/AOtzjzcbigbAVAYUb_xevEs?rlkey=banzdnevqw63t8bbfen4ne1hf&st=y7xaedxo&dl=0)

  Download the folder and place its contents under `Data/main_results/`. These files are consumed by the main-results figure scripts, the sparsity perturbation analysis (Appendix K, whose benchmark checks reference exactly this vintage), and the NLSY pipeline's GSS comparison overlays (Appendix I). Users who prefer to regenerate them from scratch can instead run `master_main_results.R` first.

---

## Reproducing Baseline Results Locally

The files included in this repository are sufficient to reproduce:

* All **baseline estimation results** for benchmark analysis and robustness check.

---

## Full Replication (Recommended: HPC)

Running the full set of robustness analyses can be computationally intensive.
We therefore recommend using **high-performance computing (HPC)** resources.

A sample SLURM-compatible submission script is provided:

* `sample.sh`

This file includes suggested settings for memory, CPU cores, and runtime.

To replicate the full set of results, run the master scripts:

* `master_main_results.R` — main text figures and Appendix D
* `master_robust_ABCE.R` — Appendices A (age), B (class typologies), C (weights), and E (alternative measures)
* `master_robust_FGHJ.R` — Appendices F (gender), G (race), H (smoothing), and J (extended cohorts)
* `master_mother_role_estimation.R` — maternal-role analyses in Appendix F (mother–son/daughter occupational pairs and the educational-attainment extension)
* `master_nlsy_pipeline.R` — Appendix I (NLSY income and occupational mobility; requires the NLSY data from the Dropbox link above)
* `master_perturbation.R` — Appendix K (sparsity perturbation analysis; requires the intermediate main-results estimates from the Dropbox link above, or a completed `master_main_results.R` run)

Typical runtime is approximately **2–4 hours per major estimation block**, depending on hardware.


---

## Acknowledgements

* The bulk of the computational implementation in this repository — including the maternal-role and educational-mobility robustness pipelines, the NLSY replication package, and the sparsity perturbation analysis — was developed with the assistance of **Claude Fable 5** (Anthropic).
* Some plotting functions were developed with reference to outputs generated using **ChatGPT-4o**.
* Functions related to **Hellinger’s dependence measure** are adapted from code provided by **Prof. Thomas Coleman**, whom the authors thank for introducing valuable perspectives on dependence measures.
* The authors thank **Luke Kushner** for his help with the NLSY occupational crosswalk.

---

##  Contact

For questions regarding the code or replication:

**Weiqi Wang**
📧 [weiqiw@uchicago.edu](mailto:weiqiw@uchicago.edu)
