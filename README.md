# Chen_Loubiere_2025

System requirements:
  - Custom scripts generated for this study were written in R (version 4.4.1) using the R studio IDE (https://www.R-project.org/).
  - No special hardware should be required.

Installation guide:
  - R and RStudio can be downloaded at https://posit.co/download/rstudio-desktop/. Installation time is approximately 20min.
  - Custom functions were wrapped into a R package that can be installed using devtools: `devtools::install_github("vloubiere/vlite@85ce0c1")`
  
Instructions for use:
  - The "main.R" file lists all the scripts that are needed for the pre-processing and annotation of ATAC-seq peaks, data preparation before deep learning, model performance evaluation, and for all downstream analyses and plots presented in the manuscript. Script names should be self-explanatory, with additional comments linking them to the corresponding figure panel when relevant.

1/ The "Functions" section contains a simple custom function used for data augmentation.

2/ The "PREPARE DATA FOR TRAINING" section contains the scripts to pre-process VISTA sequences and annotate ATAC-seq peaks, until the CV scheme where training/validation/test folds are defined.

3/ The "DATA AUGMENTATION" section is where the data augmentation was made

4/ The "EXTRA TEST SETS FOR CV" was used to subset held-out sequences that were not used for training (tiling of the whole chromosome 18, test set split on all chromosomes...)

5/ The "SELECTED SEQUENCES FOR VALIDATIONS" contains the post-processing and selection steps for Ledidi-designed synthetic enhancers

6/ The "PAPER FIGURES" section contains all the scripts required to generate the plots from the manuscript.

Running the whole set of analyses should be doable within 2 days. Further details regarding the analyses can be found in the Mehods section on the paper.

For any reasonable further request, please contact alexander.stark@imp.ac.at, shenzhi.chen@imp.ac.at or vincent.loubiere@imp.ac.at.
