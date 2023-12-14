# RNA-seq analysis and manuscript figure generation

## Overview
The R notebook contains all analysis code in sequential order required to reporduce most findings and figures from the manuscript. The survival analysis and state-transition modeling was performed using MATLAB (see [MATLAB](https://github.com/cohmathonc/CML_mRNA_state-transition/MATLAB/)).

## Setup
1. **Get data**
  The only data required to run the R notebook is the `SummarizedExperiment` object saved as an R `rds` file located at this *temporary location*: [data](https://drive.google.com/file/d/1q0UCuB5h-9lyRIZBSrIQ-joanfVms0ZK/view?usp=sharing).

2. **Run notebook**
   A single [notebook](https://github.com/cohmathonc/CML_mRNA_state-transition/R/analysis+figures_reproducibility.R) includes all analysis performed in the [manuscript](https://www.biorxiv.org/content/10.1101/2023.10.11.561908v2). The required R libraries are listed in the "Setup" section. When run sequentially, the notebook will perform all analysis, generate intermediate objects and figures, and generate all figures included in the manuscript. Successfully running the [notebook](https://github.com/cohmathonc/CML_mRNA_state-transition/R/analysis+figures_reproducibility.R) will product the following output directories:
   
    - **`manuscript_files`**:
    
      All figures and tables included in the [manuscript](https://www.biorxiv.org/content/10.1101/2023.10.11.561908v2).
    
    - **`DEG_output`**
    
      Tables (.tsv) and enrichment analysis output for all comparisons made using `DESeq2`.
    
    - **`CML_contribution`**

      Tables for each DEG comparison that summarize the Eigengene-based CML contribution of each gene.

    - **`plots`**
  
      Miscellaneous plots created during analysis.

## Contact

- [David Frankhouser](mailto:dfrankhouse@coh.org) @dfrankhouser
    
     
