## Transfer of regulatory knowledge from human to mouse for functional genomic analysis
### Abstract
Transcriptome profiling followed by differential gene expression analysis often leads to unclear lists of genes which are hard to analyse and interpret. Functional genomic tools are powerful approaches for downstream analysis, as they summarize the large and noisy gene expression space in a smaller number of biological meaningful features. In particular, methods that estimate the activity of processes by mapping transcripts level to process members are popular. However, footprints of either a pathway or transcription factor (TF) on gene expression show superior performance over mapping-based gene sets. These footprints are largely developed for human and their usability in the broadly-used model organism Mus musculus is uncertain. Evolutionary conservation of the gene regulatory system suggests that footprints of human pathways and TFs can functionally characterize mice data. In this paper we analyze this hypothesis. We perform a comprehensive benchmark study exploiting two state-of-the-art footprint methods, DoRothEA and an extended version of PROGENy. These methods infer TF and pathway activity, respectively. Our results show that both can recover mouse perturbations, confirming our hypothesis that footprints are conserved between mice and humans. Subsequently, we illustrate the usability of PROGENy and DoRothEA by recovering pathway/TF-disease associations from newly generated disease sets. Additionally, we provide pathway and TF activity scores for a large collection of human and mouse perturbation and disease experiments (2,374). We believe that this resource, available for interactive exploration and download (https://saezlab.shinyapps.io/footprint_scores/), can have broad applications including the study of diseases and therapeutics.

### Data availability
#### DoRothEA regulons
The human and mouse specific DoRothEA regulons are stored in [`data/dorothea_benchmark/regulons/`](https://github.com/saezlab/ConservedFootprints/tree/master/data/dorothea_benchmark/regulons). They are present as datatables. To convert them to the required `viper` format please use the function [`df2regulon()`](https://github.com/saezlab/ConservedFootprints/blob/3cee21853a90f78bd13b6eebcb5538fc4c129cab/src/dorothea_analysis.R#L105-L120). It exist also the complementary function allowing the transformation from viper format to datatable (see [`regulon2df()`](https://github.com/saezlab/ConservedFootprints/blob/3cee21853a90f78bd13b6eebcb5538fc4c129cab/src/dorothea_analysis.R#L122-L135)).

#### PROGENy matrices
The extended human and mouse specific PROGENy matrices are stored in [`data/progeny_benchmark/models`](https://github.com/saezlab/ConservedFootprints/tree/master/data/progeny_benchmark/models). The activity of the following pathways can be inferred from mouse and human gene expression data:

* Androgen (new)
* EGFR
* Estrogen (new)
* Hypoxia 
* JAK-STAT 
* MAPK
* NFkB 
* PI3K
* TGFb
* TNFa
* Trail
* p53
* VEGF
* WNT (new)

#### Disease sets
The newly generated disease sets are available in [`paper/material/disease_sets.csv`](https://github.com/saezlab/ConservedFootprints/blob/master/paper/material/disease_sets.csv). 

#### Resource of pathway and TF activites of perturbation and disease experiments
Pathway and TF activities of perturbation and disease experiments can be browsed in a user friendly web application available at [https://saezlab.shinyapps.io/footprint_scores](https://saezlab.shinyapps.io/footprint_scores).

#### Benchmark data
The benchmark data are not deposited in GitHub due to their file size. We are happy to share the data upon a reasonable request (christian.holland@bioquant.uni-heidelberg.de).

### Analyses
#### PROGENy benchmark (Fig. 2)
Analyses and plotting scripts available in [`analyses/progeny_benchmark`](https://github.com/saezlab/ConservedFootprints/tree/master/analyses/progeny_benchmark).

#### DoRothEA benchmark (Fig. 3)
Analyses and plotting scripts available in [`analyses/dorothea_benchmark`](https://github.com/saezlab/ConservedFootprints/tree/master/analyses/dorothea_benchmark).

#### Disease enrichment (Fig. 4)
Analyses and plotting scripts available in [`analyses/disease_enrichment`](https://github.com/saezlab/ConservedFootprints/tree/master/analyses/disease_enrichment).
