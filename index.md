## Multi-omics (and data integration)

### Websites of particular interest

 tutorials trajectory analysis:

https://rnabioco.github.io/cellar/5_trajectories.html

https://github.com/quadbiolab/scRNAseq_analysis_vignette/blob/master/Tutorial.pdf



https://www.bioconductor.org/packages/release/workflows/vignettes/TCGAWorkflow/inst/doc/TCGAWorkflow.html#Transcriptomic_analysis

https://seqqc.wordpress.com/2015/02/16/should-you-transform-rna-seq-data-log-vst-voom/



http://bioconductor.org/packages/release/bioc/vignettes/recount/inst/doc/recount-quickstart.html



https://github.com/agitter/single-cell-pseudotime

### Articles of particular interest

- [Uncovering pseudotemporal trajectories with covariates from single cell and bulk expression data](http://www.doi.org/10.1038/s41467-018-04696-6), Campbell and Yau 2018 ([reading_notes](reading_notes/Uncovering_pseudotemporal_trajectories_campbell.md))

### Books of particular interest

- 

### Reading notes

- 

### Other documents

- <a href="reading_notes/glossary.md">glossary </a>



### R packages

- **PhenoPath**. Genomic trajectories (pseudotimes) with heterogeneous genetic and environmental backgrounds. https://bioconductor.org/packages/release/bioc/html/phenopath.html.
- **dyno**. Trajectory inference workflow. https://dynverse.org/dyno/index.html.

Segmented copy number data were processed using CNtools package 

We inferred the similarity between individual cells in the samples and cell lines by calculating the pairwise Pearson correlation matrix C ={cor(i,j)} between any cell i from the Drop- seq experiments with any cell line, j in the CellAtlas leveraging the R package dmatch



Processed IlluminaHuman450Kdata(GEO accession GSE72021,described in ref. 26) from221 tumor samples (171 serous, 18 endometrioid, 14 clear cell, 9 mucinous, and 9 other histologic cancer subtypes) was downloaded, along with accompanying clinical annotations, using the GEOquery package from Bioconductor (

The InfiniumPurify R package (https://cran.r-project.org/web/ packages/InfiniumPurify/index.html) was used to estimate tumor purity levels based on DNA methylation profiles (28).

EPIC array data from JHU samples were converted to 450K meth- ylation profiles using the convertArray function in the minfi pack- age,

R package Scater

 unsupervised pathway enrichment analysis using **Reactome30**

### Python packages

### Public data

- **METABRIC** (proteomics, RNA-seq, ...)
- **TCGA** (RNA-seq, methylation, ...)
- **Protein-gene Expression Nexus**: comprehensive characterization of human cancer cell lines with proteogenomic analysis http:/combio.snu.ac.kr/pen (https://doi.org/10.1016/j.csbj.2021.08.022).
- **CancerSEA**: a cancer single-cell state atlas http://biocc.hrbmu.edu.cn/CancerSEA/ (https://doi.org/10.1093/nar/gky939)
- **MCLP**: Characterization of Human Cancer Cell Lines by Reverse-phase Protein Arrays (https://doi.org/10.1016/j.ccell.2017.01.005)
- **CCLE** The Cancer Cell Line Encyclopedia http://www.broadinstitute.org/ccle https://doi.org/10.1038/nature11003
- Kaplan–Meier Plotter is an online database that contains comprehensive clinical and microarray data for various cancers

Reverse Phase Protein Arrays (RPPAs) data

OncoKB annotation

to know the metastatic genes, HCMDB (Human Cancer Metastasis Database), an integrated database

predicted genes were checked for their differ- ential mRNA expression in oral cancer samples (OSCC/TSCC) in the Oncomine database (https ://www.oncom ine.org) and in the Expression Atlas (https ://www.ebi.ac.uk/gxa/). The analysis showed that POSTN, TNC and
ential mRNA expression in oral cancer samples (OSCC/TSCC) in the Oncomine database (https ://www.oncom ine.org) and in the Expression Atlas (https ://www.ebi.ac.uk/gxa/). The analysis showed that POSTN, TNC and FSCN1 are highly upregulated in OSCC with fold changes of 5.31,

Cancer RNA-Seq Nexus (CRN) database



. Fudan University Shanghai Cancer Center (FUSCC) presented comprehensive clinical, genomic, transcriptomic data of 465 primary TNBC19, which is the largest TNBC genomic project to date





To assign cell types to individual cells, we used a bulk RNA sequencing data-set from 95 cell lines collected by CellAtlas (Mabbott et al., 2013) that covered 33 major cell types in normal human tissue, including common immune, endothelial, epithelial, fibroblast and mesodermal cells. T



gene expression and GDSC2 drug sensitivity datasets from the Genomics of Drug Sensitivity in Cancer (GDSC) (Yang et al., 2013) . 



Transcription factor peaks (bed files) processed from cistrome.org



Recount2 - Rpackage

One issue that can be encountered when planning DEA of TCGA data is  the fact that some projects on the GDC portal do not contain normal  control samples for the comparison with the tumor samples. As explained  previously, it is now possible to query data from the *Recount2* platform to increase the pool of normal samples and apply the DEA pipelines of *TCGAbiolinks* (see [Fig 4A](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006701#pcbi-1006701-g004) for a workflow).

-> download GEO and GTEx



https://github.com/markrobinsonuzh/conquer

*conquer* (**con**sistent **qu**antification of **e**xternal **R**NA-seq data sets) repository, which provides access to consistently processed, analysis-ready single-cell RNA-seq data sets, together with quality  control and exploratory analysis reports

http://imlspenticton.uzh.ch:3838/conquer/