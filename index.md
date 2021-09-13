## Multi-omics (and data integration)

### Websites of particular interest

 tutorials trajectory analysis:

https://rnabioco.github.io/cellar/5_trajectories.html

https://github.com/quadbiolab/scRNAseq_analysis_vignette/blob/master/Tutorial.pdf



https://www.bioconductor.org/packages/release/workflows/vignettes/TCGAWorkflow/inst/doc/TCGAWorkflow.html#Transcriptomic_analysis

https://seqqc.wordpress.com/2015/02/16/should-you-transform-rna-seq-data-log-vst-voom/



http://bioconductor.org/packages/release/bioc/vignettes/recount/inst/doc/recount-quickstart.html



https://github.com/agitter/single-cell-pseudotime



https://bioconnector.github.io/workshops/ (e.g.https://bioconnector.github.io/workshops/r-survival.html)



https://yulab-smu.top/biomedical-knowledge-mining-book/



https://github.com/kevinblighe?tab=repositories; among which https://github.com/kevinblighe/EnhancedVolcano



https://github.com/kevinblighe/awesome-single-cell

List of software packages (and the people developing these methods) for single-cell data analysis, including RNA-seq, ATAC-seq, etc. Contributions welcome...

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



The clusters were tested for functional enrichment with a hypergeometric test and categories from Gene Ontology (www.geneontology.org/), KEGG (http://www.genome.jp/kegg/), and Reactome databases (http://www.reactome.org/), using the topGO [99]with a ‘weight01’ algorithm, Category and ReactomePA libraries, respectively

 https://github.com/bioFAM/MOFA. 

Sankey diagrams in R https://christophergandrud.github.io/networkD3/#sankey



recount3, a resource consisting of over 750,000 publicly available human and mouse RNA
sequencing (RNA-seq) samples uniformly processed by our new Monorail analysis pipeline. To facilitate
access to the data, we provide the recount3 and snapcount R/Bioconductor packages as well as
complementary web resources. Using these tools, data can be downloaded as study-level summaries or
queried for specific exon-exon junctions, genes, samples, or other features. Monorail can be used to
process local and/or private data, allowing results to be directly compared to any study in recount3. Taken
together, our tools help biologists maximize the utility of publicly available RNA-seq data, especially to
improve their understanding of newly collected data. recount3 is available from http://rna.recount.
bio

ReactomePA https://yulab-smu.top/biomedical-knowledge-mining-book



Data science with R https://garrettgman.github.io/

scenic, a computational method for simultaneous gene regulatory network reconstruction and cell-state identification from single-cell rna-seq data (http://scenic. aertslab.org).



Multi-omics Master- Regulator Analysis (MOMA). MOMA integrates gene expression and genomic alterations profiles to identify MR proteins and MR modules that represent the key effectors of a tumor’s mutational landscape and are thus responsible for implementing its associ- ated cancer cell identity.

The full code base for PathTurbEr is available in here: https://github.co m/Akmazad/PathTurbEr/ 

y Bayesian approach in discovering novel driver bio-markers in aberrant STPs given high-throughput gene expression (GE) data. 



sgnesR (Stochastic Gene Network Expression Simulator in R) is an R package that provides an interface to simulate gene expression data from a given gene network using the stochastic simulation algorithm (SSA)



e a method that compares multiple networks ofunlimited size at the level oflinks and nodes. Our novel method, CoDiNA (Co-expression Differencial Network Analysis), is imple-
mented as an R package that also includes an interactive tool for network visualization

### Python packages

 https://github.com/bioFAM/MOFA. 

https:// github.com/BeautyOfWeb/Multiview-AutoEncoder 

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





https://osf.io/gqrz9  ([Tatlow and Piccolo 2016](https://www.nature.com/articles/srep39259)) RNA-seq data from TCGA and CCLE -> transcript level expression estimates; can be combined to gene level expression estimates using Scater (eg. Campbell and Yau 2018)



https://github.com/mskcc/RNAseqDB



CPTAC datasets

Oslo2 landscape cohort

ubiquitination data (Kim et al 2011)

Cosmic: genes causally associated wiht cancer

protein/mRNA half-life data: Schwanhäusser et al. 2011 (mice)



Cancer-related SAAVs are cataloged in the CanProvar and COSMIC databases.



published cell line perturbation experiments from the Genomics of Drug Sensitivity in Cancer (GDSC) resource (Iorio et al., 2016; Yang et al., 2013)



Clinical Proteomic Tumor Analysis Consortium (CPTAC) CPTAC contains mass spectrometry-based proteomic analysis of tumors from TCGA. The aim of CPTAC is to create a proteogenomic resource where dysregulated proteins and phosphorylation sites can be identified and potentially connected to genomic alterations.

Proteomics Identification Database (PRIDE)

PRIDE aims to be a resource for open access sharing of mass spectrometry data, not just across cancer. They currently have over 9200 datasets available, including 297 breast cancer datasets.





GENIE GENIE combines genomic and clinical data in an attempt to associate genomic alterations with phenotypic changes



 GXB GXB compiles immunological transcriptomic data



Human Proteome Organization (HUPO) The human proteome project, run by HUPO aims to identify all the proteins in the human proteome and to begin to assess their functionalities and interactions



Transciptome Alterations in Cancer Omnibus (TACCO)

TACCO is a resource for identifying differentially regulated transcripts within different cancer types and combining these with survival data to determine prognosis



COSMIC

COSMIC contains data from over 13 million tumor samples, identifying 6 million coding mutations and over 19 million non-coding mutations. This resource collates all genes implicated in cancer through somatic mutation, of which 719 are currently listed



MOSCATO trial where druggable genomic aberrations were identified and targeted in patients (Massard et al., 2017



protein complexes, specifically the CORUM database





Sweden Cancerome Analysis Network-Breast (SCAN-B)



GO semantic similarity between genes can be calculated using tools GOSemSim (Yu et al. 2010) or GOssTo (Yu et al. 2010)



functional association between genes using data mining and text mining and can be valuable and reliable resources for generating a prior matrix, e.g., STING (Von Mering et al. 2005) and AraNet (Lee et al. 2010)



s an interactive online resource with navigable
proteomics, transcriptomics, and drug sensitivity profiles at
https://lehtio-lab.se/forall/.



the first large-scale multi-focal breast cancer proteomic study of 330 tumor regions which associated cancer cell function, pathological parameters, and spatial localization of each tumor region. 

All raw data are available via ProteomeXchange with identifier PXD024190.



MSigDB [93]website (http://software. broadinstitute.org/gsea/msigdb/) the v5.1 C7 (‘immuno- logic signatures’) collection



subcellular localization was retrieved from the
Human Protein Atlas [97](www.proteinatlas.org).



 unsupervised pathway enrichment analysis using **Reactome30**



We used ARACNe [100] to infer edges between the
hubs and the expressed genes. 



TF genes identified with ei- ther or both of two alternative annotations: (1) the human genes with a symbol annotated with the term ‘GO: 0003700’ in the Gene Ontology Consortium database (www.geneontology.org) or (2) the Ensembl gene ID re- trieved by querying the BioMart service (http://grch37. ensembl.org/) with the Gene Ontology ID ‘GO:0003700



DNase-Seq data were obtained from ENCODE (http:// hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/ wgEncodeUwDgf, sample ‘wgEncodeUwDgfTregwb7849 5824’) with DNase hypersensitivity (DHS) peaks. Protein- binding footprints (FP) within the DHS peaks were identi- fied by scanning the DHS intervals for gaps in the signals using Wellington [101],



We considered data from three different public PPI data repositories (Human Interactome database [106], STRING v10.0 [107], and iRefIndex [108]), and integrated them into a large human PPI network with



To derive association scores from the Open Targets re- source (www.opentargets.org) the version from September 2016 was used



o draw a treemap using the treemap library



TcellSubC: An Atlas of the Subcellular Proteome of Human T Cells



implemented the pypgatk package and the pgdb workflow to create proteogenomics databases based on ENSEMBL resources

The tools allow the generation of protein sequences from novel protein-coding transcripts by performing a three-frame translation of pseudogenes, lncRNAs, and other non-canonical transcripts, such as those produced by alternative splicing events. It also includes exonic out-of-frame translation from otherwise canonical protein-coding mRNAs. Moreover, the tool enables the generation of variant protein sequences from multiple sources of genomic variants including COSMIC, cBioportal, gnomAD, and mutations detected from sequencing of patient samples

pypgatk: (https://github.com/bigbio/py-pgatk/), and pgdb: (https://github.com/nf- core/pgdb ) 

the Proteomics-Genomics DataBase (pgdb - https://nf-co.re/pgdb) workflow



(Dietrich et al, 2018

a cohort of 200 patient samples of chronic lymphocytic leukaemia, profiled for somatic mutations, RNA expression, DNA methylation and ex vivo drug responses



Mertins et al. 2016 e quantitative mass-spectrometry-based proteomic and phosphoproteomic analyses of 105 genomically annotated breast cancers, of which 77 provided high-quality data





##### **PCaDB**                   http://bioinfo.jialab-ucr.org/PCaDB/  

PCaDB is a comprehensive and interactive database for transcriptomes from prostate cancer population cohorts.            We collected 50 transcriptomics datasets with 7,231 samples from public data repositories,            including TCGA, cBioPortal, GEO, and ArrayExpress. A standard bioinformatics pipeline is used to download and            process the expression data and metadata.           PCaDB provides a user-friendly interface for the comprehensive analysis of individual            genes, prognostic signatures, and the whole transcriptomes to elucidate the molecular            heterogeneity in PCa, understand the mechanisms of tumor initiation and progression,            as well as develop and validate prognostic signatures in large independent cohorts



https://cptac-data-portal.georgetown.edu/cptac/s?id=393177





Characterizing genetic influences on DNA methylation (DNAm) provides an  opportunity to understand mechanisms underpinning gene regulation and  disease. In the present study, we describe results of DNAm quantitative  trait locus (mQTL) analyses on 32,851 participants, identifying genetic  variants associated with DNAm at 420,509 DNAm sites in blood. http://mqtldb.godmc.org.uk/ 



The Human Developmental Cell Atlas (HDCA) initiative, which is part of  the Human Cell Atlas, aims to create a comprehensive reference map of  cells during development



The marked improvements in massive parallel sequencing coupled with  single-cell sample preparations and data deconvolution have allowed  single-cell RNA sequencing (scRNA-seq) to become a powerful approach to  characterize the gene expression profile in single cells ([*1*](https://www.science.org/doi/10.1126/sciadv.abh2169#pill-R1), [*2*](https://www.science.org/doi/10.1126/sciadv.abh2169#pill-R2)). The objective of the international collaborative effort Human Cell Atlas ([www.humancellatlas.org](http://www.humancellatlas.org)) takes advantage of this new technology platform to study the  distinctive gene expression profiles on RNA level across diverse cell  and tissue types and connect this information with classical cellular  descriptions, such as location and morphology ([*3*](https://www.science.org/doi/10.1126/sciadv.abh2169#pill-R3)). 

The objective of the Human Protein Atlas (HPA) ([www.proteinatlas.org](http://www.proteinatlas.org)) effort is to take advantage of these bioimaging approaches to map the  expression of all human protein-coding genes across all major human  cells, tissues, and organs.



 an effort to combine the information from these two efforts to create a publicly available HPA Single Cell Type Atlas with genome-wide  expression data from scRNA-seq experiments integrated with the spatial  antibody-based bioimaging data.  a Single Cell Type Atlas has been 
launched (www.proteinatlas.org/celltype)  A single–cell type transcriptomics map of  
human tissues



http://gent2.appex.kr/

Now, GENT2 contains more than 68,000 samples and has several new useful  functions. First, GENT2 now provides gene expression across 72 different tissues compared to 57 in GENT. Second, with increasing importance of  tumor subtypes, GENT2 provides an option to study the differential  expression and its prognostic significance based on tumor subtypes.  Third, whenever available, GENT2 provides prognostic information of a  gene of interest. Fourth, GENT2 provides a meta-analysis of survival  information to provide users more reliable prognostic value of a gene of interest.



**[OncoLnc](http://www.oncolnc.org/)**: [oncolnc.org](http://www.oncolnc.org/)

- Focus on survival analysis and RNA-seq data.
- Simple query interface across all cancers for any mRNA, miRNA, or lncRNA gene (try SERPINA1)
- Precomputed Cox PH regression for every gene, for every cancer
- Kaplan-Meier plots produced on demand

[TANRIC](http://ibl.mdanderson.org/tanric/_design/basic/index.html): focus on noncoding RNA

[MEXPRESS](http://mexpress.be/): focus on methylation and gene expression





The Database of Interacting Proteins (DIP: http://dip.doe-mbi.ucla.edu) is a database that documents experimentally determined protein–protein interactions. https://doi.org/10.1093/nar/30.1.303



 MR4Cancer: a web server prioritizing master regulators for cancer  https://doi.org/10.1093/bioinformatics/bty658 http://cis.hku.hk/MR4Cancer



Tovar et al. 2015

An important step for this algorithm is the selection of the Transcription Factors, since they will determine the rest of the cal- culation. A proper annotation of transcription factors is crucial for an accurate description of the process under investigation. Here, we used the HGU133A annotation file, in which we found 1142 TFs (Supplementary Material 1). This list was compared with other three lists. Those lists are available in Shimoni and Alvarez (2013), Vaquerizas et al. (2009) and http://www.bioguo.org/AnimalTFDB/ Animal Transcription Factor DataBase, respectively. We want to stress that all four lists show consistency among them



d fusion events (FUSs) reported by Pipeline for RNA-Sequencing Data Analysis (PRADA) (Torres- Garcı´a et al., 2014)(



VIPER has been extensively validated as an accurate method-
ology to measure a protein’s activity, on the basis of the enrichment of its tissue-specific activated and repressed tran- scriptional targets (regulon) in over and under-expressed genes (Alvarez et al., 2016)—i.e., akin to a highly multiplexed gene- reporter assay.



human protein-protein interaction networks as of the HINT database2



essential genes from the online gene essentiality data- base (OGEE)24 and the Database of Essential genes (DEG)25



1,471 manually curated sequence-specific
DNA-binding transcription factors28,29 and 501 genes from the Kinome NetworkX database27 that collects kinase information from the literature and other databases



17,511 phosphorylated proteins, 6,928 acetylated
proteins and 5,418 methylated proteins from the PhosphoSitePlus database44



genes that were annotated with a signaling function without recep-
tor domain function from Gene Ontology (GO)45 as well as 5,701 genes that carried a trans-membrane protein domain

95,722 links between 209 human transcription factor and 8,910 human genes from the TRANSFAC42 database as provided by mSigDB



a set of protein abundances30 in human cell lines. Beck, M. et al. The quantitative proteome of a human cell line. Mol Syst Biol 7, 549 (2011)

Interactome mapping :There are several widely-used databases – BioGrid,16
IntAct,17 HPRD,18 iRefWeb,19 DIP,20
MINT,21 MIPS22 and VisAnt23 – that curate both categories of interactions for humans and other model organisms.

HINT also distinguishes between interactions curated from small-scale studies and those obtained from high-throughput experiments.



WikiPathways6, Ingenuity Pathway Analysis (qiagen.com), KEGG7, and Reactome8





marina to infer master regulator



the Human Protein Reference Database (HPRD) [28], as one of the most com- prehensive protein interaction databases, only covers less than half of human protein-coding genes. Therefore, (Zhang et al. 2011)



known associations between disease phenotypes and genes extracted from the Online Mendelian Inheritance in Man (OMIM) database

he Human Protein Reference Database (HPRD) contains human protein-protein interactions that are manually extracted from the literature by expert biolo- gists [28].

Biological General Repository for Interaction Data- sets (BioGRID) contains protein and genetic interactions of major model organism specie

the Biomolecular Interaction Network Database (BIND) contains both high-throughput and manually curated interactions between biological molecules [30]. 

the IntAct molecular interaction database (IntAct) contains protein-protein interaction derived from literature [

he Molecular INTeraction database (MINT) contains information about physical interactions between pro- teins [32]. 



cataloging of protein–protein interactions across
species and conditions into databases such as STRING6

The protein interaction data were derived from the BioGRID database. 

The publicly available data sets, such as Biological General Repository for Interaction Datasets (BioGRID, https://thebiogrid.org/), Saccharomyces Gen- ome Database (SGD, http://www.yeastgenome.org/), Hu- man Protein Reference Database (HPRD, http:// www.hprd.org/), Search Tool for the Retrieval of Interacting Genes/Proteins (STRING: http://www.string- db.org/)
