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



This report proposes a framework for integrating known genetic pathways into a differential network analysis
of two populations. The framework allows any association measure to be used, and a general measure for differen- tial connectivity is considered. Statistical significance is evaluated through a permutation testing procedure. The methodology is implemented in R and is available on GitHub at https://github.com/tgrimes/dnapath



 six additional mutual information methods in the MINET R package (ARACNE, CLR, MIM, MINET, MRNET, MRNETB) [[43](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007241#pcbi.1007241.ref043)



The PAM505 classifier, which is also deployed in Genefu Bioconductor package20, makes calls based on the 50 gene centroid correlation distance to LA, LB, Basal, Her2 and normal-like centroids



Proteomic data were accessed and downloaded using the R package “TCGA-Assembler 2” (27, 28) from the CPTAC



All computations were carried out on the R platform. Package “fastICA” which implements the iterative FastICA algorithm (15) was used to extract non-Gaussian independent components with logcosh contrast function. Components were subsequently assigned to clusters using the “cluster” package. Clusters were visualized with 2d t-SNE using the R package “tsne.” The number of clusters was determined as equal to number of components extracted at each run of ICA. When the number of samples is small comparing to the number of features, which is usually the case for biological data, it is convenient to retrieve as many as independent signal sources as possible, and the number of compo- nents extracted is equal to sample size



The package “pcaMethods” was used to calculate principal com-
ponents for comparison with independent components. 



Gene Set Enrichment Analysis (GSEA) implemented in the WebGestaltR R-package (Liao et al., 2019) was used to infer signatures of approved drugs (D1, 1,202 gene sets) and kinase inhibitors (D2, 1,220 gene sets) available in the drug signatures database (DSigDB; Yoo et al., 2015; http://dsigdb.tanlab.org/DSigDBv1.0/). Based



Stemness scores were calculated as previously described (Malta et al., 2018). Firstly, we used MoonlightR (Colaprico et al., 2020)to query, download, and preprocess the pluripotent stem cell samples (ESC and iPSC) from the Progenitor Cell Biology Consortium (PCBC) dataset



[PaxtoolsR: pathway analysis in R using Pathway Commons](https://scholar.google.ch/citations?view_op=view_citation&hl=fr&user=njFSsvYAAAAJ&cstart=20&pagesize=80&sortby=pubdate&citation_for_view=njFSsvYAAAAJ:4TOpqqG69KYC)



### Python packages

 https://github.com/bioFAM/MOFA. 

https:// github.com/BeautyOfWeb/Multiview-AutoEncoder 



https://github.com/raphael-group/hotnet2

HotNet2 identifies subnetworks of a protein-protein interaction network with more mutations ("heat") than expected.



### Other tools

clustering and heatmap visualization 

(https://software.broadinstitute.org/GENE-E/index.html). Pearson

https://software.broadinstitute.org/morpheus/



GO enrichment by Gorilla (http://cbl- gorilla.cs.technion.ac.il) 



e screened for enrichment in MSigDB and chromosome positions using the hypergeometric testing implemented in the GeneListEnrich shiny tool (https://github.com/aleferna/BCLandscape_Shiny/tree/master/GeneListCompare)



Gene ontology (GO) analyses were performed using ToppGene (https//topgene.cchmc.org)(



CIBERSORT RNA-based estimates of overall I-TME provided by CI- BERSORT absolute scores (immune cell infiltration from RNA data)



A website for interactive visualization of the multi-omics dataset is available at: http://prot-shiny-vm.broadinstitute.org:3838/ CPTAC-BRCA2020. The



gephi for graph vizualization



The RNA-based tumor microenvironment inference tool ESTIMATE (Yoshihara et al., 2013) was used to derive the overall immune score and stromal score for each sample. 



Weused an established methylation-based deconvolution method, EDec (Onuchic et al., 2016) to dissect the composition of different cell types within the whole bulk tumor. I

We integrated somatic mutation, CNV, DNA methylation, RNA, protein, phosphorylation (phospho) and acetylation (acetyl) levels via iProFun (Song et al., 2019) to investigate the functional impacts of DNA alterations in GBM.



functional characterization of clustering results by single sample Gene Set Enrichment Analysis (ssGSEA)

n https://github.com/broadinstitute/ ssGSEA2.0



PANOPLY: a cloud-based platform for automated and reproducible proteogenomic data analysis. -> encapsulate the complex data processing required for proteogenomics  and provide a simple interface to deploy a range of algorithms developed for data analysis, we have developed PANOPLY—a cloud-based Platform for Automated aNd reprOducible Proteogenomic data anaLYsis



SignatureAnalyzer exploited the Bayesian variant of the NMF algorithm and enabled an inference for the optimal number of signatures from the data itself at a balance between the data fidelity (likelihood) and the model complexity (reg-



Lolliplots for CDKN2A mutations in the CPTAC head and neck squamous cell carcinoma and this lung squamous cell carcinoma cohort were generated using the ProteinPaint web application to visualize mutations (Zhou et al., 2016). Mutations



Pathway enrichment results were imported into Cytoscape (Shannon et al., 2003) using the Enrichment Map app (Merico et al., 2010) for network analysis of pathways.

BreakDancer (Chen et al., 2009) and Meerkat (Yang et al., 2013) algorithms were used to detect structural variations.

PARADIGM Integrated Pathway Analysis Integrated Pathway Levels (IPLs) mRNA expression, SCNA, and pathway interaction data for 80 UM samples were integrated using the PARADIGM software (Sedge- wick et al., 2013). Briefly, this procedure infers integrated pathway levels (IPLs) for genes, complexes, and processes, using pathway interactions, and genomic and functional genomic data from each patient sample. Normalized



hMARINa ssMARINa MARINa was run via the VIPER R package (http://www.bioconductor.org/packages/release/bioc/html/viper.html)(Alvarez et al., 2016); and hMARINa was per- formed by extending the functionality of the package



Visualization of cellular processes and pathways http://newteditor.org



We have built a microarray data analysis tool, named PATIKAmad, which can be used to associate microarray data with the pathway models in mechanistic detail, and provides facilities for visualization, clustering, querying, and navigation ofbiological graphs related with loaded microarray experiments. PATIKAmad is freely available to noncommercial users as a new module ofPATIKAweb at http://web.patika.org.



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



PPIs from HINT database



pathway enrichment analysis, using Kyoto Encyclopedia of Genes and Genomes (KEGG), Biocarta, and Reactome pathways from the Molecular Signature Database (MSigDB) (

This was further complemented by incorporating Kyoto Encyclopedia of Gene and Genomes (KEGG) ([Ogata et al., 1999](https://www.sciencedirect.com/science/article/pii/S1476927119303871?casa_token=TTrIGf-vCaAAAAAA:w4xLI5CZzXwhQlJZkG0oL7wdOFfgYDNGK7CwLXnGaxygNkR6XNnWmiDSJkkiwx33yBWzwqcbkuVd#bib0180)) PFAM protein domains ([Sonnhammer et al., 1998](https://www.sciencedirect.com/science/article/pii/S1476927119303871?casa_token=TTrIGf-vCaAAAAAA:w4xLI5CZzXwhQlJZkG0oL7wdOFfgYDNGK7CwLXnGaxygNkR6XNnWmiDSJkkiwx33yBWzwqcbkuVd#bib0240)) and INTERPRO protein domains ([Mitchell et al., 2014](https://www.sciencedirect.com/science/article/pii/S1476927119303871?casa_token=TTrIGf-vCaAAAAAA:w4xLI5CZzXwhQlJZkG0oL7wdOFfgYDNGK7CwLXnGaxygNkR6XNnWmiDSJkkiwx33yBWzwqcbkuVd#bib0160)) information of the selected genes.

list of transcription factors in the TFCheckpoint curated database (Tripathi et al., 2013).





 GRAND (https://grand.networkmedicine.org) as a database for computationally-inferred, context-specific gene  regulatory network models that can be compared between biological  states, or used to predict which drugs produce changes in regulatory  network structure. The database includes 12 468 genome-scale networks  covering 36 human tissues, 28 cancers, 1378 unperturbed cell lines, as  well as 173 013 TF and gene targeting scores for 2858 small  molecule-induced cell line perturbation paired with phenotypic  information. GRAND allows the networks to be queried using phenotypic  information and visualized using a variety of interactive tools.

(http://www.proteinatlas.org),



we examined genomic-based patterns of oncogenic pathway activity, the tumor microenvironment and other important features in human breast tumors using a panel of 52 previously published gene expression signatures (Supplementary Table - An integrated genomics approach identifies drivers of proliferation in luminal-subtype human breast cancer)



an RNAi proliferation screen in which a genome-wide shRNA library (~16,000 genes) had been used to identify essential genes (An integrated genomics approach identifies drivers of proliferation in luminal-subtype human breast cancer) Marcotte et al. 2012



Gene Active Ranking Profile (GARP)-normalized data were obtained from the COLT database



Significant pathway annotations (FDR < 0.05) from the Panther over-representa- tion test (database: GOBP complete, http://www.pantherdb. org/) were used to annotate each cluster



TCGA aggregates an extensive col- lection of omics and clinical datasets from large cohorts of patients for more than 30 types of cancers (24). It also ar- chives histopathology images for solid tumor samples from which omics data were sampled. Currently, more than 24,000 histopathology images are available and can be visualized at the Cancer Digital Slide Archive (CDSA, http://cancer. digitalslidearchive.net/). In addition, The NCI Clinical Proteomic Tumor Analysis Consortium (CPTAC) (https://proteomics. cancer.gov/programs/cptac) program also provides high- throughput proteomic data for some of the TCGA tumor spec- imens such as breast cancer, ovarian cancer, and colorectal cancer based on mass-spectrometry technology. These



To facilitate use and dissemination of the data, we have developed a web resource (https://zucchini.gs. washington.edu/BreastCancerProteome/) in which protein abundances can be queried and correlated to genomic and drug sensitivity data, as presented below. 



m the Genomics of Drug Sensitivity in Cancer (CRx) resource (Yang et al., 2013)



e https://lincs.hms.harvard.edu/db/datasets/20343

We performed quantitative proteomics on 61 human-derived breast cancer cell lines to a depth of ~13,000 proteins. The

. All datasets are freely available as public resources on the LINCS portal.

7197 proteins measured in all cell lines



single cell data breast cancer

scRNA-seq data are available for in-browser exploration and download through the Broad Institute Single Cell portal at https://singlecell. broadinstitute.org/single_cell/study/SCP1039. Processed scRNA-seq data from this study are also available through the Gene Expression Omnibus under accession number GSE176078.



known interaction pairs from CORUM database





essential survival gene datasets from The Cancer Dependency Map, the latter of which catalogs genes driving cancer progression



BioMuta and BioXpress: mutation and expression knowledgebases for cancer biomarker discovery



druggability information from DGIdb (Cotto
et al., 2018) and DEPO (Sun et al., 2018), we



Weaggregated a set of interacting proteins (e.g. kinase/phosphatase-substrate or complex partners) from OmniPath (downloaded on 2018-03-29) (T€urei et al., 2016), DEPOD(downloadedon2018-03-29) (Duanet al., 2015),CORUM(downloadedon2018-06-29) (Ruepp et al., 2010), Signor2 (downloaded on 2018-10-29) (Perfetto et al., 2016), and Reactome (downloaded on 2018-11-01) 



The Cancer Imaging Archive (TCIA) database.



a pathway enrichment analysis of Hallmark, KEGG, and Reactome. The



L1000 and P100 drug connectivity analysis

Level 4 P100 data were downloaded from the LINCS Data Portal (Stathias et al., 2019) and were used to calculate drug connec-
tivities on the phosphoprotein level as previously described (



The differentially expressed genes between gene-altered and WT samples were filtered for the 978 genes measured in the L1000
assay and then were processed using the CLUE (Subramanian et al., 2017) (summary connectivity score) and iLINCS (Pilarczyk et al., 2019) connectivity algorithms. The resulting drug connectivities were aggregated to the compound level using the summary connec- tivity score in CLUE and the Connected Perturbations Z-score in iLINCS. Target annotations for the ranked compounds were extracted from CLUE and iLINCS and combined in a single list. Level



. The PreCancer Atlas (PCA) ofthe NCI envisages a histological and multi-omic mapping strategy in time and space to provide detailedmolecular, cellular, and structural characterization of premalignant lesions and how they evolve to invasivecancers



Genes and Genomics (KEGG) (Kanehisa & Goto, 2000), Reactome (Joshi-Tope et al., 2005), WikiPathways (Slenter et al., 2018), and the SIGnaling Network Open Resource (SIGNOR) (Perfetto et al., 2016) have

Resources such as the Human Protein Reference Database (HPRD) (Peri et al., 2004; Keshava Prasad et al., 2009), the Biological General Repository for Interaction Datasets (BioGRID) (Stark et al., 2006), and the Search Tool for Retrieval of Interacting Genes/Proteins (STRING) (Szklarczyk et al., 2019) 



Popular phosphorylation focused databases such as PhosphoSitePlus (Hornbeck
et al., 2015), PHOSIDA (Gnad et al., 2007), Phospho.ELM (Diella et al., 2004) and qPhos (Yu et al., 2019) host many phosphosites and act as repositories for both low and high-throughput data.

An extensive review of phosphoproteomics resources can be found here (Savage & Zhang, 2020). 

Moving away from databases acting as phosphosite repositories, many databases holding signalling information in the form of a K-S networks exist such as RegPhos (Huang et al., 2014), PhosphoNet (Safaei et al., 2011) and Phosphonetworks (Hu et al., 2014). Though they are limited to K-S interactions, they have been proven useful in providing phosphoproteomics mechanistic insight either on their own or while integrated into other databases (Rohrs et al., 2018; McGuire et al., 2017; Tong et al., 2019).



study comparing PPI databases found 375 resources (Bajpai et al., 2020). 



Library of Integrated Network-Based Cellular Signatures (LINCS) L1000 perturbation-response signa- tures. The scores were computed using the sig_queryl1k_tool pipeline (https://hub.docker.com/u/cmap) and the LINCS L1000 Level 5 compound (trt_cp) signatures from CLUE (https://clue.io, ‘‘Expanded CMap LINCS Resource 2020 Release’’). The



Cancer-associated genes (CAG) were compiled using the Census website https://cancer.sanger.ac.uk/census. The list of genes is provided in Table S1. -> Satpathy et al. 2021



We used iProFun, an integrative analysis tool to identify multi-omic molecular quantitative traits (QT) perturbed by DNA-level varia- tions. 



Immunohistochemistry-based antibody-specific staining scores in lung tumors were obtained from the Human Protein Atlas (HPA, https://www.proteinatlas.org), in which tumor-specific staining is reported in four levels, i.e., high, medium, low, and not detected. 





We queried the PTM signatures database (PTMsigDB) v1.9.0 downloaded from https://github.com/broadinstitute/ssGSEA2.0/ tree/master/db/ptmsigdb using the flanking amino acid sequence (+/? 7 aa) as primary identifier. We used the implementation of PTM-SEA available on GitHub (https://github.com/broadinstitute/ssGSEA2.0) using the command interface R-script (ssgsea- cli.R). The



A list of known human E3 ubiquitin and ubiquitin-like ligases and DUBs was compiled from (Medvar et al., 2016; Nijman et al., 2005)



Connectivity Map (CMAP)

Dependency map (DepMap)

CT Antigen database

DepMap genetic dependency dataset (CRISPR Avana Public 20Q3) that contained 18119 genes and 789 cell lines (https://depmap.org/portal/download/ file: Achilles_gene_effect.csv). O

OncoKB database for oncoprotein and tumor suppressor classification (excluded proteins that have both annotations), and used PhosphoSitePlus and Signor databases for the activating/inhibiting classification of phosphorylation sites. 

NeoFlow (https://github.com/bzhanglab/neoflow) for neoantigen prediction (Wen et al., 2020). Specifically, Optitype (Szolek et al., 2014) was used to find human leukocyte antigens (HLA) for each sample based on WES data. Then we used netMHCpan (Jurtz et al., 2017) to predict HLA peptide binding affinity for somatic mutation–derived variant peptides with a length between 8-11 amino acids. The



We used Customprodbj (Wen et al., 2020)(https://github.com/bzhanglab/customprodbj) for customized database construction.



Remaining variant peptides were further filtered using PepQuery (http://www.pepquery.org)(Wen et al., 2019) with the p value cutoff of 0.01. Competitive filtering based on unrestricted posttranslational modification searching was enabled in PepQuery validation. The spectra of variant peptides were annotated using PDV (http://pdv.zhang-lab.org)(Li et al., 2019b



Cancer/testis (CT) antigens were downloaded from the CTdatabase (Almeida et al., 2009). C



PROGENy (Schubert et al., 2018) was used to generate activity scores for EGFR based on RNA expression data.

Pathways were obtained in BioPax Level 3 format, and included the NCIPID and BioCarta databases from http://pid.nci.nih.gov
and the Reactome database from http://reactome.org.



Creating a Curated Transcription Factor (TF) Regulome A compendium of TFs and their targets (TF regulons) were created by combining information from four databases:
(i) SuperPathway (Sedgewick et al., 2013): This is the same interaction network used in the PARADIGM analysis (above). Only links that correspond to regulation at the transcriptional level were retained for MARINa and hMARINa use.
(ii) Literome (Poon et al., 2014): The network was filtered to include only transcription links in which the regulator is a known TF. (iii) Multinet (Khurana et al., 2013): The network was reduced to links that correspond to regulation on transcriptional level. (iv) ChEA (Lachmann et al., 2010): Data from the Gene Expression Atlas (Petryszak et al., 2014) was used to filter the inferred links in the ChEA database. Specifically, the context likelihood of relatedness (CLR) method (Faith et al., 2007) was used to compute a measure of association between every pair of genes. The top 10% of gene pairs based on the CLR score were intersected with the ChEA network and the overlapping pairs were added to the final combined network.
The combined network includes 72,915 transcriptional regulatory links between 6,735 regulators and their targets. Only regulators
with at least 15 targets were considered in the final analysis, which resulted in a final network consisting of 419 TFs with 58,363 total targets (covering a set of 12,754 unique targets). 



Creating a Curated Kinase Regulome Proteins identified as kinases in Manning (Manning et al., 2002) or Uniprot (UniProt Consortium, 2014) were aggregated into a list of 546 kinases. Protein substrates were extracted from PhosphositePlus (Hornbeck et al., 2014) on March 7, 2015. Kinase-substrate interactions were retained if the kinase appeared in the Manning-Uniprot kinase list and the kinase was identified as a human protein in the PhosphositePlus database. The final compendium consisted of 5,388 links between 342 kinases and 2,260 unique substrates.



Institute for Systems Biology Regulome Explorer (http://explorer.cancerregulome.org),



databases that characterize pathways such as NetPath [13], KEGG [14], Reactome [11], and dozens of others.



Repositories such as WikiPathways [17] and Path- wayCommons [23] now contain thousands of pathways comprised of millions of interaction



Omic and Multidimensional Spatial (OMS) Atlas generated from four serial biopsies of a metastatic breast cancer patient during 3.5 years of therapy. This

the HTAN Data Coordinating Center (https://humantumoratlas.org/) as



aggregation of publicly available molecular interactions and biological pathway databases provided by the Pathway Commons (PC) resource.23 The aggregated data is represented in the standard Biological Pathway Exchange (BioPAX) language and provides the most complete and rich representation of the biological network models stored in PC. These complex biochemical reactions were reduced to pairwise relationships using rules to generate a Simple Interaction Format (SIF) representation of BioPAX interactions. The reduction of BioPAX interactions to the SIF allows for the representation of pairwise molecular interactions in the context of specific binary relationships. The



3D structures submitted directly by authors to the Protein DataBank (PDB)





Biofactoid (biofactoid.org), a web-based software system that
empowers authors to capture and share structured human- and machine-readable summaries of molecular-level interactions described in their publications



upstream regulator analysis from the commercial Ingenuity Pathway Analysis software.



Pathway Commons (https://www.pathwaycommons. org) is an integrated resource of publicly available information about biological pathways including bio- chemical reactions, assembly of biomolecular com- plexes, transport and catalysis events and physical interactions involving proteins, DNA, RNA, and small molecules 

www.pathguide.org



(https://ccr.cancer.gov/research/ cancer-moonshot).. The Human Tumor Atlas Network (HTAN), part of the National Cancer Institute (NCI) Cancer Moonshot Initiative, will establish a clinical, experimental, computational, and organizational framework to generate infor- mative and accessible three-dimensional atlases of cancer transitions for a diverse set of tumor types. This effort complements both ongoing efforts to map healthy organs and previous large- scale cancer genomics approaches focused on bulk sequencing at a single point in time. Gener- ating single-cell, multiparametric, longitudinal atlases and integrating them with clinical outcomes should help identify novel predictive biomarkers and features as well as therapeutically relevant cell types, cell states, and cellular interactions across transitions





IntAct (Orchard et al., 2014), SIGNOR (Licata et al., 2019; Perfetto et al., 2016) and Reactome (Fabregat et al., 2018)], and



