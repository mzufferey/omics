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



Data science with R https://garrettgman.github.io/



https://github.com/kevinblighe/awesome-single-cell

List of software packages (and the people developing these methods) for single-cell data analysis, including RNA-seq, ATAC-seq, etc. Contributions welcome...

### Articles of particular interest

- [Uncovering pseudotemporal trajectories with covariates from single cell and bulk expression data](http://www.doi.org/10.1038/s41467-018-04696-6), Campbell and Yau 2018 ([reading_notes](reading_notes/Uncovering_pseudotemporal_trajectories_campbell.md))





A comparison of single-cell trajectory inference methods -> benchmark of 45 trajectory inference methods (Saelens et al. 2019)



### Books of particular interest

- 

### Reading notes

- 

### Other documents

- <a href="reading_notes/glossary.md">glossary </a>



### R packages



##### Trajectory inference



- **PhenoPath**. Genomic trajectories (pseudotimes) with heterogeneous genetic and environmental backgrounds. https://bioconductor.org/packages/release/bioc/html/phenopath.html.
- **dyno**. Trajectory inference workflow. https://dynverse.org/dyno/index.html.



##### Cancer genome analysis

Segmented copy number data were processed using CNtools package 

The PAM505 classifier, which is also deployed in Genefu Bioconductor package20, makes calls based on the 50 gene centroid correlation distance to LA, LB, Basal, Her2 and normal-like centroids

y Bayesian approach in discovering novel driver bio-markers in aberrant STPs given high-throughput gene expression (GE) data. 



##### Statistics 

Processed IlluminaHuman450Kdata(GEO accession GSE72021,described in ref. 26) from221 tumor samples (171 serous, 18 endometrioid, 14 clear cell, 9 mucinous, and 9 other histologic cancer subtypes) was downloaded, along with accompanying clinical annotations, using the GEOquery package from Bioconductor (

The InfiniumPurify R package (https://cran.r-project.org/web/ packages/InfiniumPurify/index.html) was used to estimate tumor purity levels based on DNA methylation profiles (28).



All computations were carried out on the R platform. Package “fastICA” which implements the iterative FastICA algorithm (15) was used to extract non-Gaussian independent components with logcosh contrast function. Components were subsequently assigned to clusters using the “cluster” package. 

Clusters were visualized with 2d t-SNE using the R package “tsne.” The number of clusters was determined as equal to number of components extracted at each run of ICA. When the number of samples is small comparing to the number of features, which is usually the case for biological data, it is convenient to retrieve as many as independent signal sources as possible, and the number of compo- nents extracted is equal to sample size

https://www.bnlearn.com/   bnlearn - an R package for Bayesian network learning and inference  with books "Bayesian networks in R" and "Bayesian networks with examples in R"



The clusters were tested for functional enrichment with a hypergeometric test and categories from Gene Ontology (www.geneontology.org/), KEGG (http://www.genome.jp/kegg/), and Reactome databases (http://www.reactome.org/), using the topGO [99]with a ‘weight01’ algorithm, Category and ReactomePA libraries, respectively

 https://github.com/bioFAM/MOFA. 



The package “pcaMethods” was used to calculate principal components for comparison with independent components. 



Cluster reliability score (CRS) The CRS was introduced in (Alvarez et al., 2018) as a statistically sound way to assess the fit of each sample within a cluster. For each sample, a distance vector V1, representing its distance from all other samples in the same cluster and a vector V2, representing its distance from all other samples in the cohort are computed. The sample distance matrix was computed by taking the weighted VIPER scores for each sample (VIPER activity values multiplied by each MR’s MOMA Score) and calculating the pairwise Pearson correla- tions. The normalized enrichment score of V2 distances, ranked from the largest to the smallest one, in V1 distances, is then assessed using aREA. This produces a p-value that represents the tightness and separation of the cluster being considered in relation to all other samples. A cluster-wide reliability score for each cluster is assessed as the average cluster reliability (NES) of each sample in the cluster, scaled between 0 and 1. Finally, the reliability of the entire clustering solution (global cluster reliability score) is assessed as the average of the cluster-wide reliability score of all clusters in the solution.

##### Vizualization

Sankey diagrams in R https://christophergandrud.github.io/networkD3/#sankey

* treemap using the treemap library

Heatmaps were drawn using the superheat (https://github.com/rlbarter/superheat) and ComplexHeatmap R package

##### Single-cell

R package Scater



recount3, a resource consisting of over 750,000 publicly available human and mouse RNA
sequencing (RNA-seq) samples uniformly processed by our new Monorail analysis pipeline. To facilitate
access to the data, we provide the recount3 and snapcount R/Bioconductor packages as well as
complementary web resources. Using these tools, data can be downloaded as study-level summaries or
queried for specific exon-exon junctions, genes, samples, or other features. Monorail can be used to
process local and/or private data, allowing results to be directly compared to any study in recount3. Taken
together, our tools help biologists maximize the utility of publicly available RNA-seq data, especially to
improve their understanding of newly collected data. recount3 is available from http://rna.recount.
bio



##### Gene annotation

`disgenet2r` is an R package to query and expand DisGeNET data (www.disgenet.org), and to visualize the results within R framework. The disgenet2r is designed to query data for DisGeNET v7.0 (May, 2020).



##### Gene set and pathway enrichment

Gene Set Enrichment Analysis (GSEA) implemented in the WebGestaltR R-package (Liao et al., 2019) was used to infer signatures of approved drugs (D1, 1,202 gene sets) and kinase inhibitors (D2, 1,220 gene sets) available in the drug signatures database (DSigDB; Yoo et al., 2015; http://dsigdb.tanlab.org/DSigDBv1.0/). Based



We used FunRich to identify the KEGG pathways that were
significantly associated with the genes



RegEnrich: An R package for gene regulator enrichment analysis reveals key role of ETS transcription factor family in interferon signaling

aREA analysis The analytic Rank-based Enrichment Analysis (aREA) was introduced in (Alvarez et al., 2016) as an analytical methodology to assess gene set enrichment analysis statistics, producing results that are virtually identical to GSEA (Subramanian et al., 2005) without the need for time-consuming sample or gene shuffling



, there exists single sample gene set enrichment analysis [57] techniques such as gene set varia- tion analysis (GSVA [27]) and fast gene set enrichment analysis (FGSEA [53]) to estimate enrichment score for each TR in a given sample.

fgsea R package

GSVA, a non-parametric, unsupervised technique is used
to estimate TR regulon enrichment scores as a function of genes inside and outside the regulons analogously to a competitive gene settest[27].We use the ‘gsva’ function in the ‘gsva’ package

Gene set analysis was performed using Camera (Wu and Smyth, 2012
) from the limma R package. The Hallmark and cellular components gene ontology gene sets were obtained from the Molecular Signatures Database (Liberzon et al., 2015
) (downloaded June 25th, 2019) and used as the reference gene sets.
Signaling pathway impact analysis (Tarca et al., 2009
) was used to identify the list of KEGG pathways over-represented by the subgroup-specific differentially expressed genes (DEGs) using the iLINCS web-based platform (http://ilincs.org). Only pathways with adjusted p value (SPIA pajd) less than 0.05 were considered, and the status of each enriched pathway was only considered if the topology p value (Top pval) was less than 0.05.
ESTIMATE analysis was performed using the estimate R package (Yoshihara et al., 2013
). Log2-TPM values were used as the input, and the analysis was run using default parameters.
Transcription factor enrichment analysis was performed using the ChEA3 (Keenan et al., 2019
) web portal (https://maayanlab.cloud/chea3/) using default parameters and the subgroup-specific DEGs.



Gene set variation analysis was performed using the GSVA R package (Guinney and Castelo, 2019; Hänzelmann et al., 2013
) using log2-TPM values from protein-coding genes. T



##### Pathway and gene network analyses



To compare two or more networks simultaneously, we developed BioNetStat, a Bioconductor package with a user-friendly graphical interface. **BioNetStat** compares correlation networks based on the probability distribution of a feature of the graph (e.g., centrality measures).



 We develop an R packagePATHcrosstalk (available from GitHubhttps://github.com/fabotao/PATHcrosstalk) with which users can discoverpathways of interest with crosstalk effect considered



MOSClip, an R package implementing a new topo-
logical pathway analysis tool able to integrate multi-
omic data and look for survival-associated gene
modules

Pathway analysis was performed using the R package XGR[35](https://www.nature.com/articles/s41591-019-0734-6#ref-CR35) using the GOBP and Reactome databases. Induced and suppressed  transcripts were analyzed separately against the background of all  tested transcripts

ReactomePA https://yulab-smu.top/biomedical-knowledge-mining-book



e we propose an ap- proach to systematically identify the set of active recep- tor-mediated signaling pathways within any given cell, by combining PPI and gene expression data. This method is implemented as an open source packaging using the ‘R’ programming language. This open source software is called SPAGI (Signaling

iPAGE for pathway analyses



scenic, a computational method for simultaneous gene regulatory network reconstruction and cell-state identification from single-cell rna-seq data (http://scenic. aertslab.org).



Multi-omics Master- Regulator Analysis (MOMA). MOMA integrates gene expression and genomic alterations profiles to identify MR proteins and MR modules that represent the key effectors of a tumor’s mutational landscape and are thus responsible for implementing its associ- ated cancer cell identity.



The full code base for PathTurbEr is available in here: https://github.co m/Akmazad/PathTurbEr/ 



[PaxtoolsR: pathway analysis in R using Pathway Commons]



hMARINa ssMARINa MARINa was run via the VIPER R package (http://www.bioconductor.org/packages/release/bioc/html/viper.html)(Alvarez et al., 2016); and hMARINa was per- formed by extending the functionality of the package



marina to infer master regulator



e a method that compares multiple networks ofunlimited size at the level oflinks and nodes. Our novel method, CoDiNA (Co-expression Differencial Network Analysis), is implemented as an R package that also includes an interactive tool for network visualization



This report proposes a framework for integrating known genetic pathways into a differential network analysis
of two populations. The framework allows any association measure to be used, and a general measure for differen- tial connectivity is considered. Statistical significance is evaluated through a permutation testing procedure. The methodology is implemented in R and is available on GitHub at https://github.com/tgrimes/dnapath



 six additional mutual information methods in the MINET R package (ARACNE, CLR, MIM, MINET, MRNET, MRNETB) [[43](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007241#pcbi.1007241.ref043)

We used ARACNe [100] to infer edges between the hubs and the expressed genes. 

VIPER has been extensively validated as an accurate method-
ology to measure a protein’s activity, on the basis of the enrichment of its tissue-specific activated and repressed tran- scriptional targets (regulon) in over and under-expressed genes (Alvarez et al., 2016)—i.e., akin to a highly multiplexed gene- reporter assay.

The metaVIPER approach (Ding et al., 2018), available from the R VIPER package, was then used to generate two interactomes from the TCGAPRAD cohort (this manuscript) and the 2015 SU2C metastatic Castration Resistant Prostate Cancer (mCRPC) cohort (Rob-

 corto: a lightweight R package for gene network inference and master regulator analysis 



We present a new R package, CorDiffViz, that facilitates the estimation and visualization of differential  correlation networks using multiple correlation measures and inference  methods. The software is implemented in R, HTML and Javascript, and is available at https://github.com/sqyu/CorDiffViz. Visualization has been tested for the Chrome and Firefox web browsers. A demo is available at https://diffcornet.github.io/CorDiffViz/demo.html.



##### Gene expression data query and processing

We inferred the similarity between individual cells in the samples and cell lines by calculating the pairwise Pearson correlation matrix C ={cor(i,j)} between any cell i from the Drop- seq experiments with any cell line, j in the CellAtlas leveraging the R package dmatch

sgnesR (Stochastic Gene Network Expression Simulator in R) is an R package that provides an interface to simulate gene expression data from a given gene network using the stochastic simulation algorithm (SSA)

Stemness scores were calculated as previously described (Malta et al., 2018). Firstly, we used MoonlightR (Colaprico et al., 2020)to query, download, and preprocess the pluripotent stem cell samples (ESC and iPSC) from the Progenitor Cell Biology Consortium (PCBC) dataset



Recount2 - Rpackage One issue that can be encountered when planning DEA of TCGA data is  the fact that some projects on the GDC portal do not contain normal  control samples for the comparison with the tumor samples. As explained  previously, it is now possible to query data from the *Recount2* platform to increase the pool of normal samples and apply the DEA pipelines of *TCGAbiolinks* (see [Fig 4A](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006701#pcbi-1006701-g004) for a workflow). -> download GEO and GTEx



GEO2R -> download and can provide fold changes



 from TCGA obtained via the R Package curatedTCGA

##### Epigenetic data query and processing

EPIC array data from JHU samples were converted to 450K meth- ylation profiles using the convertArray function in the minfi pack- age,



##### Proteomic data query and processing

Proteomic data were accessed and downloaded using the R package “TCGA-Assembler 2” (27, 28) from the CPTAC





a novel technique, sparse multiple canonical correlation network analysis
(SmCCNet), for integrating multiple omics data types along with a quantitative phenotype of inter-
est, and for constructing multi-omics networks that are specific to the phenotype. 

The SmCCNet algorithm is written in R, and is freely available on
the web at https://cran.r-project.org/web/packages/SmCCNet/index.html.



[molnet](https://cran.r-project.org/package=molnet) v0.1.0: Implements a network analysis pipeline that enables integrative analysis of multi-omics data including metabolomics. It allows for  comparative conclusions between two different conditions, such as tumor  subgroups, healthy vs. disease, or generally control vs. perturbed. The  case study presented in the [vignette](https://cran.r-project.org/web/packages/molnet/vignettes/Molnet_Vignette.html) uses data published by [Krug (2020)](https://www.cell.com/cell/fulltext/S0092-8674(20)31400-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867420314008%3Fshowall%3Dtrue).



### Python packages

 https://github.com/bioFAM/MOFA. 

https:// github.com/BeautyOfWeb/Multiview-AutoEncoder 



https://github.com/raphael-group/hotnet2

HotNet2 identifies subnetworks of a protein-protein interaction network with more mutations ("heat") than expected.



implemented the pypgatk package and the pgdb workflow to create proteogenomics databases based on ENSEMBL resources

pypgatk: (https://github.com/bigbio/py-pgatk/), 



siVAE is implemented as a Python package and is available from PyPi (https://pypi.org/project/siVAE

I. We used LDVAE29 and scVI21 implemented in SCANPY76 package available from PyPi



Two separate Python packages were used to compute neural network feature attributions in our experiments. We used the DeepExplain Python package that implemented all feature attribution methods (Saliency Maps, Grad*Int, DeepLIFT, IntGrad, Shapley Value) included in our experiments in reverse-mode77. We used the tensorflow-forward- ad Python package for computing Saliency Maps and Grad*Int in forward-mode78

Gene Relevance, we used the published R package45. The method required the latent embeddings learned from siVAE as well as the raw count data corresponding to the embeddings.



GSEApy: Gene Set Enrichment Analysis in Python.



### Other tools

##### clustering and heatmap visualization 

(https://software.broadinstitute.org/GENE-E/index.html). Pearson

https://software.broadinstitute.org/morpheus/

Mor- pheus software (https://software.broadinstitute.org/ morpheus, Broad Institute). Resulting heatmap plots were created in GraphPad Prism software (La Jolla, California



##### General figures

bioRender.com



##### Gene set enrichment



GO enrichment by Gorilla (http://cbl- gorilla.cs.technion.ac.il) 



e screened for enrichment in MSigDB and chromosome positions using the hypergeometric testing implemented in the GeneListEnrich shiny tool (https://github.com/aleferna/BCLandscape_Shiny/tree/master/GeneListCompare)



Gene ontology (GO) analyses were performed using ToppGene (https//topgene.cchmc.org)(

functional characterization of clustering results by single sample Gene Set Enrichment Analysis (ssGSEA)  https://github.com/broadinstitute/ ssGSEA2.0



Pathway enrichment results were imported into Cytoscape (Shannon et al., 2003) using the Enrichment Map app (Merico et al., 2010) for network analysis of pathways.

GO semantic similarity between genes can be calculated using tools GOSemSim (Yu et al. 2010) or GOssTo (Yu et al. 2010)

The biological pathways for the genes was performed using ToppFun software of ToppGene suite (Chen et al., 2009). ToppGene is a one-stop portal for gene list enrichment analysis and candidate gene prioritization based on functional annotations and protein interactions network. ToppFun detects functional enrichment of the provided gene list based on transcriptome, proteome, regulome (TFBS and miRNA), ontologies (GO, Pathway), phenotype (human disease and mouse phenotype), pharmacome (Drug-Gene associations), literature co-citation, and other features. The biological pathways with FDR < 0.05 were considered significantly affected.



FunRich 3.1.3 (http://www.funrich.org) is an independent software tool, used mainly for the functional enrichment and interaction network analysis of genes and proteins



Currently, no software enables non-developers in accessing specific  pathway information. To fill this gap, we present BioPAX-Parser (BiP), a software tool to simply and efficiently perform access to the pathway  information encoded in BioPAX Level 3 (coming from any available public  or private database).



see tools mentioned in Gene Regulatory Network Inference: An Introductory Survey Huynh-Thu & Guido Sanguinetti 2019

DAVID

ConsensusPathDB allows us to perform overexpression analysis on top of differentially activated MRs to identify signif- icantly enriched molecular functions (M), cellular components (C), biological process (BP), pathways (P) and protein complexes (PC). The advantage of using ConsensusPathDB over a popular tool like DAVID [31]isthatitprovidesthe option to search through multiple databases (different types of interactions) to find enriched pathways, unlike DAVID, which only uses the KEGG database. Moreover,

Gene ontology, pathway and cytoband analysis is performed by ToppGene (https ://toppg ene. cchmc .org/)

**Deconvolution**

* CIBERSORT RNA-based estimates of overall I-TME provided by CI- BERSORT absolute scores (immune cell infiltration from RNA data)

The RNA-based tumor microenvironment inference tool ESTIMATE (Yoshihara et al., 2013) was used to derive the overall immune score and stromal score for each sample. 

Weused an established methylation-based deconvolution method, EDec (Onuchic et al., 2016) to dissect the composition of different cell types within the whole bulk tumor. 



##### Data visualization

* A website for interactive visualization of the multi-omics dataset is available at: http://prot-shiny-vm.broadinstitute.org:3838/ CPTAC-BRCA2020
* Lolliplots for CDKN2A mutations in the CPTAC head and neck squamous cell carcinoma and this lung squamous cell carcinoma cohort were generated using the ProteinPaint web application to visualize mutations (Zhou et al., 2016). Mutations



##### Networks and pathways

http://www.pathwaymapper.org/

gephi for graph vizualization



visualization: Systems Biology Graphical Notation (SBGN) project, an effort to standardise the graphical notation used in maps of biological processes. https://sbgn.github.io/



GO-Elite Pathway Analysis http://www.genmapp.org/go_elite/
GO-Elite is designed to identify a minimal non-redundant set of biological Ontology terms or pathways to describe a particular set of genes or metabolites. Default resources include multiple ontologies (Gene, Disease, Phenotype), pathways (WikiPathways, KEGG, Pathway Commons), putative regulatory targets (transcription, microRNA, domains) and cellular biomarkers. Multiple options for pathway visualization are also available. This software can be run using an intuitive graphical user interface, through command-line arguments, using a web interface and through the program GenMAPP-CS. For program details and to get answers to common questions, check out the Manual, FAQ, Tutorials or User Group.u



PathVisio https://pathvisio.github.io/
Biological pathway creation and curation softw areu

CASCADE_SCAN generates a specific pathway for a list of protein molecules using a steepest descent method. That is, the method takes the input proteins and then finds their interaction partners iteratively based on some evidences



##### Genomics and genomics data processing



We integrated somatic mutation, CNV, DNA methylation, RNA, protein, phosphorylation (phospho) and acetylation (acetyl) levels via iProFun (Song et al., 2019) to investigate the functional impacts of DNA alterations in GBM.

PROGENy (Schubert et al., 2018) was used to generate activity scores for EGFR based on RNA expression data.

We used iProFun, an integrative analysis tool to identify multi-omic molecular quantitative traits (QT) perturbed by DNA-level varia- tions. 



AltAnalyze identities alternative splicing events,impacted protein isoforms, [domain composition](http://altanalyze.readthedocs.io/en/latest/DomainAnalysis) and [microRNA targeting](http://www.altanalyze.org/help.htm#microrna).  http://www.altanalyze.org/

##### Proteomics data processing



* **PANOPLY**: a cloud-based platform for automated and reproducible proteogenomic data analysis. -> encapsulate the complex data processing required for proteogenomics  and provide a simple interface to deploy a range of algorithms developed for data analysis, we have developed PANOPLY—a cloud-based Platform for Automated aNd reprOducible Proteogenomic data anaLYsis

* the **Proteomics-Genomics DataBase** (**pgdb** - https://nf-co.re/pgdb) workflow

* **Biofactoid** (biofactoid.org), a web-based software system that
  empowers authors to capture and share structured human- and machine-readable summaries of molecular-level interactions described in their publications

ProTExA is a web-tool that provides a post-processing workflow for the  analysis of protein and gene expression datasets. Using network-based  bioinformatics approaches, ProTExA facilitates differential expression  analysis and co-expression network analysis as well as pathway and  post-pathway analysis.



##### Data processing and statistics 

* **SignatureAnalyzer** exploited the Bayesian variant of the NMF algorithm and enabled an inference for the optimal number of signatures from the data itself at a balance between the data fidelity (likelihood) and the model complexity (reg-

* **Customprodbj** (Wen et al., 2020)(https://github.com/bzhanglab/customprodbj) for customized database construction.







##### Cancer genome analysis

* **BreakDancer** (Chen et al., 2009) and **Meerkat** (Yang et al., 2013) algorithms were used to detect structural variations.

* **PARADIGM** Integrated Pathway Analysis Integrated Pathway Levels (IPLs) mRNA expression, SCNA, and pathway interaction data for 80 UM samples were integrated using the PARADIGM software (Sedge- wick et al., 2013). Briefly, this procedure infers integrated pathway levels (IPLs) for genes, complexes, and processes, using pathway interactions, and genomic and functional genomic data from each patient sample. Normalized

* interactive online resource with navigable proteomics, transcriptomics, and drug sensitivity profiles at https://lehtio-lab.se/forall/.

* fusion events (FUSs) reported by Pipeline for RNA-Sequencing Data Analysis (**PRADA**) (Torres- García et al., 2014)

* **OncoLnc**: [oncolnc.org](http://www.oncolnc.org/) online tool: link TCGA survival data to mRNA, miRNA, or lncRNA expression level



. We introduce the Molecular Oncology Almanac (MOAlmanac), a paired  clinical interpretation algorithm and knowledge base to enable  integrative interpretation of multimodal genomic data for point-of-care  decision making and translational-hypothesis generation. 



PROGENy, a method that overcomes both limitations by leveraging a large compendium of publicly available perturbation experiments to yield a common core of Pathway RespOnsive GENes. Unlike pathway mapping methods, PROGENy can (i) recover the effect of known driver mutations, (ii) provide or improve strong markers for drug indications, and (iii) distinguish between oncogenic and tumor suppressor pathways for patient survival. https://www.nature.com/articles/s41467-017-02391-6 Schubert et al. 2018



GENE-E package (https://software.broadinstitute.org/GENE-E/download.html) was used to cluster the 1× WGS samples according to the similarity of their CCF profiles using Pearson correlation.
Reporting



https://ecotyper.stanford.edu;/  EcoTyper as a broadly applicable framework for
high-throughput identification of cell states and multicellular

communities from primary tissue specimens. It consists of three

key steps: digital purification of cell-type-specific gene expres-
sion profiles from bulk tissue transcriptomes, identification and
quantitation of transcriptionally defined cell states, and co-

assignment of cell states into multicellular communities 



##### Pathway and gene network analyses



**GraphWeb** is a public web server for graph-based analysis of biological  networks  https://biit.cs.ut.ee/graphweb/welcome.cgi?t=help 

Visualization of cellular processes and pathways http://newteditor.org



We have built a microarray data analysis tool, named PATIKAmad, which can be used to associate microarray data with the pathway models in mechanistic detail, and provides facilities for visualization, clustering, querying, and navigation ofbiological graphs related with loaded microarray experiments. PATIKAmad is freely available to noncommercial users as a new module ofPATIKAweb at http://web.patika.org.



pathway visualization tools ChiBE29,30 or Newt.



functional association between genes using data mining and text mining and can be valuable and reliable resources for generating a prior matrix, e.g., STING (Von Mering et al. 2005) and AraNet (Lee et al. 2010)





##### Drawing tools

* biorender.com
* Google drawing
* Canvas
* Vectr
* http://www.celldesigner.org
* http://www.cellillustrator.com
* https://bioicons.com/
* https://app.diagrams.net/
* https://smart.servier.com/



### Public data



##### Genomics

* CellAtlas To assign cell types to individual cells, we used a bulk RNA sequencing data-set from 95 cell lines collected by CellAtlas (Mabbott et al., 2013) that covered 33 major cell types in normal human tissue, including common immune, endothelial, epithelial, fibroblast and mesodermal cells. T
* Transcription factor peaks (bed files) processed from cistrome.org
* GENIE GENIE combines genomic and clinical data in an attempt to associate genomic alterations with phenotypic changes
*  GXB GXB compiles immunological transcriptomic data

* http://gent2.appex.kr/ Now, GENT2 contains more than 68,000 samples and has several new useful  functions. First, GENT2 now provides gene expression across 72 different tissues compared to 57 in GENT. Second, with increasing importance of  tumor subtypes, GENT2 provides an option to study the differential  expression and its prognostic significance based on tumor subtypes.  Third, whenever available, GENT2 provides prognostic information of a  gene of interest. Fourth, GENT2 provides a meta-analysis of survival  information to provide users more reliable prognostic value of a gene of interest.

[TANRIC](http://ibl.mdanderson.org/tanric/_design/basic/index.html): focus on noncoding RNA

known associations between disease phenotypes and genes extracted from the Online Mendelian Inheritance in Man (OMIM) database

We generated a single-cell tumor immune atlas, jointly analyzing  published data sets of >500,000 cells from                     217 patients and 13 cancer types, providing the  basis for a patient stratification based on immune cell compositions. (Nieto et al. 2021 - A single-cell tumor immune atlas for precision oncology)

miRTarBase: miRNA target database

publicly available transcriptomic data from the Genomics of Drug Sensitivity in Cancer (GDSC) database

, an additional validation on the set ofeight datasets
(BLCA, BRCA, COAD, GBM, HNSC, LUAD, OV and SKCM can- cers) obtained from the PRECOG repository was conducted. F

##### Epigenomics

* ubiquitination data (Kim et al 2011)
* DNase-Seq data were obtained from ENCODE (http:// hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/ wgEncodeUwDgf, sample ‘wgEncodeUwDgfTregwb7849 5824’) with DNase hypersensitivity (DHS) peaks. 
* Protein- binding footprints (FP) within the DHS peaks were identi- fied by scanning the DHS intervals for gaps in the signals using Wellington [101],

* Characterizing genetic influences on DNA methylation (DNAm) provides an  opportunity to understand mechanisms underpinning gene regulation and  disease. In the present study, we describe results of DNAm quantitative  trait locus (mQTL) analyses on 32,851 participants, identifying genetic  variants associated with DNAm at 420,509 DNAm sites in blood. http://mqtldb.godmc.org.uk/ 

  [MEXPRESS](http://mexpress.be/): focus on methylation and gene expression

  

##### Proteomics

* Reverse Phase Protein Arrays (RPPAs) data
* protein complexes, specifically the CORUM database
* subcellular localization was retrieved from the Human Protein Atlas [97](www.proteinatlas.org).
* protein/mRNA half-life data: Schwanhäusser et al. 2011 (mice)
* Clinical Proteomic Tumor Analysis Consortium (CPTAC) CPTAC contains mass spectrometry-based proteomic analysis of tumors from TCGA. The aim of CPTAC is to create a proteogenomic resource where dysregulated proteins and phosphorylation sites can be identified and potentially connected to genomic alterations. https://cptac-data-portal.georgetown.edu/cptac/s?id=393177
* Proteomics Identification Database (PRIDE). PRIDE aims to be a resource for open access sharing of mass spectrometry data, not just across cancer. They currently have over 9200 datasets available, including 297 breast cancer datasets.
* Human Proteome Organization (HUPO) The human proteome project, run by HUPO aims to identify all the proteins in the human proteome and to begin to assess their functionalities and interactions

* TcellSubC: An Atlas of the Subcellular Proteome of Human T Cells

* All raw data are available via ProteomeXchange with identifier PXD024190.

* The objective of the Human Protein Atlas (HPA) ([www.proteinatlas.org](http://www.proteinatlas.org)) effort is to take advantage of these bioimaging approaches to map the  expression of all human protein-coding genes across all major human  cells, tissues, and organs.

  

  17,511 phosphorylated proteins, 6,928 acetylated proteins and 5,418 methylated proteins from the PhosphoSitePlus database44

a set of protein abundances30 in human cell lines. Beck, M. et al. The quantitative proteome of a human cell line. Mol Syst Biol 7, 549 (2011)



A list of known human E3 ubiquitin and ubiquitin-like ligases and DUBs was compiled from (Medvar et al., 2016; Nijman et al., 2005)



32TCGA studies that provide proteomic and phospho- proteomic measurements of tumor biopsies from various types of cancer patients. Those studies provide RPPA profiles of a total of 7,694 patients using 259 antibodies.



Resources such as the Human Protein Reference Database (HPRD) (Peri et al., 2004; Keshava Prasad et al., 2009), the Biological General Repository for Interaction Datasets (BioGRID) (Stark et al., 2006), and the Search Tool for Retrieval of Interacting Genes/Proteins (STRING) (Szklarczyk et al., 2019) 

Popular phosphorylation focused databases such as PhosphoSitePlus (Hornbeck
et al., 2015), PHOSIDA (Gnad et al., 2007), Phospho.ELM (Diella et al., 2004) and qPhos (Yu et al., 2019) host many phosphosites and act as repositories for both low and high-throughput data.

An extensive review of phosphoproteomics resources can be found here (Savage & Zhang, 2020). 

Moving away from databases acting as phosphosite repositories, many databases holding signalling information in the form of a K-S networks exist such as RegPhos (Huang et al., 2014), PhosphoNet (Safaei et al., 2011) and Phosphonetworks (Hu et al., 2014). Though they are limited to K-S interactions, they have been proven useful in providing phosphoproteomics mechanistic insight either on their own or while integrated into other databases (Rohrs et al., 2018; McGuire et al., 2017; Tong et al., 2019).



The human PPIs were collected from five public databases (i.e., BioGRID20, DIP21, IntAct22, MINT23 and
MIPS24

http://www.proteinatlas.org



ProtMapper identifies valid reference positions with high precision and reasonable recall, making it possible to filter out machine reading errors from text mining and thereby assemble a corpus of 29,400 regulatory annotations for 13,668 sites, a 2.8-fold increase over PhosphoSitePlus, the current gold standard.



* CT Antigen database

NeoFlow (https://github.com/bzhanglab/neoflow) for neoantigen prediction (Wen et al., 2020). Specifically, Optitype (Szolek et al., 2014) was used to find human leukocyte antigens (HLA) for each sample based on WES data. Then we used netMHCpan (Jurtz et al., 2017) to predict HLA peptide binding affinity for somatic mutation–derived variant peptides with a length between 8-11 amino acids. The



Immunohistochemistry-based antibody-specific staining scores in lung tumors were obtained from the Human Protein Atlas (HPA, https://www.proteinatlas.org), in which tumor-specific staining is reported in four levels, i.e., high, medium, low, and not detected. 



Remaining variant peptides were further filtered using PepQuery (http://www.pepquery.org)(Wen et al., 2019) with the p value cutoff of 0.01. Competitive filtering based on unrestricted posttranslational modification searching was enabled in PepQuery validation. The spectra of variant peptides were annotated using PDV (http://pdv.zhang-lab.org)(Li et al., 2019b



Cancer/testis (CT) antigens were downloaded from the CTdatabase (Almeida et al., 2009). C

PTMsigDB 	PTM signature (for phosphorylation-driven signature analysis)

PhosphoSitePlus phosphorylation site annotation

Signor  phosphorylation site annotation





##### Gene sets, gene ontologies, gene signatures, gene pathways and gene networks

* MSigDB [93]website (http://software. broadinstitute.org/gsea/msigdb/) the v5.1 C7 (‘immuno- logic signatures’) collection

* TF genes identified with ei- ther or both of two alternative annotations: (1) the human genes with a symbol annotated with the term ‘GO: 0003700’ in the Gene Ontology Consortium database (www.geneontology.org) or (2) the Ensembl gene ID re- trieved by querying the BioMart service (http://grch37. ensembl.org/) with the Gene Ontology ID ‘GO:0003700

Finally,we perform downstream analysis of the MRs specific to ICR-L using ConsensusPathDB [35] 



Each set of marker genes were systematically tested for enrichment in a chemical and genetic perturbation sig- nature (the CGP collection) curated by the Molecular Signature Database (MSigDB)





To derive association scores from the Open Targets re- source (www.opentargets.org) the version from September 2016 was used

* unsupervised pathway enrichment analysis using **Reactome**



The Database of Interacting Proteins (DIP: http://dip.doe-mbi.ucla.edu) is a database that documents experimentally determined protein–protein interactions. https://doi.org/10.1093/nar/30.1.303

WikiPathways6, Ingenuity Pathway Analysis (qiagen.com), KEGG7, and Reactome8



Tovar et al. 2015: An important step for this algorithm is the selection of the Transcription Factors, since they will determine the rest of the cal- culation. A proper annotation of transcription factors is crucial for an accurate description of the process under investigation. Here, we used the HGU133A annotation file, in which we found 1142 TFs (Supplementary Material 1). This list was compared with other three lists. Those lists are available in Shimoni and Alvarez (2013), Vaquerizas et al. (2009) and http://www.bioguo.org/AnimalTFDB/ Animal Transcription Factor DataBase, respectively. We want to stress that all four lists show consistency among them

genes that were annotated with a signaling function without receptor domain function from Gene Ontology (GO)45 as well as 5,701 genes that carried a trans-membrane protein domain

95,722 links between 209 human transcription factor and 8,910 human genes from the TRANSFAC42 database as provided by mSigDB




DNA-binding transcription factors28,29 and 501 genes from the Kinome NetworkX database27 that collects kinase information from the literature and other databases



 GRAND (https://grand.networkmedicine.org) as a database for computationally-inferred, context-specific gene  regulatory network models that can be compared between biological  states, or used to predict which drugs produce changes in regulatory  network structure. The database includes 12 468 genome-scale networks  covering 36 human tissues, 28 cancers, 1378 unperturbed cell lines, as  well as 173 013 TF and gene targeting scores for 2858 small  molecule-induced cell line perturbation paired with phenotypic  information. GRAND allows the networks to be queried using phenotypic  information and visualized using a variety of interactive tools.





We assessed the overlap of these prior relations with the ‘‘canonical pathways’’ gene sets in MSigDB to understand its coverage. This collection has 2,815 gene sets curated from the databases BioCarta, KEGG, NCI-PID, Reactome, and WikiPath- ways. 



pathway enrichment analysis, using Kyoto Encyclopedia of Genes and Genomes (KEGG), Biocarta, and Reactome pathways from the Molecular Signature Database (MSigDB) (

This was further complemented by incorporating Kyoto Encyclopedia of Gene and Genomes (KEGG) ([Ogata et al., 1999](https://www.sciencedirect.com/science/article/pii/S1476927119303871?casa_token=TTrIGf-vCaAAAAAA:w4xLI5CZzXwhQlJZkG0oL7wdOFfgYDNGK7CwLXnGaxygNkR6XNnWmiDSJkkiwx33yBWzwqcbkuVd#bib0180)) PFAM protein domains ([Sonnhammer et al., 1998](https://www.sciencedirect.com/science/article/pii/S1476927119303871?casa_token=TTrIGf-vCaAAAAAA:w4xLI5CZzXwhQlJZkG0oL7wdOFfgYDNGK7CwLXnGaxygNkR6XNnWmiDSJkkiwx33yBWzwqcbkuVd#bib0240)) and INTERPRO protein domains ([Mitchell et al., 2014](https://www.sciencedirect.com/science/article/pii/S1476927119303871?casa_token=TTrIGf-vCaAAAAAA:w4xLI5CZzXwhQlJZkG0oL7wdOFfgYDNGK7CwLXnGaxygNkR6XNnWmiDSJkkiwx33yBWzwqcbkuVd#bib0160)) information of the selected genes.

list of transcription factors in the TFCheckpoint curated database (Tripathi et al., 2013).

IntAct (Orchard et al., 2014), SIGNOR (Licata et al., 2019; Perfetto et al., 2016) and Reactome (Fabregat et al., 2018)], and

Searching for these patterns in Pathway Commons generated 28,517 prior relations in four different types (listed below). To increase coverage, we added relations from several other databases (PhosphoNetworks,13 iPTMnet,14 TRRUST,15 and TFactS),16 which are not in Pathway Com- mons, and increased our relationships to 39,232:



Institute for Systems Biology Regulome Explorer (http://explorer.cancerregulome.org),

databases that characterize pathways such as NetPath [13], KEGG [14], Reactome [11], and dozens of others.

Repositories such as WikiPathways [17] and Path- wayCommons [23] now contain thousands of pathways comprised of millions of interaction

aggregation of publicly available molecular interactions and biological pathway databases provided by the Pathway Commons (PC) resource.23 The aggregated data is represented in the standard Biological Pathway Exchange (BioPAX) language and provides the most complete and rich representation of the biological network models stored in PC. These complex biochemical reactions were reduced to pairwise relationships using rules to generate a Simple Interaction Format (SIF) representation of BioPAX interactions. The reduction of BioPAX interactions to the SIF allows for the representation of pairwise molecular interactions in the context of specific binary relationships. The



Significant pathway annotations (FDR < 0.05) from the Panther over-representa- tion test (database: GOBP complete, http://www.pantherdb. org/) were used to annotate each cluster



we examined genomic-based patterns of oncogenic pathway activity, the tumor microenvironment and other important features in human breast tumors using a panel of 52 previously published gene expression signatures (Supplementary Table - An integrated genomics approach identifies drivers of proliferation in luminal-subtype human breast cancer)



a pathway enrichment analysis of Hallmark, KEGG, and Reactome. The

Genes and Genomics (KEGG) (Kanehisa & Goto, 2000), Reactome (Joshi-Tope et al., 2005), WikiPathways (Slenter et al., 2018), and the SIGnaling Network Open Resource (SIGNOR) (Perfetto et al., 2016) have

Pathways were obtained in BioPax Level 3 format, and included the NCIPID and BioCarta databases from http://pid.nci.nih.gov
and the Reactome database from http://reactome.org.

We queried the PTM signatures database (PTMsigDB) v1.9.0 downloaded from https://github.com/broadinstitute/ssGSEA2.0/ tree/master/db/ptmsigdb using the flanking amino acid sequence (+/? 7 aa) as primary identifier. We used the implementation of PTM-SEA available on GitHub (https://github.com/broadinstitute/ssGSEA2.0) using the command interface R-script (ssgsea- cli.R). The





Creating a Curated Transcription Factor (TF) Regulome A compendium of TFs and their targets (TF regulons) were created by combining information from four databases:
(i) SuperPathway (Sedgewick et al., 2013): This is the same interaction network used in the PARADIGM analysis (above). Only links that correspond to regulation at the transcriptional level were retained for MARINa and hMARINa use.
(ii) Literome (Poon et al., 2014): The network was filtered to include only transcription links in which the regulator is a known TF. (iii) Multinet (Khurana et al., 2013): The network was reduced to links that correspond to regulation on transcriptional level. (iv) ChEA (Lachmann et al., 2010): Data from the Gene Expression Atlas (Petryszak et al., 2014) was used to filter the inferred links in the ChEA database. Specifically, the context likelihood of relatedness (CLR) method (Faith et al., 2007) was used to compute a measure of association between every pair of genes. The top 10% of gene pairs based on the CLR score were intersected with the ChEA network and the overlapping pairs were added to the final combined network.
The combined network includes 72,915 transcriptional regulatory links between 6,735 regulators and their targets. Only regulators
with at least 15 targets were considered in the final analysis, which resulted in a final network consisting of 419 TFs with 58,363 total targets (covering a set of 12,754 unique targets). 



Creating a Curated Kinase Regulome Proteins identified as kinases in Manning (Manning et al., 2002) or Uniprot (UniProt Consortium, 2014) were aggregated into a list of 546 kinases. Protein substrates were extracted from PhosphositePlus (Hornbeck et al., 2014) on March 7, 2015. Kinase-substrate interactions were retained if the kinase appeared in the Manning-Uniprot kinase list and the kinase was identified as a human protein in the PhosphositePlus database. The final compendium consisted of 5,388 links between 342 kinases and 2,260 unique substrates.





upstream regulator analysis from the commercial Ingenuity Pathway Analysis software.



Pathway Commons (https://www.pathwaycommons. org) is an integrated resource of publicly available information about biological pathways including bio- chemical reactions, assembly of biomolecular com- plexes, transport and catalysis events and physical interactions involving proteins, DNA, RNA, and small molecules 

www.pathguide.org



see tables in: The status of causality in biological databases: data resources and data retrieval possibilities to support logical modeling; Briefings in Bioinformatics, 22(4), 2021, 1–15 https://doi.org/10.1093/bib/bbaa390



see tools mentioned in Computational methods to dissect gene regulatory networks in cancer Iyer et al. 2017

see tools in Computational methods for Gene Regulatory Networks reconstruction and analysis: A review Delgado et al. 2019 (AI)

see ref in A census of pathway maps   in cancer systems biology Kuenzi and Ideker 2020 https://www.nature.com/articles/s41568-020-0240-7.pdf



Here, we worked within the framework of the TCGA PanCancer Atlas initiative ([Cancer Genome Atlas Research Network et al., 2013c](https://www.sciencedirect.com/science/article/pii/S0092867418303593#bib20)) to build a uniformly processed dataset and a unified data analysis  pipeline aimed at exploring similarities and differences in canonical  cancer pathway alterations across 33 cancer types. Oncogenic Signaling Pathways in The Cancer Genome Atlas Sanchez-Vega et al. 2018



Biocarta, CORUM, Innate DB, KEGG, WikiPathways, Reactome, Nepath, PIC and PINdb, all of which are available in Consensus- Pathdb, for

##### Protein-protein interactions

PrePPI (http://bhapp.c2b2.columbia.edu/PrePPI) is a database that  combines predicted and experimentally determined protein-protein  interactions (PPIs) using a Bayesian framework.

* Human Interactome database: PPI 
* iRefIndex: PPI

* databases  that curate both categories of interactions for humans and other model organisms: BioGrid, IntAct, HPRD, iRefWeb, DIP, MINT, MIPS and VisAnt 

DEPOD
CORUM
Signor2
Reactome

DIP PPI https://dip.doe-mbi.ucla.edu/dip/Main.cgi [11]
MINT PPI https://mint.bio.uniroma2.it/ [12]
IntAct PPI https://www.ebi.ac.uk/intact/ [13]
BioGRID PPI https://thebiogrid.org/ [14]
STRING Various https://string-db.org/ [15]
COXPRESdb Co-expression https://coxpresdb.jp/ [16]
SEEK Co-expression http://seek.princeton.edu/ [17



e PPIs, we used the high confidence human PPIs (Karagoz et al., 2016), comprising 147,923 interactions among 13,213 proteins. Karagoz and coworkers assembled and integrated physical PPIs of Homo sapiens from six publicly available databases including BioGRID (Chatr-Aryamontri et al., 2015), DIP (Salwinski, 2004), IntAct (Orchard et al., 2014), HIPPIE (Schaefer et al., 2012), HomoMINT (Persico et al., 2005), and HPRD (Prasad et al., 2009). 



PPIs from the open-access IntAct database which adopts a merging algorithm and a scoring system to provide richly annotated molecular interaction data. IntAct



Database Species Number of entries URL Reference
BIND About 1500 species 200 000 interactions http://bind.ca [26]
DIP 10 species 25 612 proteins and 75 400 interactions http://dip.doe-mbi.ucla.edu [27]

MINT 6 species 35 511 proteins and 241 458 interactions http://mint.bio.uniroma2.it/mint [28]

IntAct Over 15 species 65 200 proteins and 312 217 interactions http://www.ebi.ac.uk/intact [29]

BioGRID 51 species 471 829 interactions http://thebiogrid.org [30]

HPRD Homo Sapiens 30 047 proteins and 41 327 interactions http://hprd.org [31]

STRING Over 1100 species 5 million proteins and 200 million interactions http://string-db.org [32]

HIPPIE Homo Sapiens 72 916 interactions http://cbdm.mdc-berlin.de/tools/hippie [33]





public databases have cataloged networks of known and predicted PPIs, such as STRING (Szklarczyk et al., 2019), IntAct (Orchard et al., 2014), CellCircuits (Mak et al., 2007), and PINA (Cowley et al., 2012) [more comprehensive lists are described by Huang et al. (2018) and Miryala et al. (2018)].



PPI network from PINA 2.0 contained approximately 16,000 nodes
and 170,000 interactions. Data from PINA comes from six different databases: IntAct, MINT, BioGRID, DIP, HPRD, and MIPS MPact [32].



The database-driven network graph is taken from the GeneMania
database (https:// GeneMania.org/) [36]. GeneMania has a large num- ber of interactions and incorporates both gene-to-gene and protein-to- protein interactions. Since

some tissue- and cell-specific interaction databases are now available such as TissueNet (Basha et al., 2017), the Integrated Interactions Database (Kotlyar et al., 2019), and HumanBase (Greene et al., 2015). Tissue-specificity



* Human Protein Reference Database (HPRD): human protein-protein interactions that are manually extracted from the literature by expert biologists;  one of the most com- prehensive protein interaction databases, only covers less than half of human protein-coding genes

* Biological General Repository for Interaction Datasets (BioGRID): protein and genetic interactions of major model organism species

* Biomolecular Interaction Network Database (BIND): both high-throughput and manually curated interactions between biological molecules

* IntAct: molecular interaction database (IntAct) protein-protein interactions derived from literature 

* Molecular INTeraction database (MINT): information about physical interactions between proteins

* STRING: catalog of protein–protein interactions across species and conditions into databases 

* BioGRID: protein interaction data 

* Biological General Repository for Interaction Datasets (BioGRID, https://thebiogrid.org/), Saccharomyces Gen- ome Database (SGD, http://www.yeastgenome.org/), Hu- man Protein Reference Database (HPRD, http:// www.hprd.org/), Search Tool for the Retrieval of Interacting Genes/Proteins (STRING: http://www.string- db.org/)

* HINT: PPI database; also distinguishes between interactions curated from small-scale studies and those obtained from high-throughput experiments.

* CORUM: known interaction pairs from 

* DEPOD
* CORUM
* Reactome
* Signor2
*  OmniPath 
* study comparing PPI databases found 375 resources (Bajpai et al., 2020). 



3D structures submitted directly by authors to the Protein DataBank (PDB)

##### Cancer-related

review with a huge list of tools and references: Vlachavas et al. 2021

A Detailed Catalogue of Multi-Omics Methodologies for Identification of Putative Biomarkers and Causal Molecular Networks in Translational Cancer Research A Detailed Catalogue of Multi-Omics Methodologies for Identification of Putative Biomarkers and Causal Molecular Networks in Translational Cancer 

A list of targetable drivers was provided by the PCAWG driver group (Rheinbay et al., 2020).

*Datasets*

* Fudan University Shanghai Cancer Center (FUSCC) presented comprehensive clinical, genomic, transcriptomic data of 465 primary TNBC19, which is the largest TNBC genomic project to date

* https://osf.io/gqrz9  ([Tatlow and Piccolo 2016](https://www.nature.com/articles/srep39259)) RNA-seq data from **TCGA** and **CCLE** -> transcript level expression estimates; can be combined to gene level expression estimates using Scater (eg. Campbell and Yau 2018)

* CPTAC datasets
* Oslo2 landscape cohort
* Sweden Cancerome Analysis Network-Breast (SCAN-B)




 Gene expression-based Outcome for Breast cancer Online (GOBO database) [[24](https://www.mdpi.com/2079-7737/10/3/247/htm#B24-biology-10-00247)]. GOBO (http://co.bmc.lu.se/gobo, accessed date 26 February 2021) enables a rapid assessment of gene  expression levels, the identification of co-expressed genes and  association with the outcome for single genes, gene sets, or gene  signatures in an 1881-sample breast cancer data set, generated on  Affymetrix U133A microarrays. 



the PDX Encyclopedia (Gao et al., 2015), in which 399 PDX tumors of varying tissue types had been screened against a total of 40 different monotherapies and 27 combination thera- pies. The genotypes of each PDX had also been established (Gao et al., 2015), which

TIMER (Tumour Immune Estimation Resource) web server is a comprehensive  resource for the systematic analysis of immune infiltrates across  diverse cancer types. (https://cistrome.shinyapps.io/timer/)

TISIDB is also a web portal for tumour and immune system interaction, which integrates multiple heterogeneous data types. (http://cis.hku.hk/TISIDB/index.php) 

Gene Expression Profiling Interactive Analysis (GEPIA) (http://gepia.cancer-pku.cn/index.html) [[20](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04379-y#ref-CR20)], which analyses the RNA sequencing expression from the TCGA and GTEx projects of 9736 tumours and 8587 normal samples



Oncoprotein-specific molecular interaction maps (SigMaps) for cancer network analyses Broyde et al. 2021

. We introduce SigMaps as context-specific networks, comprising modulators, effectors and cognate binding-partners of a specific oncoprotein. SigMaps are reconstructed de novo by integrating diverse evidence sources—including protein structure, gene expression and mutational profiles—via the OncoSig machine learning framework



The multiple myeloma data was obtained from the CoMMpass study by the Multiple Myeloma Research Foundation (MMRF) (https://themmrf.org) and



*Variant and mutation calling* 

* Cancer-related SAAVs are cataloged in the **CanProvar** and **COSMIC** databases.

* 

*Proteomics*

- **MCLP**: Characterization of Human Cancer Cell Lines by Reverse-phase Protein Arrays (https://doi.org/10.1016/j.ccell.2017.01.005)
- 

*Genomics and gene annotation*

* **HCMDB** (Human Cancer Metastasis Database), an integrated database for known metastatic genes

* 
* predicted genes were checked for their differ- ential mRNA expression in oral cancer samples (OSCC/TSCC) in the **Oncomine** database (https ://www.oncom ine.org) and in the **Expression Atlas** (https ://www.ebi.ac.uk/gxa/). 
* **Transciptome Alterations in Cancer Omnibus** (**TACCO**) -TACCO is a resource for identifying differentially regulated transcripts within different cancer types and combining these with survival data to determine prognosis

*Single cell*

- **METABRIC** (proteomics, RNA-seq, ...)
- **TCGA** (RNA-seq, methylation, ...)
- **Protein-gene Expression Nexus**: comprehensive characterization of human cancer cell lines with proteogenomic analysis http:/combio.snu.ac.kr/pen (https://doi.org/10.1016/j.csbj.2021.08.022).
- **CancerSEA**: a cancer single-cell state atlas http://biocc.hrbmu.edu.cn/CancerSEA/ (https://doi.org/10.1093/nar/gky939)
- **CCLE** The Cancer Cell Line Encyclopedia http://www.broadinstitute.org/ccle https://doi.org/10.1038/nature11003
- **Kaplan–Meier Plotter** is an online database that contains comprehensive clinical and microarray data for various cancers

* **OncoKB** annotation

* 

* **Cancer RNA-Seq Nexus** (CRN) database; http://syslab4.nchu.edu.tw/CRN). CRN has a user-friendly web interface  designed to facilitate cancer research and personalized medicine. It is  an open resource for intuitive data exploration, providing  coding-transcript/lncRNA expression profiles to support researchers  generating new hypotheses in cancer research and personalized medicine; a database of phenotype-specific transcriptome profiling in cancer cells

* 

* (Dietrich et al, 2018) a cohort of 200 patient samples of chronic lymphocytic leukaemia, profiled for somatic mutations, RNA expression, DNA methylation and ex vivo drug responses

* Mertins et al. 2016 e quantitative mass-spectrometry-based proteomic and phosphoproteomic analyses of 105 genomically annotated breast cancers, of which 77 provided high-quality data

* https://github.com/mskcc/RNAseqDB -> **TCGA** & GTEX **data** processed for comparing normal vs tumor

* **COSMIC** - COSMIC contains data from over 13 million tumor samples, identifying 6 million coding mutations and over 19 million non-coding mutations. This resource collates all genes implicated in cancer through somatic mutation, of which 719 are currently listed

  

* genomic variants: sources of genomic variants including **COSMIC**, **cBioportal**, **gnomAD**, and mutations detected from sequencing of patient samples

* **PCaDB** is a comprehensive and interactive database for transcriptomes from prostate cancer population cohorts.            We collected 50 transcriptomics datasets with 7,231 samples from public data repositories,            including TCGA, cBioPortal, GEO, and ArrayExpress. A standard bioinformatics pipeline is used to download and            process the expression data and metadata.           PCaDB provides a user-friendly interface for the comprehensive analysis of individual            genes, prognostic signatures, and the whole transcriptomes to elucidate the molecular            heterogeneity in PCa, understand the mechanisms of tumor initiation and progression,            as well as develop and validate prognostic signatures in large independent cohorts
* **MR4Cancer**: a web server prioritizing master regulators for cancer  https://doi.org/10.1093/bioinformatics/bty658 http://cis.hku.hk/MR4Cancer



* (https://ccr.cancer.gov/research/ cancer-moonshot).. The **Human Tumor Atlas Network** (HTAN), part of the National Cancer Institute (NCI) Cancer Moonshot Initiative, will establish a clinical, experimental, computational, and organizational framework to generate infor- mative and accessible three-dimensional atlases of cancer transitions for a diverse set of tumor types. This effort complements both ongoing efforts to map healthy organs and previous large- scale cancer genomics approaches focused on bulk sequencing at a single point in time. Gener- ating single-cell, multiparametric, longitudinal atlases and integrating them with clinical outcomes should help identify novel predictive biomarkers and features as well as therapeutically relevant cell types, cell states, and cellular interactions across transitions

* **TCGA** aggregates an extensive col- lection of omics and clinical datasets from large cohorts of patients for more than 30 types of cancers (24). It also ar- chives histopathology images for solid tumor samples from which omics data were sampled. Currently, more than 24,000 histopathology images are available and can be visualized at the **Cancer Digital Slide Archive** (**CDSA**, http://cancer. digitalslidearchive.net/). In addition, The NCI **Clinical Proteomic Tumor Analysis Consortium** (**CPTAC**) (https://proteomics. cancer.gov/programs/cptac) program also provides high- throughput proteomic data for some of the TCGA tumor spec- imens such as breast cancer, ovarian cancer, and colorectal cancer based on mass-spectrometry technology. 

To facilitate use and dissemination of the data, we have developed a web resource (https://zucchini.gs.washington.edu/BreastCancerProteome/) in which protein abundances can be queried and correlated to genomic and drug sensitivity data, as presented below. 

* The Cancer Imaging Archive (TCIA) database.

* 

https://lincs.hms.harvard.edu/db/datasets/20343; we performed quantitative proteomics on 61 human-derived breast cancer cell lines to a depth of ~13,000 proteins.  All datasets are freely available as public resources on the LINCS portal. 7197 proteins measured in all cell lines

The LINCS project is arguably the most comprehensive systematic perturba- tion dataset currently available, in which multiple cell lines were treated with over tens of thousands perturbagens (e.g., small molecules or single gene knockdowns), followed by monitoring gene expression profiles using a new technology known as the L1000 assay, which utilizes ~1000 (978) landmark genes to infer the entire transcriptome7



* **BioMuta** and **BioXpress**: mutation and expression knowledgebases for cancer biomarker discovery

* Cancer-associated genes (CAG) were compiled using the **Census** website https://cancer.sanger.ac.uk/census. The list of genes is provided in Table S1. -> Satpathy et al. 2021

* The **PreCancer Atlas** (PCA) ofthe NCI envisages a histological and multi-omic mapping strategy in time and space to provide detailedmolecular, cellular, and structural characterization of premalignant lesions and how they evolve to invasivecancers

* **Dependency map** (DepMap) - The goal of the Dependency Map (DepMap) portal is to empower the  research community to make discoveries                    related to cancer vulnerabilities by providing open  access to key cancer dependencies analytical and visualization tools.      DepMap genetic dependency dataset (CRISPR Avana Public 20Q3) that contained 18119 genes and 789 cell lines (https://depmap.org/portal/download/ file: Achilles_gene_effect.csv). O

* **OncoKB** database for oncoprotein and tumor suppressor classification (excluded proteins that have both annotations), and used PhosphoSitePlus and Signor databases for the activating/inhibiting classification of phosphorylation sites. 

* **Omic and Multidimensional Spatial (OMS) Atlas** generated from four serial biopsies of a metastatic breast cancer patient during 3.5 years of therapy. This

* The **Human Tumor Atlas** (HTA) (https://humantumoratlas.org/) 



. For substitutions, structural variants and copy number events, these included a set of genes compiled from the TARGET database from the Broad Institute and multiple sequencing datasets for EAC1

##### Drug drug sensitivity and perturbation datasets

* gene expression and GDSC2 drug sensitivity datasets from the **Genomics of Drug Sensitivity in Cancer** (GDSC) (Yang et al., 2013) . 
* published cell line perturbation experiments from the **Genomics of Drug Sensitivity in Cancer** (GDSC) resource (Iorio et al., 2016; Yang et al., 2013)
* **MOSCATO** trial where druggable genomic aberrations were identified and targeted in patients (Massard et al., 2017

* L1000 and P100 drug connectivity analysis; Level 4 P100 data were downloaded from the **LINCS** Data Portal (Stathias et al., 2019) and were used to calculate drug connec-
  tivities on the phosphoprotein level as previously described (

* **Library of Integrated Network-Based Cellular Signatures** (LINCS) L1000 perturbation-response signa- tures. The scores were computed using the sig_queryl1k_tool pipeline (https://hub.docker.com/u/cmap) and the LINCS L1000 Level 5 compound (trt_cp) signatures from CLUE (https://clue.io, ‘‘Expanded CMap LINCS Resource 2020 Release’’). The

The differentially expressed genes between gene-altered and WT samples were filtered for the 978 genes measured in the L1000
assay and then were processed using the **CLUE** (Subramanian et al., 2017) (summary connectivity score) and **iLINCS** (Pilarczyk et al., 2019) connectivity algorithms. The resulting drug connectivities were aggregated to the compound level using the summary connec- tivity score in CLUE and the Connected Perturbations Z-score in iLINCS. Target annotations for the ranked compounds were extracted from CLUE and iLINCS and combined in a single list. Level





the Comparative Toxicogenomics Database (CTD) (http://ctdbase.org/; accessed on 1 August 2021), which is a free online database that  provides information including chemical–gene/protein interactions. 

drug IC50 data of breast cancer cell lines and the gene expression data  of these cell lines were obtained from the Genomics of Drug Sensitivity  in Cancer database (GDSC) (https://www.cancerrxgene.org/;

* **Genomics of Drug Sensitivity in Cancer** (CRx) resource (Yang et al., 2013)

* essential survival gene datasets from The Cancer Dependency Map, the latter of which catalogs genes driving cancer progression

* druggability information from **DGIdb** (Cotto et al., 2018) and **DEPO** (Sun et al., 2018), 

* **Connectivity Map** (CMAP)   genome-scale library of cellular signatures that catalogs  transcriptional responses to chemical, genetic, and disease  perturbation. To date, the library contains more than 1 Million profiles resulting from perturbations of multiple cell types. - https://clue.io/cmap



review with list of drug responses database: Deep learning for drug response prediction in cancer
Delora Baptista, Pedro G. Ferreira and Miguel Rocha

two large cancer
drug screening resources: the Cancer Therapeutics Response Portal (CTRP) v2 and the Genomics of Drug Sensitivity in Cancer (GDSC) database (Seashore-Ludlow et al., 2015; Yang et al., 2013).



drug-protein interactions from STITCH database



DeepSynergy database (Preuer et al., 2018), in which all pairs of 25 drugs had
been tested across a panel of 39 cell lines (Figure

PRISM (Profiling Relative Inhibition Simultaneously in Mixtures) 
GDSC (Genomics of Drug Sensitivity) 

**Dependency**

DepMap	a compendium of genomics, proteomics, shRNA and CRISPR based gene dependency datasets of hundreds of cell lines.

##### Single cell

* https://github.com/markrobinsonuzh/conquer - ***conquer*** (**con**sistent **qu**antification of **e**xternal **R**NA-seq data sets) repository, which provides access to consistently processed, analysis-ready single-cell RNA-seq data sets, together with quality  control and exploratory analysis reports http://imlspenticton.uzh.ch:3838/conquer/
* The marked improvements in massive parallel sequencing coupled with  single-cell sample preparations and data deconvolution have allowed  single-cell RNA sequencing (scRNA-seq) to become a powerful approach to  characterize the gene expression profile in single cells ([*1*](https://www.science.org/doi/10.1126/sciadv.abh2169#pill-R1), [*2*](https://www.science.org/doi/10.1126/sciadv.abh2169#pill-R2)). The objective of the international collaborative effort **Human Cell Atlas** ([www.humancellatlas.org](http://www.humancellatlas.org)) takes advantage of this new technology platform to study the  distinctive gene expression profiles on RNA level across diverse cell  and tissue types and connect this information with classical cellular  descriptions, such as location and morphology ([*3*](https://www.science.org/doi/10.1126/sciadv.abh2169#pill-R3)). 

*  an effort to combine the information from these two efforts to create a publicly available HPA **Single Cell Type Atlas** with genome-wide  expression data from scRNA-seq experiments integrated with the spatial  antibody-based bioimaging data.  a Single Cell Type Atlas has been 
  launched (www.proteinatlas.org/celltype)  A single–cell type transcriptomics map of  
  human tissues
* single cell data breast cancer; scRNA-seq data are available for in-browser exploration and download through the **Broad Institute Single Cell portal** at https://singlecell. broadinstitute.org/single_cell/study/SCP1039. Processed scRNA-seq data from this study are also available through the Gene Expression Omnibus under accession number GSE176078.

##### Other

* **Human Developmental Cell Atlas** (HDCA) initiative, which is part of  the Human Cell Atlas, aims to create a comprehensive reference map of  cells during development

* an **RNAi proliferation screen** in which a genome-wide shRNA library (~16,000 genes) had been used to identify essential genes (An integrated genomics approach identifies drivers of proliferation in luminal-subtype human breast cancer) Marcotte et al. 2012

* essential genes from the online **gene essentiality data- base** (OGEE) and the **Database of Essential genes** (DEG)

* Gene Active Ranking Profile (GARP)-normalized data were obtained from the **COLT** database http://colt.ccbr.utoronto.ca/cancer Genome-wide pooled shRNA  screens enable global identification of the genes essential for cancer  cell survival and proliferation and provide a ‘functional genetic’ map  of human cancer to complement genomic studies. Using a lentiviral shRNA  library targeting approximately 16 000 human genes and a newly developed scoring approach, we identified essential gene profiles in more than 70 breast, pancreatic and ovarian cancer cell lines. We developed a  web-accessible database system for capturing information from each step  in our standardized screening pipeline and a gene-centric search tool  for exploring shRNA activities within a given cell line or across  multiple cell lines. The database consists of a laboratory information  and management system for tracking each step of a pooled shRNA screen as well as a web interface for querying and visualization of shRNA and  gene-level performance across multiple cancer cell lines. COLT-Cancer  Version 1.0 is currently accessible at http://colt.ccbr.utoronto.ca/cancer.



TCIA	digital histopathology slides


Achilles shRNADEMETER knockout scores were downloaded from The Broad Institute for all cell lines in CCLE for all TFs a



gene essentiality data from: McFarland JM, Ho ZV, Kugener G, Dempster JM, Montgomery PG, Bryan JG, et al. Improved estima- tion of cancer dependencies from large-scale RNAi screens using model-based normalization and data integration. Nat Commun. 2018; 9(1):4610. Epub 2018/11/06. https://doi.org/10.1038/s41467-018- 06916-5 PubMed Central PMCID: PMC6214982. PMID: 30389920
39.



DisGeNET is one of the largest public collections of gene-disease asso- ciations.



CHASM (Carter et al., 2009) and GISTIC2.0 (Mermel et al., 2011) putative driver events versus



ISwine = an online comprehensive knowledgebase in which we incorporated almost all the published swine multi-omics data.





Software tools, databases and resources in metabolomics: updates from 2018 to 2019 O'Shea and Misra 2020





 CR2Cancer: a database for chromatin regulators in human cancer 



GDA scores provided by Dis- GeNET (https ://www.disge net.org/)



very nice review with list of ressources From correlation to causation: analysis of metabolomics data using 
systems biology approaches Rosato et al. 2018

LitPathExplorer, a visual text analytics tool that integrates advanced  text mining, semi-supervised learning and interactive visualization, to  facilitate the exploration and analysis of pathway models using  statements (i.e. events) extracted automatically from the literature and organized according to levels of confidence.



we developed the Biofactoid (biofactoid.org) software suite, which crowdsources structured knowledge in articles from authors. Biofactoid is a web-based system that lets scientists draw a network of interactions between genes, their products, and chemical compounds and employs smart-automation to translate user input into a structured language using the expressive power of a formal ontology. The



see refs in                         The status of causality in biological databases: data resources and data retrieval possibilities to support logical  modeling                    (Touré et al. 2020)