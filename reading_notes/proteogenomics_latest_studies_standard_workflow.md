### Proteogenomic and metabolomic characterization of human glioblastoma - Wang et al. 2021

* genomic data analysis	
  * copy number variant calling
    * CNVEX
  * somatic variant calling
    * somatic variants were called using Sentieon TNScope
    * ) AnnotationusingVEP
    * Additional annotations were attachedtothevariantsusingVcfAnno including allele frequencies inGNOMAD(v2.0.1),ClinVar (downloadedDec.2019),Cosmic(down- loaded Dec. 2019), dbSNP (20180418
  * germline variant calling
    * DNAScope
  * TERT promoter  mutation calling
  * structural variant calling
  * DNA methylation
  * classification of MGMT promoter DNA methylation status
  * telomere length quantification and genotyping (estimate the telomere length using WXS and WGS tumor and blood normal paired BAMs)
* RNA analysis
  * RNA quantification
  * RNA fusion detection
  * miRNA quantification
  * circular RNA prediction and quantification
* snRNA-seq analysis
* MS data
  * quantification of phosphopeptides
  * quantification of acetylated peptides
  * metabolome and lipidome

* other analyses
  * ancestry prediction using SNPs from 1000 genomes project
  * multi-omics subtyping using NMF
  * expression-based TCGA subtyping
  * unsupervised clustering of DNA methylation
  * MolecularNeuroPathology (MNP) DNA Methylation Classification of Central Nervous System (CNS) Tumors
  * unsupervised clustering of miRNA expression Unsupervised
  * unsupervised clustering of individual data type
  * determination of stemness score
  * multi-omics cis association analysis using iProFun (integrate somatic mutation, CNV, DNA methylation, RNA, protein, phosphorylation (phospho) and acetylation (acetyl) levels via iProFun (Song et al., 2019) to investigate the functional impacts of DNA alterations in GBM.)
  * mutation impact on the RNA, proteome, phosphoproteome, lipidome and metabolome
  * protein and RNA marker identification for multi-omics mixed subtype
  * kinase-substrate pairs regression analysis
  * differential proteomic, phosphosite, metabolome and lipidome analysis
  * copy number impact on transcriptome and proteome
  * cell type enrichment deconvolution using gene expression
  * immune clustering using cell type enrichment scores
  * cell type enrichment deconvolution using protein abundance
  * deep learning histopathology image analysis
  * gene set enrichment analysis
  * histone protein and acetylation calculation
  * histone acetylation association with HATs and HDACs
  * pathway enrichment analysis along histones H2B and H3/H4 acetylation axes
  * causative pathway interaction discovery using CausalPath
  * L1000 and P100 drug connectivity analysis



<u>CausalPath</u>

*Application of CausalPath (Babur et al., 2018) to the protein and phosphoprotein expression data (Figure S7A, Table S5) showed upregulation of the hypoxia pathway in mesenchymal tumors, evidenced by significant activation of multiple HIF-1 downstream targets (network permutation p = 0.0012).* 

To discover the causative pathway interactions in our proteomic and phosphoproteomic data, we took the **normalized expression of protein** with < 10% missing values **and phosphoprotein** with < 25% missing values across all tumor and normal samples as the input to CausalPath (commit 7c5b934). We ran CausalPath in the mode that **tests the mean values between test and control groups** (value-transformation = significant-change-of-mean), where the test group being the tumors of one subtype and control group being the rest of the tumors. The pathway interaction discovery data source was Pathway Commons v9 (built-in- network-resource-selection = PC). Additionally, we enabled the causal reasoning if all the downstream targets of a gene were active or inactive (calculate-network-significance = true, use-network-significance-for-causal-reasoning = true, permutations- for-significance = 10000). The causative interactions with FDR < 0.05 were extracted and visualized (fdr-threshold-for-data-signif- icance = 0.05 phosphoprotein, fdr-threshold-for-data-significance = 0.05 protein, fdr-threshold-for-network-significance = 0.05). Full result tables were available in Table S5.





## A proteogenomic portrait of lung squamous cell carcinoma - Satpathy et al. 2021

* genomic data analysis
  * copy number calling
  
    * CNVEX
  
  * somatic variant calling
  
    * somatic variants were called using Sentieon TNScope
    * AnnotationusingVEP
    * Additional annotations were attachedtothevariantsusingVcfAnno including allele frequencies inGNOMAD(v2.0.1),ClinVar (downloadedDec.2019),Cosmic(down- loaded Dec. 2019), dbSNP (20180418
  
  * germline variant calling
  
    * DNAScope
  
  * to identify signif- icantly amplified or deleted focal-level and arm-level events
  
    * GISTIC
  
  * to evaluate the significance of mutated genes and es-timate mutation densities of sample
  
    *  MutSig2CV 
  
* RNA data analysis - RNAseq and miRNAseq quantification
  * RNA quantification
    * featureCounts
  * Isoform specific RNA quantificatio
    * CIRI (v2.0.6) was used to call circular RNA
    * andBWA(version 0.7.17-r1188) was used as the mapping tool.
    * RSEM (version 1.3.1) with Bowtie2 to quantify gene, isoform, and circular RNA expression based on the mixed transcripts.
  * DNA methylation data preprocessing
    * detect QC issues: Mclust
  * miRNA-seq data analysis miRNA-seq
  
* proteomics data analysis

* proteogenomic analysis
  * unsupervised multi-omic clustering using NMF (non-negative matrix factorization (NMF)-based multi-omic clustering using protein, phosphosite, acetylsite, RNA transcript and gene copy number variants (CNV) as previously described); weights determined with ssGSEA (normalized enrichment scores (NES) of cancer-relevant gene sets by projecting the matrix of signed multi-omic feature weights (Wsigned) onto Hallmark pathway gene sets (Liberzon et al., 2015) using ssGSEA );  The entire workflow described above has been implemented as a module for PANOPLY 

  * integrative analysis with Stewart et al. (describe a way to integrate their protein data with protein data from another dataset !)

  * RNA subtyping

  * chromosomal instability index (CIN) -> from the CNA data

  * fusion detection and analysis

    * using a combination of CRISP, CODAC MI-ONCOSEQ pipeline (Robinson et al., 2017; Wu et al., 2018), fusioncatcher_v1.10 (Nicorici et al.) and arriba_v1.1.0 (https://github.com/suhrig/arriba/).

  * mRNA and Protein correlation

  * CNA-driven cis and trans effects

  * CMAP analysis

    * Candidate genes driving response to copy number alterations were identified using large-scale Connectivity Map (CMAP) queries. The CMAP (Lamb et al., 2006; Subramanian et al., 2017) is a collection of about 1.3 million gene expression profiles from cell lines treated with bioactive small molecules (?20,000 drug perturbagens), shRNA gene knockdowns (?4,300) and ectopic expression of genes. The CMAP dataset is available on GEO (Series GSE92742). For this analysis, we use the Level 5 (signatures from aggregating replicates) TouchStone dataset with 473,647 total profiles, containing 36,720 gene knock-down profiles, with measurements for 12,328 genes. See https://clue.io/GEO-guide for more information. To identify candidate driver genes, proteome profiles of copy number-altered samples were correlated with gene knockdown mRNA profiles in the above CMAP dataset, and enrichment of up/downregulated genes was evaluated.

  * Differentially expressed genes between the 5 NMF LSCC subtypes were identified using the limma package

  * Protein abundance comparisons were performed between all 5NMFsubtypes using the Wilcoxon rank-sum test and p values were adjusted using the Benjamini & Hochberg method. S

  * LINCS analysis

    * These NMF-specific signatures were used as input to calculate normalized weighted connectivity scores (WTCS) against the Library of Integrated Network-Based Cellular Signatures (LINCS) L1000 perturbation-response signa- tures.

  * defining cancer-associated genes

    * compiled using the Census website https://cancer.sanger.ac.uk/census. The

  * CpG island methylator phenotype

    * gene-level methylation score

  * iProFun based cis association analysis (iProFun = an integrative analysis tool to identify multi-omic molecular quantitative traits (QT) perturbed by DNA-level variations)

    * five functional molecular quantitative traits (gene expression, protein, phosphoprotein, acetylprotein, ubiqui- tylproteome levels) for their associations with DNA methylation, accounting for mutation, copy number variation, age, gender, tumor purity and smoking status
    * Tumor purity was determined using TSNet from RNA-seq data. The

  * EMT-derived cluster and fibroblast enrichment

    * The ssGSEA was performed on a protein dataset using GSVA (v3.11) to calculate normalized enrichment score for EMT and fibroblast proliferation using MSigDB (v6.1)

  * differential marker analysis 

    * a Wilcoxon signed rank test was performed on TMT-based global proteomic data between tumor and matched normal samples to determine differential abundance of proteins between tumor and NAT samples

  * survival analyses

    * Normalized RNA expression was categorized into low and high expression group based on mean.; The hazard ratio [exp(cox coefficient)] was used to compare poor survival in low and high-expression group.
    * We used Kaplan-Meier analysis to explore survival differ- ences associated with tumor grade, mutation burden, mutation status for SMGs, CUL3-NFE2L2-KEAP1 combined mutation status, CNA for KAT6A/SOX2/TP63/FGFR1/CDKN2A, immune subtype, NMF cluster, CIN, tumor grade and ploid

  * continuous smoking score

    * Non-negative matrix factorization (NMF) was used in deciphering mutation signatures in cancer somatic mutations stratified by 96 base substitutions in tri-nucleotide sequence contexts. To obtain a reliable signature profile, we used SomaticWrapper to call mu- tations from WGS data. SignatureAnalyzer exploited the Bayesian variant of the NMF algorithm and enabled an inference for the optimal number of signatures from the data itself at a balance between the data fidelity (likelihood) and the model complexity (reg- ularization). After decomposing into three signatures, signatures were compared against known signatures derived from COSMIC (Tate et al., 2019) and cosine similarity was calculated to identify the best match.

  * Identification of differentially regulated events in NRF2 mutant tumors

    * normalized levels ofmRNA/protein/phosphoprotein (log transformeddata)werefit intoa linear regressionmodel. Inaddi- tion to NFR2 mutation status, gender, tumor purity and ethnicity were also included in the mode

  * immune cluster identification based on cell type composition

    * abundance of 64 different cell types were computed via xCell based on transcriptomic profiles (Aran
    * Consensus clustering was performed based on immune cells, fibroblasts, endothelial and epithelial cells from xCell using the R package ConsensusClusterPlus (Wilkerson
    * signatures were partitioned into three major clusters using the Partitioning Around Medoids (PAM) algorithm

  * TCGA pan-cancer immune subtyping

    * Gene expression data (log2 FPKM; lscc-v3.0-rnaseq-uq-fpkm-log2-NArm) was input to the ImmuneSubtypeClassifier R package

  * ranking tumors by inferred activity of IFN-g pathway

    * true biological activity of a pathway is regulated by collective changes of expression levels of majority of proteins
      involved in this pathway.

    * difference of a pathway activity between tumors can be assessed by a difference in positioning of
      expression levels of proteins involved in this pathway in ranked list of expression levels of all proteins in each of tumors. 

    * relative positioning of pathway proteins between tumor by determining two probabilities: 

      1. a probability of pathway
         proteins to occupy by random the observed positions in a list of tumor proteins ranked by expression levels from the top to the bottom
      2. a probability to occupy by random the observed positions in a list of expression levels ranked from the bottom to the
         top (Reva et al., 2020). 

    * the inferred relative activation of a given pathway across tumors assessed as a negative logarithm of
      the ratio of the above ‘‘top’’ and ‘‘bottom’’ probabilities. 

      => Thus, for a pathway of a single protein, its relative activity across tumors was
      assessed as a negative log of ratio of two numbers: 

      1. a number of proteins with expression level bigger than an expression level of
         given protein, 
      2.  a number of proteins with expression levels less than an expression level of given protein.

    *  For pathways of multiple
      proteins, the ‘‘top’’ and ‘‘bottom’’ probabilities were computed as geometrically averaged P values computed for each of proteins
      using Fisher’s exact test, given protein’s ranks in a list of pathway proteins and in a list of ranked proteins of a tumor, a number of
      proteins in a pathway, and the total number of proteins with the assessed expression level in a given tumor (Reva et al., 2020).

    *  the scoring function is positive, when expression
      levels of pathway’s proteins are overrepresented among top expressed proteins of a tumor, and it is negative, when pathway’s pro-teins are at the bottom of expressed proteins of a tumor; the scoring function is close to zero, when expression levels are distributed
      by random or equally shifted toward top or bottom. (thermodynamic interpretation)

  * estimation of tumor purity, stromal and immune scores 

    * Besides xCell, we utilized ESTIMATE (Yoshihara et al., 2013) to infer immune and stromal scores based on gene expression data 
    * Cibersort absolute immune scores were obtained by evaluating upper-quartile normalized RNA-seq FPKM data
    * To infer tumor purity, TSNet was uti- lized (Petralia et al., 2018)(Table S6). 

  * differentially expressed genes and pathway analysis

    * For each immune cluster, considering the set of genes up(down)-regulated with Benjamini-Hochberg adjusted p value lower than 10%, a Fisher’s exact test was implemented to derive enriched pathways.
    * ssGSEA (Barbie et al., 2009) was utilized to obtain pathway scores based on RNA-seq and global prote- omics data using the R package GSVA (Ha¨nzelmann et al., 2013). For this analysis, pathways from the Reactome (Fabregat et al., 2018), KEGG (Kanehisa et al., 2017) and Hallmark (Liberzon et al., 2015) databases were

  * deriving RTK CBPE scores

    * A total of 42 human RTKs were present in our proteomics dataset. For each phosphosite in our dataset we computed a linear association with each of the RTKs

  * independent component analysis

    * Independent component analysis
      ICA was performed with a workflow modified from previously described (Liu et al., 2019). Decomposition was run for 100 times on the
      matrix of protein abundance difference between tumor/NAT pairs (n = 99). Independent components were in the form of vectors
      comprised with weights of all genes in the original data. Components extracted from each run were clustered using HDBSCAN al-gorithm (McInnes et al., 2017) with cosine distance as dissimilarity metric, min_cluster_size
    * (correlation between the extracted signatures and known clinical characteristics were examined by regressing the corresponding mixing scores for all members of a component cluster against 64 sample annotations to obtain within-cluster average of log10 p values)
    * Signifi- cance was controlled for multiple testing at 0.01 level (log10 (p value) = ?5.3). Each signature vector (cluster centroid) was submitted to GSEA pre-ranked test for functional annotations

  * mutation-based cis- and trans-effects 

    * (we examined the cis- and trans-effects of 22 genes with somatic mutations that were significant in a previous large-scale TCGA LSCC study (Bailey et al., 2018) on the RNA, proteome, and phosphoproteome of known interactome DBs including Omnipath, Phosphositeplus, DEPOD, Signor, and CORUM)

  * germline quantitative trait loci (QTL) analysis 

    * (to identify germline genetic variants that explain variation in tumor gene (eQTL) and protein (pQTL) expression, we utilized the gold- standard mapping pipeline https://github.com/molgenis/systemsgenetics/wiki/eQTL-mapping-analysis-cookbook-(eQTLGen). 
    * To adjust for population stratification, we identified the multi-dimensional scaling (MDS) components using the genotype data using PLINK v1.07

  * miRNA analysis presented

    * targets of miRNAs were downloaded from the miRNA targets database miRTarBase and only the miRNA/target pairs with strong experimental evidence were retained

  * pathway projection using ssGSEA

    * e single sample Gene Set Enrichment Analysis (ssGSEA) implementation available on https://github.com/broadinstitute/ssGSEA2.0 was used to project log2(FPKM) mRNA abundances to MSigDB cancer hallmark gene sets using the following parameters

  * phosphorylation-driven signature analysis 

    * We performed phosphosite-specific signature enrichment analysis (PTM-SEA) (Krug et al., 2019) to identify dysregulated phosphor- ylation-driven pathways
    * To adequately account for both magnitude and variance of measured phosphosite abundance, we used p values derived from application of the Wilcoxon rank-sum test to phosphorylation data as ranking for PTM-SEA
    * PTM-SEA relies on site-specific annotation provided by PTMsigDB and thus a single site-centric data matrix data is required such
      that each row corresponds to a single phosphosite. We
    * the heuristic method introduced by Krug et al. (Krug et al., 2019) to deconvolute multiple phosphorylated peptides
      to separate data points (log-transformed and signed p values). Briefly, phosphosites measured on different phospho-proteoform peptides were resolved by using the p value derived from the least modified version of the peptid
    * We queried the PTM signatures database (PTMsigDB) v1.9.0 downloaded from https://github.com/broadinstitute/ssGSEA2.0/ tree/master/db/ptmsigdb using the flanking amino acid sequence (+/? 7 aa) as primary identifier. We used the implementation of PTM-SEA available on GitHub (https://github.com/broadinstitute/ssGSEA2.0) using the command interface R-script (ssgsea- cli.R). T

  * CDKN2A and RB1 annotations and pathway analysis

    * Comprehensive tumor annotation for CDKN2A and RB1 genomic status was carried out using multiple molecular features for each patient. 
    * 1) mutation types (missense/in-frame indels mutations, nonsense (stop gain rameshift indels) mutations, and splice site (splice donor, splice acceptor) mutations) as separate categories and 2) copy number data.

  * Association analysis between KGG-site abundances and E3 ligases and DUBs 

    * A list of known human E3 ubiquitin and ubiquitin-like ligases and DUBs was compiled from (Medvar et al., 2016; Nijman et al., 2005). We then fit a linear model using limma in R with the formula kgg_site_abundance ?protein_abundance, followed by empirical bayes shrinkage. The

  * cluster and pathway analysis of significantly modulated K-GG sites in tumors

    * Genes in these site-wise clusters were used for pathway enrichment analysis against the KEGG, Reactome, and WikiPathways databases using g:profiler (Raudvere et al., 2019). Pathway enrichment was performed using a gene background containing all observable genes in the K-GG dataset. Pathway

  * PTM CLUMPS analysis 

    * We employ two methods to select tumor-specific sites to include in structural analysis. 
    * first, we take PTM-sites for solely tumor-
      derived samples and binarize modifications by negative versus positive normalized signal. These are considered as individual ‘‘events’’ as done with mutations in CLUMPS v1 (Kamburov et al., 2015).
    *  For a more robust approximation of tumor-specific acet- ylation or ubiquitination, we perform differential expression using Limma between NAT and Tumor samples (; We then binarize these tumor sites as done in the first approach. We first map the PTM-sites to corresponding UNIPROT ids (ID) with available PDBs. For each crystal structure, we compute an initial WAP score and randomly sample sites as done in CLUMPS.We generate an empirical p value based on a random sampling of lysines in each crystal structure to limit the selection to residues capable of ubiquitination or acetylation. The null hypothesis we define is each random sampling of lysines will have a WAP score less than or equal to the initially computed WAP score. We run 1e6 permu- tations to generate an empirical p valu

  * DepMap genetic dependency and drug response analysis

    * Cell lines annotated as ‘‘NSCLC_squamous’’ from DepMap were considered as LSCC for which molecular profiles, dependency scores, and drug response were used.
    * To determine copy number status in LSCC cell lines, copy number segmentation files for LSCC cell lines were processed with GIS-
      TIC2.0
    * For identifying if top protein biomarkers (502 proteins significantly overexpressed (log2(FC) > 2, FDR < 0.01) in tumors
      relative to their matched NATs, most with coherent overexpression in multi-omic analysis) also conferred altered dependencies in LSCC cell lines, we leveraged DepMap genetic dependency dataset (CRISPR Avana Public 20Q3) that contained 18119 genes and 789 cell lines (https://depmap.org/portal/download/ file: Achilles_gene_effect.csv)
    * Continuous log2 copy number data from DepMap was used in Figure S4M to correlate SOX2 copy number with EZH2 shRNA de- pendency data (Combined shRNA screen from DEMETER2 Data v6 in DepMap). Drug

  * CausalPath analysis 

    * (CausalPath (Babur et al., 2018) searches for known biological mechanisms that can explain correlated proteomic changes in terms of causal hypotheses)
    * We used the OncoKB database for oncoprotein and tumor suppressor classification (excluded proteins that have both annotations), and used PhosphoSitePlus and Signor databases for the activating/inhibiting classification of phosphorylation sites. In the phosphorylation regulation networks, we included only the targetable regulators (activated proteins) and excluded the untargetable regulators (inactivated proteins).

  * variant peptide identification 

    * Weused NeoFlow (https://github.com/bzhanglab/neoflow) for neoantigen prediction
    * pecifically, Optitype (Szolek et al., 2014) was used to find human leukocyte antigens (HLA) for each sample based on WES data. Then we used netMHCpan (Jurtz et al., 2017) to predict HLA peptide binding affinity for somatic mutation–derived variant peptides with a length between 8-11 amino acids. T
    * . We used Customprodbj (Wen et al., 2020)(https://github.com/bzhanglab/customprodbj) for customized database construction.

  * cancer/testis antigen prediction (cancer/testis (CT) antigens were downloaded from the CTdatabase (Almeida et al., 2009))

    * Cancer/testis (CT) antigens were downloaded from the CTdatabase (Almeida et al., 2009).

  * LSCC, HNSCC and LUAD integrative analysis

    * data acquisition from published manuscripts
    * Differential expression analysis
    * Copy number drivers for
    * Spearman correlation was performed for these genes between CNA and RNA and between CNA and protein. Proteins were considered drivers if the correlation between both the CNA and RNA and CNA and protein were significantly positively associated (
    * To identify genes associated with the immune score, correlation between CNA, RNA, protein, and the immune score was per- formed for
    * The immune score was calculated as the z-score transformation of the ESTIMATE Immune score, which was calculated for all three cohorts as described in the method Estimation of Tumor Purity, Stromal and Immune Scores. To

  * PROGENy scores 

    * (PROGENy (Schubert et al., 2018) was used to generate activity scores for EGFR based on RNA expression data)

<u>CausalPath</u>

CausalPath (Babur et al., 2018) **searches for known biological mechanisms that can explain correlated proteomic changes in terms of causal hypotheses**. We set CausalPath parameters to **compare tumors and NATs** [NAT=normal adjacent tissue] with a paired t test, used 0.1 as FDR threshold **for proteomic change significance and network significance**, and detected **5917 potential causal relations between proteins**. We repeated the **same analysis for each NMF subtype separately** and identified 4378 (basal-inclusive), 5334 (classical), 3048 (EMT-enriched), 3744 (inflamed-secretory), and 4332 (proliferative-primitive) relations. We used these CausalPath network results in the preparation of Figure 7C, **identifying potential upstream regulators of oncogenic phosphoproteomic changes**. Here **an oncogenic phosphoproteomic change** can be any of the following 4 events: increase of activating phosphorylation of an oncoprotein, decrease of inactivating phosphorylation of an oncoprotein, decrease of activating phosphorylation of a tumor suppressor protein, and increase of inactivating phosphorylation of a tumor suppressor protein. We used the OncoKB database for oncoprotein and tumor suppressor classification (excluded proteins that have both annotations), and used PhosphoSitePlus and Signor databases for the activating/inhibiting classification of phosphorylation sites. In the phosphorylation regulation networks, we included only the targetable regulators (activated proteins) and excluded the untargetable regulators (inactivated proteins).

Figure 7: For identifying if top protein biomarkers (502 proteins significantly overexpressed (log2(FC) > 2, FDR < 0.01) in tumors
relative to their matched NATs, most with coherent overexpression in multi-omic analysis) also conferred altered dependencies in LSCC cell lines, we leveraged DepMap genetic dependency dataset (CRISPR Avana Public 20Q3) that contained 18119 genes and 789 cell lines (https://depmap.org/portal/download/ file: Achilles_gene_effect.csv). Only 16 cell lines were classified as LSCC (Cell Line Sample Info.csv). Median dependencies were calculated and for every gene and density plot in Figure 7C shows dependency score of the 502 genes corresponding to 502 protein biomarkers relative to other genes.

Figure 7(C): Genetic dependencies of 502 proteins (log2 FC > 2, FDR < 0.01 and NAs < 50%) in LSCC cell lines (n = 16) profiled as part of the Achilles Dependency-Map project

*[...] we leveraged the TCGA LSCC dataset (Hammerman et al., 2012) to examine the association of these 502 candidate tumor biomarkers with overall survival (OS) or disease-free survival (DFS). Expression of four of the most highly differential genes showed significant association with poor OS and another 15 genes with poor DFS (Figure 7B; Figure S7C). Furthermore, knockdown of these tumor biomarker proteins reduced fitness across 16 LSCC cell lines (https://depmap.org/), suggesting crit- ical roles in key cellular transformation and proliferation pro- cesses (Figure 7C; Table S7)*



## A proteogenomic portrait of lung squamous cell carcinoma - Satpathy et al. 2021





**Proteogenomic landscape of LSCC**

separation between tumors and NATs: genomic landscape, PCA on proteomic and PTM 

annotate the CNA with protegenomic data

detect genes with CNA that could be driver oncogene

investigate the impact of CNAs on noncognate gene products by matching patterns of these significant trans-effects (vertical stripes in Figure 1D) to perturbation profiles from the Connectivity Map (CMap)

association of CNA with clinical metadata, xCell immune score, protein level

LSCC vs. tumor: hypermethylation; CIMP clusters

‘‘cascading’’ promoter methylation cis-effects acrosscognatemRNA, protein,andPTMabundances

**Multi-omic clustering identifies five LSCC molecular subtypes, including one that is EMT-enriched**

multi-omic NMF to define subtypes

association of subtypes with cohort metadata

comparison clusters with TCGA-derived RNA clusters

subtypes association with pathway enrichment, CIMP, CNA, phosphorylation, protein level, xCell signature, immunohistochemical co-staining

a Library of Integrated Network-Based Cellular Signatures- based (LINCS) query for compounds that reversed subtype specific signature showed enrichment for TGFb inhibitors (F

integrate with prior proteome data to recapitulate the multi-omic clusters

**NMF EMT-E subtype tumors show phosphorylation- driven PDGFR and ROR2 signaling**

other modes of RTK activation, inferred from phosphoproteomic data, may

serine/threonine-predomi- nant correlation-based phosphosite enrichment score (CBPE score) for all RTKs in our tumor cohort. Of nine RTKs with high CBPE scores in LSCC tumors, seven were significantly associ- ated with NMF subtypes (Kruskal-Wallis

CBPE scores across subtypes; phosphosite- based evidence for upregulation

**Loss of CDK4/6 pathway inhibitors is a universal feature of LSCC but Rb1 expression is variable**

investigated the impact of recurrent mutations on cognate RNA, proteins, and PTMs (cis-effects) and on a set of CAGs (trans-effects)

e.g. association between mutation and RNA expression/protein/phosphoprotein levels

**NRF2 pathway activation in tumors with and without NRF2 pathway mutations**

differential expression in mRNA, protein and phosphoprotein levels in NRF2 pathway mutated tumors

ssGSEA-based NRF2 pathway score using the proteogenomic signature

compare with phosphorylation level and mutations

**Proteogenomic analysis of chromosome 3 prioritizes therapeutic targets in LSCC**

identify copy number alterations

CNA correlation with RNA and protein expression (association copy number - transcript - protein abundance - histology)

effect of inhibitors in cell lines

association gene expression - DNA methylation

association gene expression between other genes gene expression or at protein level

abundance of a biomarker and its regulator by subtype

find druggable targets

**Crosstalk between lysine acetylation and ubiquitylation impacts cancer metabolism**

Consensus clustering of K-GG peptide abundances

detect differentially modified proteins between clusters

To identify candidate enzymes driving Ub and UbL modifica-
tions in LSCC, we correlated E3 ubiquitin ligases or deubiquity- lases (DUBs) to KGG-sites

consistency in change in protein abundance

investigate crosstalk between
lysine PTMs. We employed a modified version of CLUMPS (Kamburov et al., 2015) to detect clustering of either acetylation or ubiquitylation sites within protein 3D structures

change in protein and ubiquitylation within subtypes

acetylation sites with differential expression in tumors relative to NATs

**Immune landscape and regulation in LSCC**

Consensus clustering based on the xCell signatures 

for the clusters: enrichment analysis (expression, acetylproteome, immunohistochemical staining, CBPE scores)

pQTLs in proteins differentially abundant between clusters

**Proteomic biomarker candidates for prognosis, diagnosis, and treatment**

different protein expression and pathway enrichment in tumor vs. NAT

detection of neoantigens

assess effect of tumor biomarker knockdown in cell lines (https://depmap.org/)

comparison immune score with other type of cancers

identify tumor-specific phosphorylation events compared to other types of cancers

## Proteogenomic and metabolomic characterization of human glioblastoma - Wang et al. 2021

**Proteogenomic and metabolomic features delineate molecular subtypes of glioblastoma**

GBM samples with unmatched normal GTEX samples

compare genomic properties with TCGA

multi-omic clustering

pathway enrichment by subtypes

subtype association with clinical data

subtype association with omics traits

**Driver genetic alterations influence oncogenic protein abundance and phosphorylation**

association genetic alteration with omics traits and signaling

completed with immunohistochemistry staining patterns

**RTK signaling cascades are activated in GBM**

compy number alteration matching other omics

pathway activation (phosphorylation levels) consistency across omics

kinase-substrate study

**Distinct immune marker expression and epigenetic events characterize GBM immune subtypes**

immune-based subtype definition with cell-type immune enrichment scores using xCell signature

association immune subtypes with GBM subtypes

association with immune gene expression

identified differentially expressed proteins (DEPs) and phosphoproteins (DEPPs) in known immune targets

identify overrepresented pathways among DEPs and DEPPs

analyzed morphologic differences between immune subtypes by applying a deep learning model using sampled tumor tiles from H&E-stained sections.

**Mesenchymal tumor and microenvironment characteristics**

CausalPath on protein and phosphoprotein expression data to detect upregulated pathways in subtypes

snRNA-seq data enabled identifying mesen- chymal features in tumor and immune cells

subtype specific differences in marker expression  (snRNA-seq + bulk RNA-seq)

**Differential acetylation of histone proteins is associated with specific subtypes and pathways**

histone acetylation

clustering: detect tumors with differentially acetylated histones

upregulation histone acetylation in cancer vs. normal

Lasso linear regression between histone acetylation sites and the protein and acetylation abundances of histone acetyltransferases (HATs), bromodomain-containing proteins (BRDs), and deacetylases (HDACs). It revealed potential connec- tions between HATs and BRDs and H2B acetylation sites,

association histone acetylation and pathways

**Lipid composition and metabolomic features associated with GBM subtypes**

lipid differentially abundant across tumor subtypes

connection between selected lipid abundance and metabolomically related proteins (by subtypes and vs. normal)

association with histone acetylation and phosphorylation

metabolite abundance comparison in IDH-mutant vs IDH-wt

**Key oncogenic pathways and therapeutic opportunities in GBM**

integrate genetic alterations and the RNA, protein, and phosphosite levels per expression subtype to examine three important oncogenic signaling pathways i

n conjunction with druggability information from DGIdb (Cotto
et al., 2018) and DEPO (Sun et al., 2018), we conducted kinase- substrate and outlier analyses to identify druggable pairs

phosphorylation level and signaling pathway

similarity between alteration-specific RNA or phosphoprotein signatures with corresponding transcriptional (L1000 assay) and phosphoproteomics LINCS signatures







NB: PROGENY

Aberrant cell signaling can cause cancer and other diseases and is a focal point of drug research. A common approach is to infer signaling activity of pathways from gene expression. However, mapping gene expression to pathway components disregards the effect of post-translational modifications, and downstream signatures represent very specific experimental conditions. Here we present PROGENy, a method that overcomes both limitations by leveraging a large compendium of publicly available perturbation experiments to yield a common core of Pathway RespOnsive GENes. Unlike pathway mapping methods, PROGENy can (i) recover the effect of known driver mutations, (ii) provide or improve strong markers for drug indications, and (iii) distinguish between oncogenic and tumor suppressor pathways for patient survival. Collectively, these results show that PROGENy accurately infers pathway activity from gene expression in a wide range of conditions.



NB: iProFun - Song et al. 2019

 we propose iProFun, an integrative analysis tool to screen for  proteogenomic functional traits perturbed by DNA copy number alterations (CNAs) and DNA methylations. The goal is to characterize functional  consequences of DNA copy number and methylation alterations in tumors  and to facilitate screening for cancer drivers contributing to tumor  initiation and progression. Specifically, we consider three functional  molecular quantitative traits: mRNA expression levels, global protein  abundances, and phosphoprotein abundances. We aim to identify those  genes whose CNAs and/or DNA methylations have cis-associations with  either some or all three types of molecular traits. 

 integrative analysis pipeline—iProFun—for revealing dynamic cis  regulatory patterns in tumors. Briefly, iProFun takes as input the  association summary statistics from associating CNAs and methylations of genes to each type of cis-molecular trait, aiming to detect the joint  associations of DNA variations and molecular traits in various  association patterns. Of particular interest are the genes with  “cascading effects” on all cis molecular traits of interest and the  genes whose functional regulations are unique at global/phospho protein  levels. iProFun can incorporate prior biological knowledge through a  filtering procedure and can identify significant genes with calculated  posterior probabilities exceeding a threshold while assessing the  empirical FDR (eFDR) through permutation. Downstream enrichment analyses were also embedded into our pipeline to allow for more direct  interpretations of different association patterns.

INPUT: molecular quantitative traits (e.g. mRNA, phosphoprotein measurements) + DNA alterations (e.g. CNA, methylation) + covariate (e.g. 1st PC for population stratification) 

=> obtain association summary statistics via gene-level multiple regression => association with pattern recognition  => FDR assessment => enrichment analysis