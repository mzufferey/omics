* keywords
  * master regulator
  * gene regulatory network
  * causal network analysis
* popular tools: SCENIC, ARACne, VIPER, MARINa, DISCERN, GENIE3, Hotnet2, MOMA; lasso methods
* knwoledge databases: TF databases, PPI databases, pathways/GOs/signature 
* most current workflow : Aracne -> VIPER -> GO analysis

* knowledge integration
  * bayesian method: prior on gene-gene weights
  * NN: architecture of NN matching gene pathways
* network inference
* differential network analysis (e.g. KEDDY, Co-expression Differential Net- work Analysis (CoDiNA))
* differential rewiring analysis
* data integration: transcriptomics, proteomics, mutations
* AI-based: NN (lung example with prior knowledge architecture); VAE (GLUE: example with multi-omics data integration; Cao and Gao 2021)
* single-cell development
* issues:
  * uncomplete knowledge databases
* most interesting Bio examples: Paull et al. 2021;Walsh et al. 2017; Matsubara et al. 2019; Tapia-Carrillo et al. 2019; Tovar et al. 2015)

---



### Zhang et al. 2011

a linear model to explain similarities between disease phenotypes using gene proximities that are quantified by diffusion kernels of one or more PPI networks. We solve this model via a **Bayesian approach**, and we derive an analytic form for Bayes factor that naturally measures the strength of association between a query disease and a candidate gene and thus can be used as a score to prioritize candidate genes

This method is intrinsically capable of integrating multiple PPI networks

The **Bayesian regression approach** can achieve much higher performance than the existing CIPHER approach and the ordinary linear regression method.

Most existing computational methods for inferring causative genes from candidates are formulated as a one-class novelty learning problem that is usually solved with the guilt-by-association principle,  which suggests to compute a score from functional genomics data to quantify the strength of association between a query dis- ease and a candidate gene, and then rank candidate genes according to their scores to facilitate the selection of susceptibility genes [4].

Pearson’s correlation coeffi-
cient of similarities between phenotypes and closeness of genes in a single protein-protein interaction (PPI) network can be used as a concordance score to facilitate the priori- tization of candidate genes

However, PPI networks are far from complete.

we propose a **Bayesian**
**regression approach that can be used with either a single PPI network or multiple networks to prioritize candidate genes**. We adopt a linear model to explain disease simi- larity using gene proximity, and we solve this model via a Bayesian approach, which yields an analytic form of Bayes factor for measuring the strength of association between a query disease and a candidate gene. We then use Bayes factors as scores to prioritize candidate genes. We show the validity of assumptions of this approach, and we demonstrate the effectiveness of this approach on five PPI networks via large scale leave-one-out cross- validation experiments and comprehensive statistical analysis. We further show the capability of our approach in integrating multiple PPI networks

known associations between disease phenotypes and genes extracted from the Online Mendelian Inheritance in Man (OMIM) database

he Human Protein Reference Database (HPRD) contains human protein-protein interactions that are manually extracted from the literature by expert biolo- gists [28].

Biological General Repository for Interaction Data- sets (BioGRID) contains protein and genetic interactions of major model organism specie

the Biomolecular Interaction Network Database (BIND) contains both high-throughput and manually curated interactions between biological molecules [30]. 

the IntAct molecular interaction database (IntAct) contains protein-protein interaction derived from literature [

he Molecular INTeraction database (MINT) contains information about physical interactions between pro- teins [32]. From

We adopt a linear regression model to explain disease similarities in the phenotype similarity profile using gene similarities in one or more gene proximity profiles, and we solve this regression model via a Bayesian approach [34].



### Neapolitan et al. (2013)

the **expression level of genes in signaling pathways can be modeled using a causal Bayesian network (BN)** that is altered in tumorous tissue. 

the expression levels of genes that code for proteins on a signal transduction network (STP) are causally related and that this causal structure is altered when the STP is involved in cancer.



. We identified all the genes related to each of the 22 pathways and developed separate gene expression datasets for each pathway. We obtained significant results indicating that the causal structure of the expression levels of genes coding for proteins on STPs, which are believed to be implicated in both breast cancer and in all cancers, is more altered in the cases relative to the controls than the causal structure of the randomly chosen pathways.

To address the second issue (the prediction of how stimu-
lations and inhibitions will affect the overall activity of the STP), we need a causal model of the variables related to an STP

STPs can be modeled as causal Bayesian networks (BNs) if each node in the network represents the phosphorylation activity of a protein. A strength of BNs is that they represent probabilis- tic relationships, and therefore they can manage the noise in biological data. A second strength is that they can model the natural causal relationships in biology

the central
hypothesis to be investigated in this paper is that the expres- sion levels of genes that code for proteins on a given STP are causally related, and that this causal structure is altered when the STP is involved in a particular cancer. If this hypothesis is correct, using the ample gene expression datasets and BN learning algorithms, we can learn the causal network struc- ture of the gene expression levels in an STP that is altered in a given cancer, and then identify driver genes based on the topology of the network.

A BN consists of a directed acyclic graph (DAG) G = (V, E) whose nodeset V contains random variables and whose edges E represent relationships among the random variables, the prior prob- ability distribution of every root variable in the DAG and the conditional probability distribution of every non-root variable given each set of values of its parents. Often the DAG is a causal DAG, which is a DAG containing the edge X → Y only if X is a direct cause of Y.29 The probability dis- tribution of the variables in a BN must satisfy the Markov condition, which states that each variable in the network is probabilistically independent of its non-descendents condi- tional on its parents.

In the constraint-based approach,30 we learn a DAG
model from the conditional independencies that the data sug- gest are present in the generative probability distribution P. In the score-based approach, we assign a score to a DAG based on how well the DAG fits the data.

The Bayesian score is the probability of the Data given the DAG model.31 A popular variant of this score is the Bayesian Dirichlet equivalent uni- form (BDeu) score.32 If

Kyoto Encyclopedia of Genes and Genomes (KEGG)43 has a collection of manually drawn pathways representing our knowledge of about 136 pathways. STPs are not thought to be stand-alone networks, but rather they have inter-pathway communication.4



If we represent the phosphorylation level of each pro-
tein in an STP by a random variable and draw an arc from X to Y if there is an edge from protein X to protein Y in the STP, then we are modeling the STP as a BN. For this BN to represent the joint probability distribution of the random variables, the Markov condition must be satisfied. Woolf et al.22
argue that the steady-state concentrations
should satisfy this condition



to investigate whether the causal structure of the expression levels of genes on an STP is altered when the STP is involved in cancer, we compared results obtained using the breast cancer data for 5 STPs implicated in breast cancer, 7 STPs implicated in other cancers, and 10 random pathways. 









### Tovar et al. 2015

880 microar-
ray expression profiles from several experimental datasets that are available on the Gene Expression Omnibus site (http://www. ncbi.nlm.nih.gov/geo/GEO)



First, we **generated a network for every known human TF in the primary breast cancer gene expression dataset by using the Algorithm for the Reconstruction of Accurate Cellu- lar Networks (ARACNE)** (Basso et al., 2005; Margolin et al., 2006). ARACNE is a computational algorithm widely used to identify statistical relationships among genes, by calculating the mutual information (MI) between gene pairs from microarray expression data (Basso et al., 2005; Margolin et al., 2004)



**MARINais designed to infer transcription factors that control the**
**transition between two phenotypes A and B, as well as the mainte- by the activation or repression of specific TFs then their targets nance of the latter phenotype.** If the A→B transition is supported should be among the most differentially expressed genes between the two cellular phenotypes, with activated and repressed targets at opposite ends of the expression range. **MARINa estimates the importance and biological relevance of a TF on a given phenotype by computing the statistical significance of the overlap between its regulon and the gene expression signature using sample per- mutation to estimate the distribution of the enrichment score (ES) (Subramanian et al., 2005) in the null condition (Lefebvre et al. 2010)**



An important step for this algorithm is the selection of the Transcription Factors, since they will determine the rest of the cal- culation. A proper annotation of transcription factors is crucial for an accurate description of the process under investigation. Here, we used the HGU133A annotation file, in which we found 1142 TFs (Supplementary Material 1). This list was compared with other three lists. Those lists are available in Shimoni and Alvarez (2013), Vaquerizas et al. (2009) and http://www.bioguo.org/AnimalTFDB/ Animal Transcription Factor DataBase, respectively. We want to stress that all four lists show consistency among them



shadowingand synergy ofTMRsover their target gene sets (Lefebvre et al., 2010). 



**Causal network (CN) analysis**

(the Ingenuity Knowledge Base
(IKB)). IKB reports a series of experimentally observed cause- effect relationships related to transcription, expression, activation, molecular modification, binding events and transport processes



Network informa- tion was supplied with differential expression analysis between tumoral and healthy samples, which defined the input for the CN study. Differentially



Causal network (CN) analysis was performed with the Ingenuity Pathway Analysis method ((IPA®,



### Chen et al. 2015

Modeling complex diseases with pathology involving interacting tissues is difficult when using reverse-engineered regulatory networks due to the presence of different networks in a biopsy sample containing multiple cell types. Chen et al. introduce an approach for deconvolving regulatory modules that originate uniquely in one tissue and

We use context-specific regulatory networks to deconvolve and identify skin-specific regulatory modules with

master regu- lators (MRs). MRs represent the minimal number of transcription factors (TFs) that are predicted to specifically activate or repress a target module and, by extension, the associated physiological behavior

The inference of MRs is made possible through the reverse engineering of context-specific regulatory networks us- ing computational algorithms such as ARACNe (Margolin et al., 2006a).

However, this type of computational approach is only starting
to be implemented to target pathogenic, non-cell autonomous interactions between different tissues such

inferring MRs cannot be done directly using typical ARACNe-based analysis because of fundamental as- sumptions made during the generation of a regulatory network: (1) that the samples used are relatively pure or represent the one underlying transcriptional network; and (2) the underlying molecular behavior of a data set exists at a steady state such that each sample can be treated as a ‘‘snapshot’’ of regulatory dependency within the overall network

A contaminated sample, particularly by a tissue that exhibits a different context-specific regulatory network, can impair the accuracy of regulatory predictions. Further, when pathogenesis is depen- dent on the interaction between the two tissues, there will always be an artifact correlation between contaminant gene signatures and the molecular modules that recruit them, but are expressed in the other tissue.

This makes it difficult to clearly define modules exclusive to one tissue or the other when analyzing gene expression data generated from a mixture of the two tissues.

we leverage context-specific regulatory networks for the regulatory deconvolution of a mixed- signature gene expression profile of

to develop a framework capable of separating mixed AA tissue biopsy gene expression data into skin-specific mod- ules of AA pathology and infiltrate recruitment. We

First, we created a molecular signature comparing AA patients to controls to generate a molecular representation of AA

We created an overall gene expression signature by com-
paring patients of two distinct clinical presentations, patchy AA (AAP) and totalis and universalis (AT/AU) all against unaffected controls.

generate a context-specific transcriptional interaction network for scalp skin, we employed the ARACNe algorithm (Margolin

MR Analysis MRs for a specific gene expression signature were defined as TFs whose direct ARACNe-predicted T (regulon) are statistically enriched in the gene expression signature. Each TF’s regulon was tested for enrichment of the AAGS using Fisher’s exact test, FDR = 0.05. This analysis allows for the ranking and determination of the minimum number of TFs required to specif- ically cover a gene expression module associated with a physiological trait. http://wiki.c2b2.columbia.edu/califanolab/index.php/Software/MARINA



### Mertins et al. 2016

Here we describe quantitative mass-spectrometry-based proteomic and phosphoproteomic analyses of 105 genomically annotated breast cancers, of which 77 provided high-quality data

To provide greater analytical breadth, the NCI Clinical Proteomic Tumor Analysis Consortium (CPTAC) is using mass spectrometry to analyse the proteomes of genome-annotated TCGA tumour samples5,

The cohort included a balanced rep- resentation of PAM50-defined intrinsic subtypes9 including 25 basal- like, 29 luminal A, 33 luminal B, and 18 HER2 (ERBB2)-enriched tumours, along with 3 normal breast tissue samples. Samples were analysed by high-resolution accurate-mass tandem mass spectrom- etry (MS/MS) that included extensive peptide fractionation and phosphopeptide enrichment (Extended Data Fig. 1a)

12,553 proteins (10,062 genes) and 33,239 phosphosites, with their relative abundances quantified across tumours, were used in subsequent analyses in this study.

clustering analyses were first restricted to the reduced set of PAM50 genes. When RNA data for the 50 PAM50 genes were clustered directly (without using a classifier), the clus- tering was similar to the TCGA PAM50 annotation (second anno- tation bar in Fig. 3a). Restricting both the RNA and proteome data to the set of 35 PAM50 genes observed in the proteome produced a similar result (bottom two annotation bars in Fig. 3a), and all of the major PAM50 groups were recapitulated in the proteome almost as well as in the RNA data

Global proteome and phosphoproteome data were then used to identify proteome subtypes in an unsupervised manner. Consensus clustering identified basal-enriched, luminal-enriched, and stromal-enriched clusters (Extended Data Figs 6a–d, 7a). Unlike the clustering observed with PAM50 genes, mRNA-defined HER2- enriched tumours were distributed across these three proteomic sub- groups. 

The basal-enriched and luminal-enriched groups showed a strong overlap with the mRNA-based PAM50 basal-like and luminal subgroups, whereas stromal-enriched proteome subtype represented a mix of all PAM50 mRNA-based subtypes, and has a significantly enriched stromal signature (Extended Data Fig. 3e). Among the stromal-enriched tumours there was strong representation of reac- tive type I tumours, as classified by RPPA (Supplementary Table 12), showing agreement between the RPPA and mass-spectrometry-based protein analyses for the detection of a tumour subgroup characterized by stromal gene expression1

As the basal- and luminal-enriched proteome subgroups are
coherent, pathway analyses were conducted on these two subtypes, using the stromal-enriched subgroup as a control to assess spec- ificity (Fig. 3c, Extended Data Fig. 7b, Supplementary Table 13). The luminal-enriched subgroup was exclusively enriched for oestradiol- and ESR1-driven gene sets. In contrast, multiple gene sets were enriched and upregulated specifically in the basal-like tumours. Particularly extensive basal-like enrichment was seen for MYC target genes; for cell cycle, checkpoint, and DNA repair pathways including regulators AURKA/B, ATM, ATR, CHEK1/2, and BRCA1/2; and for immune response/inflammation, including T-cell, B-cell, and neutro- phil signatures. The complementarity of transcriptional, proteomic, and phosphoproteomic data was also highlighted in these analyses (Extended Data Fig. 7c, d)



Using phosphorylation status as a proxy for activity, phosphop-
roteome profiling can theoretically be used to develop a signalling- pathway-based cancer classification. K-means consensus clustering was therefore performed on pathways derived from single sample gene set enrichment analysis (GSEA) of phosphopeptide data (Methods, Supplementary Tables 14 and 15). Of four robustly segregated groups, subgroups 2 and 3 substantially recapitulated the stromal- and luminal-enriched proteomic subgroups, respectively (Fig. 3d, Extended Data Fig. 8a). Subgroup 4 included a majority of tumours from the basal-enriched proteomic subgroup, but was admixed par- ticularly with luminal-enriched samples. This subgroup was defined by high levels of cell cycle and checkpoint activity. All basal and a majority of non-basal samples in this subgroup had TP53 mutations. Consistent with high levels of cell cycle activity, a multivariate kinase– phosphosite abundance regression analysis highlighted CDK1 as one of the most highly connected kinases in this study (Extended Data Fig. 8b, Supplementary Table 16). Subgroup 1 was a novel subgroup defined exclusively in the phosphoproteome pathway activity domain, with no enrichment for either proteomic or PAM50 subtypes. It was defined by G protein, G-protein-coupled receptor, and inositol phos- phate metabolism signatures, as well as ionotropic glutamate signal- ling (Fig. 3d). Co-expression patterns among genes/proteins across different subgroups were also analysed using a Joint Random Forest method25 that identified network modules, such as an MMP9 mod- ule, with different interaction patterns between basal-enriched and luminal-enriched subgroups. These latter patterns appeared specific to the proteome-level data (Extended Data Fig. 8c–f, Supplementary Table 17 and Supplementary Methods)



A central goal in breast cancer research has been the identification
of druggable kinases beyond HER2







### Alvarez et al. 2016

We have previously shown that regulon analysis, using the master
regulator inference algorithm (MARINa), can help identify aberrantly activated tumor drivers17–21. However, this requires multiple samples representing the same tumor phenotype and cannot be used to assess aberrant protein activity from individual samples. To



VIPER infers protein activity by systematically analyzing expression of the protein’s regulon, which is strongly tumor-context-dependent20 (Fig. 1b). We used the algorithm for the reconstruction of accurate cellular networks (ARACNe23; Online Methods) to systematically infer regulons from tissue-specific gene expression data (

We extended ARACNe to detect maximum information path targets (Online Methods), as originally proposed in ref. 21, to allow identification of regulons that report on the activity of proteins representing indirect regu- lators of transcriptional target expression, such as signaling proteins

VIPER is based on a probabilistic framework that directly
integrates target ‘mode of regulation’, that is, whether targets are activated or repressed (Fig. 1b and Supplementary Figs. 1 and 2), statistical confidence in regulator-target interactions (Fig. 1b) and target overlap between different regulators (pleiotropy) (Fig. 1d) to compute the enrichment of a protein’s regulon in differentially expressed genes (Online Methods). Several

VIPER uses a fully probabilistic yet efficient enrichment analysis framework, supporting seamless integration of genes with different likelihoods of representing activated, repressed or undetermined targets, and probabilistic weighting of low vs. high-likelihood protein targets. To achieve this, we introduce analytic rank-based enrichment analysis (aREA) a statistical analysis based on the mean of ranks (Fig. 1c and Online Methods). Differential protein activity is quantitatively inferred as the normalized enrichment score computed by aREA.



### Aibar et al. 2017

**scenic, a computational method for simultaneous gene regulatory network reconstruction and cell-state identification from single-cell rna-seq data** (http://scenic. aertslab.org). on a compendium of single-cell data from tumors and brain, we demonstrate that cis-regulatory analysis can be exploited to guide the identification of transcription factors and cell states. scenic provides critical biological insights into the mechanisms driving cellular heterogeneity.



SCENIC workflow consists of three steps (Fig. 1a, Supplementary Fig. 1 and see Online Methods). In the first step, sets of genes that are coexpressed with TFs are identified using GENIE3 (ref. 8) (Supplementary Fig. 1a). Since the GENIE3 modules are only based on coexpression, they may include many false positives and indirect targets. To identify putative direct-binding targets, each coexpression module is subjected to cis-regulatory motif analysis using RcisTarget (Supplementary Fig. 1b and see Online Methods). Only modules with significant motif enrichment of the correct upstream regulator are retained, and they are pruned to remove indirect targets lacking motif support. We refer to these processed modules as regulons.

we developed the AUCell algorithm to score
the activity of each regulon in each cell (Supplementary Figs. 1c and 2, and see Online Methods). For a given regulon, comparing AUCell scores across cells makes it possible to identify which cells have significantly higher subnetwork activity. The resulting binary activity matrix has reduced dimensionality, which can be useful for downstream analyses.

**SCENIC workflow. SCENIC is a workflow based on three new R/bioconductor packages: (i) GENIE3, to identify potential TF targets based on coexpression; (ii) RcisTarget, to perform the TF- motif enrichment analysis and identify the direct targets (regu- lons); and (iii) AUCell, to score the activity of regulons (or other gene sets) on single cells.** We

**GENIE3. GENIE3 (ref. 8) is a method for inferring gene regula- tory networks from gene expression data. In brief, it trains ran- dom forest models predicting the expression of each gene in the data set and uses as input the expression of the TFs. The different models are then used to derive weights for the TFs, measuring their respective relevance for the prediction of the expression of each target gene.** The highest weights can be translated into TF-target regulatory links8. Since GENIE3 uses random-forest regression, it has the added value of allowing complex (e.g., non- linear) coexpression relationships between a TF and its candidate targets

The input to GENIE3 is an expression matrix. The preferred
expression values are gene-summarized counts (which might or might not use unique molecular identifiers, UMIs25). Other meas- urements, such as counts or transcripts per million (TPM) and FPKM/RPKM are also accepted as input. However, note that the first network-inference step is based on coexpression, and some authors recommend avoiding within-sample normalizations (i.e., TPM) for this task because they may induce artificial covaria- tion26

The output of GENIE3 is a table with the genes, the potential regulators, and their ‘importance measure’ (IM), which represents the weight that the TF (input gene) has in the prediction of the target. We

GRNBoost. GRNBoost is based on the same concept as GENIE3: inferring regulators for each target gene purely from the gene expression matrix. However, GRNBoost does so using the gra- dient-boosting machines (GBM)28 implementation from the XGBoost library29. A GBM is an ensemble learning algorithm that uses boosting30 as a strategy to combine multiple weak learners, like shallow trees, into a strong one. This contrasts with random forest, the method used by GENIE3, which uses bag- ging (bootstrap aggregation) for model averaging to improve regression accuracy. GRNBoost uses gradient-boosted stumps (regression trees of depth 1)31 as the base learner. GRNBoost’s main contribution is casting this multiple regression approach into a Map/Reduce32 framework based on Apache Spark2

RcisTarget. RcisTarget is a new R/Bioconductor implementation of the motif enrichment framework of i-cisTarget and i**Regulon RcisTarget identifies enriched TF-binding motifs and candidate transcription factors for a gene list.** In brief, RcisTarget is based on two steps. First, it selects DNA motifs that are significantly over- represented in the surroundings of the transcription start site (TSS) of the genes in the gene set. This is achieved by applying a recovery-based method on a database that contains genome-wide cross-species rankings for each motif. The motifs that are anno- tated to the corresponding TF and obtain a normalized enrich- ment score (NES) > 3.0 are retained. Next, for each motif and gene set, RcisTarget predicts candidate target genes (i.e., genes in the gene set that are ranked above the leading edge). 

AUCell. **AUCell is a new method that allows researchers to iden- tify cells with active gene regulatory networks in single-cell RNA- seq data.** The input to AUCell is a gene set, and the output is the gene set ‘activity’ in each cell. In SCENIC, these gene sets are the regulons, which consist of the TFs and their putative targets. AUCell calculates the enrichment of the regulon as an area under the recovery curve (AUC) across the ranking of all genes in a par- ticular cell, whereby genes are ranked by their expression value. This method is therefore independent of the gene expression units and the normalization procedure

Cell clustering based on gene regulatory networks. The cell regulon activity is summarized in a matrix in which the columns represent the cells and the rows the regulons. In the binary regu- lon activity matrix, the coordinates of the matrix that correspond to active regulons in a given cell will contain a “1,” and “0” oth- erwise. The equivalent matrix, which contains the continuous AUC values for each cell regulon, is normally referred to as the AUC activity matrix. Clustering of either of the regulon activ- ity matrices reveals groups of regulons (jointly, a network) that are recurrently active across a subset of cells. The binary activity matrix tends to highlight higher order similarities across cells (and therefore highly reduces batch effects and technical biases); on the other hand, the AUC matrix allows researchers to observe more subtle changes. 





### Gao et al. 2017

i**dentifying cancer pathway modules through coordination between coverage and exclusivity.**

**CovEx, to predict the specific patient oriented modules by 1) discovering candidate modules for each considered gene, 2) extracting significant candidates by harmonizing coverage and exclusivity and, 3) further selecting the patient oriented modules based on a set cover model.**

https:// sourceforge.net/projects/cancer-pathway/files/

a driver pathway module usually exhibits two combinatorial properties: high coverage and high exclusivity

High coverage means that most patients have at least one mutated gene in the module. High

High exclusivity means that most patients have only one mutated gene in that module. 

RME [16] calculated the exclusivity weight as the percentage of covered patients that contain exactly one mutation within a gene set. Another exclusive metric, Dendrix weight [17], was defined as the difference of coverage and coverage overlap of a gene set (see Methods Section for details). Based on Dendrix weight, Vandin et al. proposed a greedy algorithm and a Markov chain Monte Carlo (MCMC) algorithm [17], and Zhao et al. introduced MDPFinder including a binary linear programming (BLP) method and a genetic algorithm [18] to identify large weight modules. Compared to MCMC algorithm, the BLP method is more efficient. A multi-objective optimization model based on a Genetic Algorithm (MOGA) was introduced to adjust the trade-off between coverage and exclusivity [19]. Multi- Dendrix [20] designed a new metric as the sum of Dendrix weights and adopted a new programming model to identify multiple modules simultaneously. CoMDP proposed an exact mathematical programming method to identify co-occurring mutated driver pathways [21]. 

genes of which some are mutated frequently and rarely for others, making the exclusive metrics employed in these methods unable to handle this broad spectrum of mutational frequencies properly

Gene sets identified without consideration of pathways or PPI networks may not be correlated and thus not necessarily to form a driver module

the approaches for identification of driver genes or pathway modules over PPI networks should be developed. To

To integrate two datasets together, one fast and reliable method was presented in [23]. A combinatorial model was proposed for global module detection in complex networks [24]. HotNet2 [25] is a network based method that delves into the long tail of rarely mutated genes and finds mutated subnetworks. MEMo [26] and MEMCover [27] are both network based methods to systematically identify mutually exclusive network modules. MEMo outputs the significantly exclusive modules evaluated by a random permutation testing method. MEMCover evaluates the mutual exclusivity degree for gene pairs with random permutation testing method. Both MEMo and MEMCover only consider those gene pairs representing interactions in a PPI network, which restricts the discovered networks to existing interaction networks. 

the probabilistic methods are too complicated to efficiently calculate the exclusive scores for gene modules

A combinatorial evaluation metric overcoming the limitation of current combinatorial metrics would be desirable and may be applied to search for mutually exclusive modules to a much larger scale.

we present a network-based algorithm to identify exclusive network modules. The algorithm consists of three phases. In the first phase, we exhaustively enumerate gene sets for each considered gene in a local constructed influence network.

* we search for the candidate modules by optimizing Dendrix weight on local networks each roots at a node across an influence network which measures the topological relationships between genes in the dataset. As we previously mentioned, the candidate modules identified in the first phase may not be exclusive at all. In

In the second phase, by a new designed combinatorial coverage and exclusivity evaluation metric, CovEx, which overcomes the limitation of Dendrix weight (see Methods Section for details), we filter out those candidate modules with poor coverage and exclusivity property. After the first two phases, only significant candidate modules were left.

Biologically, each patient should have at least one driver module, and the driver module for different patients may be different. To identify the specific driver module for each patient, we employ a minimum set cover model to select the specific patient oriented driver module in the third phase. Obviously, the modules selected in the third phase may be critical and are most likely to be the desired pathway modules.



### Fu et al. 2017

d probabilistic graphical model to develop a new method that integrates genetic interaction and protein interaction data and infers exquisitely detailed pathway structure. **We modeled the pathway network as Bayesian network** and applied this model to infer pathways for the coherent subsets of the global genetic interaction profiles, and the available data set of endoplasmic reticulum genes.

The protein interaction data were derived from the BioGRID database. Our

. In order to automatically identify detailed pathway structures using high-throughput genetic interaction data, the activity pathway network (APN) was developed [26]. However, these available approaches cannot fully take advantages of the comple- mentarity between protein and genetic interaction data to infer the biological pathway structures.

we present **a Bayesian model that inte-**
**grates high-throughput protein and genetic interaction data to reconstruct detailed biological pathway structures. The model can organize related genes into the corre- sponding pathways, arrange the order of genes within each pathway, and d**ecide the orientation of each intercon- nection. Based

Based on protein interaction network, the model predicts detailed pathway structures by using genetic interaction information to delete redundancy edges and reorient the kept edges in the network.

Similar to APN [26], our model represents a biological pathway network as a Bayesian network [27], in which each node presents the activity of a gene product. Different from APN that drew network sample from complete network, our method introducing protein interaction networks as underlying pathway structures. 

a scoring func- tion is defined by gene pairwise score, which can avoid the unadjusted balance between gene pairwise score and edge score in the APN. T

our model is able to improve computational efficiency of stochastic simulation algo- rithm and overcome the limitation of APN that some edges in the results are difficult to interpret. In our model, each edge in the network can capture physical docking, and represent functional dependency



### Walsh et al. 2017

A model of transcriptional regulation in breast carcinoma is assembled with **ARACNe**

Interrogating this network reveals a transcriptional hierarchy underlying metastasis



we reverse engineered and inter- rogated a breast cancer-specific transcriptional interaction network (interactome) to define tran- scriptional control structures causally respon- sible for regulating genetic programs underlying breast cancer metastasis in individual patients.

we ap- proached metastatic progression as a transition between two cellular states defined by the differential gene expression signa- ture of these states in the same patient, using patient-matched primary and metastatic samples. 

we investigated the specific transcriptional regulators responsible for initiating this transition, and ultimately for maintaining the stability of the metastatic state, based on their mechanistic ability to regulate differentially expressed genes (i.e., the transition signature).

the specific transcriptional regulators that determine such cancer-related state transitions can be efficiently and systematically elucidated by interrogating tu- mor-specific transcriptional networks (henceforth interac- tomes) with representative differential gene expression signa- tures, **using the virtual inference of protein activity by regulon enrichment analysis (VIPER) algorithm (Alvarez et al., 2016), which further extends the master regulator inference algorithm (MARINa)** (Lefebvre et al., 2010) to the analysis of single samples

we first assembled a breast cancer-specific regulatory network, using the algorithm for the reconstruction of accurate cellular networks (ARACNe) (Basso et al., 2005; Margolin et al., 2006b). 

Then we identified differen- tial gene expression signatures representing same-patient cell state transitions from primary tumors to lymph node metas- tases, using both ER-positive (ER+) and triple-negative breast cancer (TNBC) samples. 

Finally, we **used the VIPER algorithm to prioritize transcriptional regulators that are the most likely causal determinants of these metastasis-related signatures**



we first characterized the transcriptional signature representative of breast carcinoma metastatic pro- gression (metastasis [MET] gene expression signature [MET- GES]). This was achieved by differential expression analysis of patient-matched primary tumors and lymph node metastases from 20 ER+ and 11 TNBC patients

To identify the genes that causally implement the MET-GES,
thus representing candidate causal determinants of breast carci- noma metastatic progression, we utilized the VIPER algorithm

To assemble a breast carcinoma-specific interactome, we analyzed 851 The Cancer Genome Atlas (TCGA) breast carcinoma gene expression pro- files using the ARACNe algorithm

**ARACNe is an information theory-based approach to infer mechanistic interactions between transcription factors (TFs) and target genes based on large sets of gene expression data, which has proven very effective in assembling interac- tomes for VIPER analysis** (Alvarez

The **ARACNe-inferred breast cancer interactome included 1,748 TFs, regulating 18,783 target genes through 365,634 transcriptional interactions.**

**Finally, we inferred the regulatory proteins that are candi-**
**date drivers of the MET-GES by VIPER analysis of the breast carcinoma interactome.** The

The algorithm prioritizes the regulatory proteins that are the most likely determinants of an observed differential expression signature, and thus of the associated cell state transition, by assessing the enrichment of their direct targets (regulons) in differentially expressed signature genes (i.e., in genes that are over- or underexpressed during metastatic progression, in this case). Thus, the ARACNe-inferred regulon of each regulatory protein (Table S2) is used as a highly multi- plexed, endogenous reporter for its role in physically controlling metastatic progression.

We used the statistical significance, estimated by sample permutation analysis, to rank-sort the list of putative MRs of the metastatic phenotype in ER+ breast cancer and TNBC (Table S3). Given the strong concordance between ER+ and TNBC MET-GES and between candidate MRs of the two subtypes (Figures S1B–S1D), we inferred the MRs of metastatic progression of breast carcinoma regardless of subtype by VIPER analysis of a subtype-agnostic MET-GES. This generated a single ranked list of breast cancer metastatic progression candidate MRs (breast cancer [BRCA]), indepen- dent of hormonal status

A breast carcinoma context-specific network model of transcriptional regula- tion was assembled with the ARACNe, based on 851 RNA-seq expression pro- files obtained from TCGA. ARACNe was run with 100 bootstraps, a p value threshold of 10?8, and 0 data processing inequality (DPI) tolerance, generating a network of 1,748 TFs associated with 18,783 target genes by 459,569 inter- actions. The regulatory models were generated from the ARACNe results using the VIPER package from Bioconductor (http://bioconductor.org/packages/ release/bioc/html/viper.html)

The gene expression signatures for 20 ER+ and 11 TNBC metastases (MET-
GES) were computed with paired Student’s t test by comparing their profiles against the matching primary tumor ones. Then, the enrichment of each regula- tory protein regulon on the MET-GESs was inferred by the VIPER algorithm (Al- varez et al., 2016; Aytes et al., 2014), as implemented in the VIPERpackage for R available fromBioconductor (https://www.bioconductor.org/packages/release/ bioc/html/viper.html). Statistical

For single patient-based analysis, gene expression signatures were computed by comparing each MET expression profile with the matching pri- mary tumor expression profile. A null model for statistical testing was gener- ated by permuting the samples uniformly at random 1,000 times



### Guo et al. 2018

identify the personalized-sample driver genes from the cancer omics data due to the lack of samples for each individual. To circumvent this problem, here we present a novel **single-sample controller strategy (SCS)** to identify personalized driver mutation profiles from net- work controllability perspective

SCS integrates mutation data and expression data into a reference molecular network for each patient to obtain the driver mutation profiles in a personalized-sample manner

The key idea of SCS is to detect those mutated genes which can achieve the transition from the normal state to the disease state based on each individual omics data from network controllability perspective

The MATLAB-package for our SCS is freely available from http:// sysbio.sibcb.ac.cn/cb/chenlab/software.htm

d Single-sample Controller Strategy (SCS) to assess the im- pact potential of gene mutations on the changes in gene expression patterns. Intuitively, we consider mutations as controllers and gene expression profiles as states in a network, and thus SCS aims to de- tect a small number of mutation genes (i.e. driver genes) which can achieve the transition from the normal state to the disease state from the network controllability viewpoint, based on each individual gene expression data. Our

integrates the mutation data and expression data to a gene-gene regulation network for each pa- tient. I

we aim to identify the minimal number of Individual mutations to control the Individual differentially ex- pressed genes (DEGs) in the Individual gene network

The main steps of SCS method include: 

(1) to obtain personalized DEGs for each patient by comparing the expression profile of the tumor sample with that of the corresponding normal sample, and then extract the individual gene regulations from the expert-curated databases; 

(2) to identify the minimal number of in- dividual mutations with network controllability on the maximal coverage of individual DEGs in the individual gene network and 

(3) based on the dynamic network control theory, to rank and select the driver genes from individual mutations according to the uncovered consensus modules consisting of confidence-weighted paths from the driver genes to the target genes on the gene network. 

SCS considers the cancer- specific or patient-specific mutated network (topology) for predict- ing individual driver genes.

SCS ranks poten- tial driver genes based on their influence on the overall differential expressions of its downstream genes in the individual molecular net- work instead of the collective molecular network.

apply the Condorcet method (Pihur et al., 2008) to determine the sum- mary ranking of genes in a patient population

SCS classifies driver genes regardless of mutation fre- quency, it allows us possibly to discover rare (infrequent) driver genes. Finally,

to bridge the personalized driver mutation discovery problem and the structural network controllability prob- lem. 



### Saelens et al. 2018

numerous alternative module detection methods have been proposed, which improve upon clustering by handling co-expression in only a subset of samples, modelling the regulatory network, and/or allowing overlap between modules. In this study we use known regulatory networks to do a comprehensive and robust evaluation of these different methods

decomposition methods outperform all other strategies, while we do not find a clear advantage of biclustering and network inference-based approaches on large gene expression datasets. Using

Modules in this context are defined as groups of genes with similar expression profiles, which also tend to be functionally related and co-regulated. 

The most popular approach, clustering, has been used since the first gene expression datasets became available and is still the most widely used to this day6–8,10. However, in the context of gene expression, clustering methods suffer from three main drawbacks.



First, clustering methods only look at co-expression among all samples. As transcriptional regulation is highly context specific12, clustering potentially misses local co-expression effects which are present in only a subset of all biological samples. 

Second, most clustering methods are unable to assign genes to multiple modules. The issue of overlap between modules is especially problematic given the increasing evidence that gene regulation is highly combina- torial and that gene products can participate in multiple path- ways13,14. A

A third limitation of clustering methods is that they ignore the regulatory relationships between genes. As the	variation in target gene expression can at least be partly explained by variation in transcription factor expression15, including this information could therefore boost module detection.

Decom- position methods16 and biclustering17
try to handle local co-
expression and overlap. These methods differ from clustering because they allow that genes within a module do not need to be co-expressed in all biological samples, but that a sample can influence the expression of a module to a certain degree (decomposition methods) or not at all (biclustering methods).

Two other alternative methods, direct network inference15 (direct NI) and iterative NI18, use the expression data to additionally model the regulatory relationships between the genes.

However, because clustering methods do not detect local co-expression effects, they could potentially miss relevant modules or exclude important genes from a module. In use cases where it can be desirable that all modules are discovered in a dataset, e.g., to generate signatures for disease, therapy and prevention4,11,32,orto find a set of genes responsible for a biological function, methods that detect such local co-expression and/or overlapping modules could therefore provide a substantial advantage. Consistent with this, we found that decomposition methods based on ICA were better at reco- vering known modules consistently across datasets.

a third major application of co- expression modules is in the inference of gene regulatory net- works, where modules can be used to improve the network by combining information from several genes33 but can also improve the ease of interpretation. When we assessed the accuracy of the inferred network by combining a state-of-the-art network infer- ence methods with different module detection methods, we found that ICA-based decomposition methods lead to the highest improvement in accuracy (primarily



### Li et al. 2018

We propose a network-based method to detect cancer specific driver modules (CSDM) in a certain cancer type to other cancer types.

We construct **the specific network of a cancer by combining specific coverage and mutual exclusivity** in all cancer types, to catch the specificity of the cancer at the pathway level.

In this work, we propose a network-based method to detect cancer specific driver modules
(CSDM), which can catch the specificity of a certain cancer type to other cancer types at the pathway level. A cancer specific driver module must have high coverage and high exclusivity in a certain cancer, and a higher percentage of samples in this cancer than other cancer types. We first construct the specific network for a certain cancer type by integrating specific coverage and mutual exclusivity in all cancer types. Then, we use a greedy algorithm to detect all of the specific driver modules in the specific network. 

We first construct the specific network for a certain cancer by integrating specific coverage
and mutual exclusivity. Then, we use a greedy search to detect the cancer specific driver modules. The overview of our method is shown in Figure 6. Finally, we utilize the specific coverage and the significance of mutual exclusivity to evaluate the cancer specific driver modules

The specific coverage for each gene pair in a certain cancer type is proposed to catch the specificity
of the cancer specific driver module. Then, the mutual exclusivity is used to quantify the exclusivity of each gene pair in the same driver module. Finally, a cancer specific network is constructed by combining the specific coverage and mutual exclusivity

Cancer Specific Driver Module Detection We first define the internal coverage, external coverage, and specific coverage of the specific
driver module in cancer k.

We use the greedy algorithm to detect the specific driver modules from the specific network N(k)
of cancer k. 

We apply this greedy algorithm on each gene in specific network N(k) and consider the modules
that have at least three genes as the specific driver modules. This greedy algorithm protects the gene pairs with high external coverage in cancer k, which guarantees high specific coverage of driver modules in a certain cancer to other cancer types.





### Boltz et al. 2019

determined **collective influencers (CI)**, defined as network nodes that damage the integrity of the underlying networks to the utmost degree

Morone and Makse22 introduced **an optimization method to determine an optimized set of nodes termed collective influencers (CI). Such nodes were obtained via optimal percolation theory through the investigation of their propensity to damage the underlying network, strongly emphasizing the role of weakly con- nected nodes**. Wondering if collective influencers in protein-protein interaction network carry biological signif- icance, we expected that collective protein influencers were enriched with e.g. disease or essential genes. Within the human interactome, CI proteins were indeed enriched with essential, regulatory, signaling and disease genes as well as drug targets, strongly suggesting that such well-defined protein groups have significance. Furthermore, we found that CI proteins were evolutionarily conserved as CI proteins in evolutionarily conserved networks in different organisms.
Results

, a set of collective influencers (CI) is defined as the minimum set of nodes that, upon deletion, destroy the largest connected component of the underlying network22

. Calculating a score for each protein that reflects its propensity to damage the underlying largest connected component, proteins with the largest score were removed in each step. The procedure stopped when the largest connected component disappeared, providing a list of removed nodes as collective influencers (CI). W



The determination of collective influencers is based on optimal percolation, aiming at the determination of a minimum set of nodes that fragments the underlying network.



The collective influence theory for optimal percolation is based on the message passing equations of the per-
colation process. For

ref22 = Morone, F. & Makse, H. A. Influence maximization in complex networks through optimal percolation. 



### Yuan et al. 2019

we propose an approach, which use **biweight midcorrelation to measure the correlation between factors and make use of nonconvex penalty based sparse regression for gene regulatory network inference (BMNPGRN)**. BMNCGRN incorporates multi-omics data (including DNA methylation and copy number variation) and their interactions in gene regulatory network model. The experimental results on synthetic datasets show that BMNPGRN outperforms popular and state-of-the-art methods (including DCGRN, ARACNE and CLR) under false positive control. Furthermore, we applied BMNPGRN on breast cancer (BRCA) data from The Cancer Genome Atlas database and provided gene regulatory network.



we propose an approach, which use
Biweight Midcorrelation to measure the correlation between factors and make use of Nonconvex Penalty based sparse regression for Gene Regulatory Network inference (BMNPGRN). BMNPGRN integrates heterogeneous multi-omics data to infer gene regulatory network with gene expression data, DNA methylation data and copy number variation data. In order to infer gene regulatory network. Firstly, we combine biweight midcorrelation coefficient algorithm, which is an efficient algorithm for computing correlation coefficient, with ‘differential correlation strategy’ to learn associations among DNA methylation sites. Then, nonconvex penalty based sparse regression is used to find gene-related biological factors, and the parameter of method is determined by cross-validation. Meanwhile, nonconvex penalty based sparse regression is used under stability selection which can control false positives effectively [29]. Finally, BMNPGRN identifies gene regulatory network based on the probabilities of biological factors (i.e. gene, DNA methylation sites, and CNV)



Our proposed approach BMNPGRN has advantages
over existing gene regulatory network inference methods. Firstly, BMNPGRN can find more DNA methylation sites which are associated with gene regulatory network. Such method provide deeper insight into gene regulation mechanism. Secondly, BMNPGRN can effectively control false positives using stability selection strategy. Furthermore, BMNPGRN can more accurately find biological factors using nonconvex penalty based sparse regression. Finally, to the best of our knowledge, BMNPGRN is the first method which is applied to breast cancer data obtained by high-throughput sequencing technology







### Matsubara et al. 2019

merges proteome and transcriptome data 

**spectral clustering** of the network and map resulting eigenvectors from the 2nd and 3d smallest eigenvalues into 2D space which preserves topological distance and cluster structure

map gene expression profiles onto this mapping->generate a set of image-like representation of protein networks processed by **deep deconvolutional layers**

PPIs from HINT database



### Wani et al. 2019 - review

Network based integration and inference Network based integration approaches provide the easiest and
simplest way to integrate different types of network data into a single unified representation. It

**Similarity Network Fusion (SNF) (**Wang et al., 2014) is a typical
example of constructing merged networks by employing network based data integration principal. SNF integrates expression data (mRNA, miRNA) and DNA methylation data and creates a network of samples (i.e., patients) for each data type and then builds an integrated network from the fusion of individual networks.

Integrative gene regulatory network (iGRN) (Zarayeneh et al. 2017 , is another similar **network integration technique** that combines gene expression data with CNV and DNA methylation data.



### Grimes et al. 2019

A general framework for integrating known gene regulatory pathways into a differential network analysis between two populations is proposed. The framework importantly allows for any gene-gene association measure to be used, and inference is carried out through permutation testing. A simulation study investigates the performance in identifying differentially connected genes when incorporating known pathways, even if the pathway knowledge is partially inaccurate



network with PCC The resulting networks will often be very dense, and the natural interpretation of the edges - that of a direct connection - is not valid; the edges represent marginal, not direct or causal connections

the partial correlation does reveal direct connections

If the gene expression profiles follow a multivariate Gaussian distribution, then two genes have non-zero partial correlation if and only if they are conditionally dependent given the other genes16. The resulting network of conditionally dependent nodes is often called a Gaussian graphical model (GGM). 

The incorporation of pathway information has been shown to improve performance in differential expres-
sion (DE) analyses. Researchers have used KEGG pathways29 to construct a Markov random field (MRF) that improves performance in finding DE genes30. Others used KEGG pathways to inform a spatially correlated mix- ture model for DE analysis31. Reactome pathways32 have also been utilized in DE analysis; in one example altered pathways in lung adenocarcinoma and colon cancer were identified33. The integration of pathway information into DE analysis is an on-going area of research3

This report proposes a framework for integrating known genetic pathways into a differential network analysis
of two populations. The framework allows any association measure to be used, and a general measure for differen- tial connectivity is considered. Statistical significance is evaluated through a permutation testing procedure. The methodology is implemented in R and is available on GitHub at https://github.com/tgrimes/dnapath

1. Given a pathway G, compute the estimated association networks, ˆS and ˆS , for the two groups. 
2. Evaluate the differential connectivity scores, δE, for eachEG∈ ?() of interest.

3. Assess the statistical significance of these scores using a permutation testing procedure.
4. Repeat steps 1–3 for each pathway
   G ∈ ?

Estimation of association networks. The first step is to estimate the gene-gene association network within each group.

Differential connectivity score. The differential network analysis measures the change in a set of connections
EG |∈ <
⊂= ?
() {( ,) ,, } for a given pathwayG ∈ ?. We generalize the differential connectivity score
ij ij Gi j
proposed in earlier work44 by the p-norm of the difference in connectivity scores in E

tests for significance: This sets up a permutation testing procedure that can be used to estimate a p-value for d under the null, whereby permutations of the group labels are used, i.e. the observations are shuffled between groups. The total number of distinct permutations will often be quite large even for moderate sample sizes.



### Gysi et al. 2020

a method for the systematic comparison of an unlimited number of networks, with unlimited number of transcripts: **Co-expression Differential Net- work Analysis (CoDiNA).**

CoDiNA detects links and nodes that are common, specific or different among the networks

. Several meth- ods for pairwise network comparisons exist, for example: CoXpress [21], CSD [22], DCGL [20, 23], DICER [14], DiffCoEx [24], DiffCorr [15, 25], Gain [26], MIMO [27], ModMap [28], NetAlign [29], SAGA [30], the discordant method [31] and QNet [32]. I

It is often aimed for a comparison ofmore than two networks simultaneously, such as gene co-expression networks arising from different species, tissues or diseases, or co-existence net- works from different environments. Existing methods for contrasting multiple networks focus on identifying modules ofdifferentially co-expressed genes [1, 20, 23, 34–36], thereby allowing the identification ofgene groups that are functionally related. 



CoDiNA requires as input a set ofindependently constructed undirected networks to be com-
pared (Fig 1a, 1b and 1c). These could be any kind ofundirected—weighted or unweighted— networks (e.g., protein-protein interaction networks, metagenomics networks, co-occurrence networks etc.)

The CoDiNA method can be divided into five steps

1. Remove nodes that were not measured in all networks to be compared; 
2. Define a minimum cutoff for the weight a link needs to have to be considered present
3. Remove links that are absent in all network; 
4. Classify (and sub-classify) and score links as common, specific or different between net- works; denoted as F and ~F, respectively;
5. Classify the nodes.

input: correlation network or or can also be Consensus Networks (CN) derived from similar networks to achieve higher confidence in the network inference



### Trilla-Fuertes et al. 2020

n previous works, our group defined a new hormonal receptor positive subgroup, the TN-like subtype, which had a prognosis and a molecular profile more similar to triple negative tumors. In this study, prote- omics and Bayesian networks were used to characterize protein relationships in 96 breast tumor samples. Components obtained by these methods had a clear functional structure.

**Undirected probabilistic graphical models (PGM),** based on a Bayesian approach, allow characterizing differences between tumor samples at functional level [2, 3, 6, 7]. In this study we explored the utility ofBayesian networks in the molecular characterization ofbreast cancer.

The main feature oftargeted Bayesian networks is that they provide a hierar- chical structure and targeted relationships between proteins



PGM are graph-based representations ofjoint probability distributions where nodes represent random variables and edges (directed or undirected) represent stochastic dependencies among the variables. In particular, we have used a type ofPGM called Bayesian networks (BN)

With these models, the dependences between the variables in our data are specified by a directed acyclic graph (**DAG**). The obtained networks will indicate causality i.e. if protein A and B are connected and protein A changes its expression value, protein B changes its expres- sion value as well

we find the BN that best explains our data [23]. There are different algorithms to
learn a DAG from data but we have selected the well-known PC algorithm (named as its inven- tors Peter Spirtes and Clark Glymour), a constraint-based structure learning algorithm [24] based on conditional independence tests. The PC algorithm was shown to be consistent in high-dimensional settings [25]. Moreover, an order-independent version ofthe PC algorithm, called PC-stable, was proposed in [26]. All these procedures are implemented in R within packages pcalg [25] and graph [27]. We used protein expression data without other a priori information.

our data are represented by a large graph that can be partitioned into several connected components. Then, we focused on finding suitable subgraphs that give us a much clearer understanding ofthe interrelations therein

Using proteomics data, directed acyclic graphs (DAG) were performed. Altogether, it was pos- sible to establish 789 edges ofwhich 662 were guided and 127 are undetermined

We characterized components from DAG analysis. Components including less than 9
nodes were dismissed because they were little informative. All components were named with the number ofnodes included by the DAG analysis. Afterwards, components were interrogated for biological function

Component activities were calculated for each node. There were significant differences between ER-true, TN-like and TNBC tumors in the component activity for component 23: mitochondria, component 17: RNA binding, component 13: extracellular matrix, and compo- nent 10: extracellular exosome (Fig

we used proteomics and DAG to characterize relationships between proteins in breast cancer tumor samples

our DAG method supplies directed relationships between proteins and a hierarchical structure. Tradi- tionally, protein-protein interaction (PPI) networks, such as STRING, are based in relation- ships described in the literature. However, we built a directed network, i.e. a graph formed by edges with a direction, using protein expression data without other a priori information, so it was possible to propose new hypotheses about protein interactions. We used probabilistic graphical models (PGM) because they offer a way to relate many random variables with a com- plex dependency structure



arrows in directed networks indicate causality
between two proteins. This approach allows making hypotheses about causal relationships between proteins and proposes a hierarchical structure.



### Li et al. 2020

To explore the dynamics of the mammalian cellular aging
network, we employ non-linear differential equations (Tyson and Novák, 2010) to describe the dynamics of each genes expression in the network. A

A biological system is naturally subject to intrinsic and extrinsic fluctuations. Therefore, we added an additional fluctuation term in the ODEs to characterize the stochastic behaviors of the mammalian cellular aging process. We use the Langevin dynamic approach to simulate the gene circuit dynamics. From the resulting dynamic trajectories of the gene expressions, we collected the statistics and quantified the underlying potential landscape

, we presented a mathematical model to describe the dynamic features of the mammalian cellular aging process. We built the underlying gene regulatory network by integrating the information from previous experimental studies. The genes and wirings in the gene regulatory network were formed, and the dynamics of gene expression was described by nine non- linear ordinary differential equations. 

we have provided a framework to reveal the
underlying mechanism of fast-aging and slow-aging in mammals based on **landscape and flux theory**. We predict the key genes and interactions in the fast-aging and slow-aging processes.



### Lee et al. 2020

**heterogeneous multi-layered network (HMLN) has proven successful in integrating diverse biological data for the representation of the hierarchy of biological system**

relation inference
(or link prediction) can be formulated as a problem in a homogenous network, a multiplex network, or a heterogeneous multi-layered network (HMLN), as shown in Figure 1.In a homogenous network (Figure 1A), all nodes from different domains, as well as intra-domain and cross-domain relations, are treated equally. In contrast, multiplex and multi-layered networks separate different types of nodes and relations. A multiplex network is often used to represent homogeneous nodes that have different types of characterizations (a.k.a. views). For example, a gene can be characterized by multiple measurements ofgene expression, essentiality, literature citation, phylogenetic profile, neighborhood in the interaction network, biological pathway involved, Gene Ontology annotation, protein domain profile etc. (Hwang et al., 2019). Each type of measurement can form a unique type of link between genes (Figure 1B). In a HMLN (Figure 1C), multiple types of heterogeneous nodes are involved. The nodes from each type are grouped into a single layer and treated separately. In the same vein, different types of intra-domain and cross-domain relations are marked differently in a multi-layered network. We note that more complex network representations, such as multiplex multi- layered network, may be needed in real applications. In this review, we focus on the cross-domain relation inference (or link prediction) problem for the HMLN. Readers can refer other excellent reviews of the multiplex networks (Chauvel et al., 2019)

Compared to a homogeneous single-layered network, a unique topological characteristic of a multi-layered network lies in its cross-layer relation or dependency structure in addition to intra- layer connectivity.

The cross-layer relation inference problem is conceptually related to collaborative filtering (Goldberg et al., 1992). Commonly used collaborative filtering methods can be classified into two groups: neighborhood methods (Breese et al., 1998) and latent factor methods (Koren et al., 2009). As the latent factor approach is generally more effective in capturing the implicit cross-layer relations, many variants ofthis methodology, such as recommended systems (Portugal et al., 2018), have been proposed to address relation inference problems in a two-layered network (Gao et al., 2019; Xuan et al., 2019). However, few methods have been developed for the multi-layered network.



Besides the traditional algorithms, like matrix factorization, random walk, and meta-path, introduced in previous sections, the embedding of HMLN can also benefit from Deep Learning techniques, especially the Neural Networks (NNs). Though NNs are initially proposed to learn the embedding of data, such as texts, images, and videos, they have shown powerful performance when dealing with graph structured data, which exist in non-Euclidean domain. Due to the growing interests and demands in recent years, Graph Neural Networks (GNNs) have been proposed to learn the embedding of graphs



A GNN consists of a number of hidden layers that employ iterative, propagation procedures in order to transform different node and edge features. Each layer takes the output of the previous layer as the input. With graph structured data, GNNs adopt element (node or edge) features X and the graph structure A as input to learn the representation of each element hi,or graph hG, for different tasks. Each hidden layer employs the
“aggregation” functions and the “update” functions (Battaglia et al., 2018). Each aggregation function r takes a set of node or edge features as input and reduces it to a single element which represents the aggregated information. The aggregations usually operate on the nearest neighbors or the local subgraphs of each element to capture local information gradually. Since the permutation invariance of the input holds in graph data, the r functions must also have the same property. These functions can take variable numbers of arguments. Commonly used r functions include sum (Xu et al., 2019), mean (Kipf and Welling, 2017), max-pooling (Hamilton et al., 2017)and attention mechanism (Velickovic et al., 2018; Wang et al., 2019; Fan et al., 2019). Update functions f are applied across all elements to compute per-element updates after the aggregations. In the final layer, the generated embedding can be fed into the classification/prediction layer, and the whole model is trained for different (e.g. node classification, link prediction) tasks



the meta-path based GNN shares the same limitations of other meta-path based algorithms. New types of GNNs, those that explicitly take different types of relations into consideration, are needed for the link prediction problem in HMLN



Yao et al. integrated multi-omics data to construct a three-
layered network model MetPriCNet, which consists of metabolite network, gene network, phenotype network, metabolite-phenotype network, metabolite-gene network, and gene-phenotype network (Yao et al., 2015). Afterwards, an RWR algorithm is applied to prioritize metabolites associated with diseases



the RWRmethod
has been used to identify other molecular dysregulations that are associated with diseases based on the multi-layered network model. To

### Baur et al. 2020 - review

there are two major approaches to model gene regulatory networks [15]. 

1- The first approach predicts regulatory proteins, for example, TFs, chromatin remodelers, and signaling proteins that determine a gene’s expression level (Figure 1b). 

2- The second approach identifies CREs affecting gene expression levels (Figure 1c), that is, TF binding sites in or near the promoter, or distal enhancers interacting with the promoter through chromatin loop- ing [2].

CRE=cis-regulatory element

A context-specific regulatory network is defined as a network that captures interactions between regulatory proteins and target genes in a particular context, such as a developmental stage, tissue, disease, or environment

While there are a large number of mathematical formalisms to represent regulatory net- works [16], we focus on statistical models, which include information-theoretic methods [17], dependency net- works [18,19], Bayesian networks [20], and Gaussian graphical models [21]. Given

*<u>Incorporating prior to constrain network structure and estimate regulator activity</u>* 

Prior-based regulatory network inference can use either 

* a structure prior approach
  * prior probability distribution on the graph structure, such that an edge with motif support is more likely to be included in the regulatory network model
  * examples: For example, Chasman et al. [23] used DNase-seq filtered motifs as a structure prior for a network inference algorithm [18] and context-specific gene expression data of hindbrain and spinal cord development to infer a regulatory network in neuroepithelial stem cells. 
* a parameter prior approach 
  * regularized regression setting [24], which learns a sparse network by penalizing edge addition while optimizing prediction of a target’s expression level
  * the prior can be used to reduce the penalty for addition of an edge in the network [24]



Newer methods have used the prior network to addi- tionally estimate hidden TFactivity (TFA) to overcome the assumption that mRNA levels must correlate with the regulator’s protein activity on a gene’s promoter [24,27]. **TFA estimation takes in an initial, noisy regu- latory network and expression and models the expres- sion as a product of the hidden TFAs and a refined network structure** [27]. The TFA can be estimated independently followed by network inference, for example, modified least absolute shrinkage and selec- tion operator with stability approach to regularization selection (mLASSO-StARS) [24], or in a single iterative algorithm, for example, Network Reprogramming using EXpression (NetRex) [27]. mLASSO-StARS [24] learns a dependency network by solving a set of independent Least Absolute Shrinkage and Selection Operator (LASSO) regression problems, one per gene. Both mRNA and TFA are used as potential predictors of a target gene’s expression. NetRex [27] jointly estimates both the TFA and the network to explain the expression data by using regularization to penalize the mismatch in the number of edges of the inferred and prior networks. NetRex updates the prior network until the rewired network and the estimated TFA optimally explains the expression data. Both NetRex and mLASSO-StARS need to specify the hyperparameters to obtain a final network structure. mLASSO-StARS uses StARS to select the hyperparameters by controlling the average instability of edges. NetRex uses a grid search over hyperparameters and obtains the final result from the consensus of multiple nearly optimal settings. NetRex has additional constraints that can better capture core- gulatory relationships in gene modules but can make network inference computationally expensive. These approaches were shown to improve the quality of inferred networks and were useful for identifying key regulators and targets in mammalian [23,24] and insect [27] systems



<u>*Network inference across multiple contexts*</u> 

Often data from multiple contexts are available and a key question is to define regulatory networks for each context and identify similarities and differences. Although one approach would be to infer networks for each context independently [28], examining multiple contexts simultaneously could improve the quality of the inferred networks [21,29], make comparisons easier [30,31], and help in scenarios with low sample sizes [32]. To enable simultaneous network inference across multiple conditions, several approaches have used multi-task learning (MTL) (Figure 2a) [33]. In MTL, related tasks are jointly solved while sharing the infor- mation across tasks. In one class of methods, a regulatory network is modeled as a Gaussian Graphical Model (GGM), which represents the genome-wide expression levels of genes as a multivariate Gaussian and the network structure by a precision matrix. One GGM approach, JRmGRN [21], models each condition’s network as sum of two sparse precision matrices, one for condition-specific components and the second for shared component across all conditions. JRmGRN additionally regularized the shared network to favor more hub genes and was more effective at identifying shared hubs and context-specific network components than existing approaches. While JRmGRN models each condition with one precision matrix, another GGM- based approach, NETI2 [31], additionally models the intrasample heterogeneity, which is common in tumor samples because they represent a mix of cancerous and non-cancerous cells. NETI2 infers cancer subtypee specific networks and a network for non-cancerous cells shared across the subtypes. The proportion of cancerous and non-cancerous cells for each sample is known; however, the expression of each cell type is not known and is estimated using an expectationemaximization algorithm. On The Cancer Genome Atlas (TCGA) data [34] for different breast cancer subtypes, NETI2 captured subtype-specific gene hubs that exhibited reduced connections in the non-cancerous network



### Mahajan et al. 2021

We combined these epigenomic datasets and integrated them with the reference human protein interactome using a novel **network propagation** approach.

a network-based approach to integrate multi-omic epigenetic datasets.

At a conceptual level, tumor mutations and features derived from epigenomic datasets drive large-scale cell behavior in similar ways. Individual mutations or chromatin features may exert a weaker, more localized influence, but groups of mutations or epigenetic marks will exert a stronger, systemic influence by influencing multiple components of the underlying molecular network. Thus, both mutations and epigenetic features represent biological priors that are weak or uninformative when considered in isolation but are strong when analyzed together in the context of the underlying topological network structure. A network-based approach can be used to filter out isolated noisy signals and hone in on biologically significant pathways of genes that mutations or epigenomic changes would target at multiple points.

a novel network-based approach to propagate epigenomic signals
derived from our ATAC-seq and CUT&RUN studies and thus uncover network modules that are likely to drive Tfh cell differentiation. Our approach is based on the premise that biologically relevant epigenetic signatures reflect the need for coordinated regulation of genes that encode interacting proteins. A similar conceptual premise has been used to identify subnetworks enriched for mutations in cancer genomics and network-based GWAS (Cowen et al., 2017; Leiserson et al., 2013, 2015; Reyna et al., 2018; Vandin et al., 2011). The

he underlying technique of network propagation represents a powerful way to combine a wide range of signals taking into account the structure of the underlying network. **We implemented network propagation using random walk with restart, which is equivalent to an insulated heat diffusion process** ((Cowen et al., 2017; Leiserson et al., 2013, 2015; Reyna et al., 2018; Vandin et al., 2011)). **We used this approach to unify features derived from the different epigenomic datasets and discover networks driving Tfh cell differentiation.** 

we incorporated higher-order relationships between proteins encoded by these genes using the high-quality (i.e., each edge is validated experimentally using multiple independent assays) reference human protein interactome network (Cusick et al., 2009; Das and Yu, 2012).



Each of the epigenetic signals analyzed exhibited distinct trends across the four stages of Tfh cell differentiation. These differences likely reflect complementary regulatory processes underlying the differentiation of this polarized cell type (Fig. 1d). We expect that these processes all involve regulating the expression of an underlying set of components that drive the transition from naive to GC Tfh, and multi-omic integration strategies will be required to uncover these core genes. We focused on the discovery of trajectories using the epigenomic data. The transcriptomic data was used as an orthogonal dataset to validate the identified trajectories.

Next, for each differentiation pattern, peaks from each dataset were aggregated into gene-centric scores using a proximity-weighted count heuristic (Fig. 2b) based on first principles underlying the organization of promoters and enhancers. Gene-centric scores from each dataset were then combined by weighting the datasets equally, as we did not have any a priori information regarding relative information context across the datasets. We refer to each gene’s score as its “peak proximity score” (PPS), which we can calculate for each possible differentiation pattern and use as a surrogate for gene regulation

**We ran network propagation through HotNet2 using the combined PPS’s for each protein’s gene as seed values (Figs. 3a, 3b). HotNet2 helps us incorporate the topological structure of the underlying protein interaction network as an informative prior on how proteins are functionally related to one another**

**The use of network propagation is motivated by the hypothesis that proteins encoded by relevant genes relevant to Tfh cell differentiation are more likely to interact with one another and make up cohesive network submodules (Fig. 3b). Further, our approach can help filter out genes that may not play a role in Tfh cell differentiation but nevertheless show similar epigenetic patterns to relevant genes**

**network propagation constitutes the final step in a three-step pipeline: identification of the most likely patterns of Tfh cell differentiation, generation of gene-centric scores using a proximity-weighted count heuristic (i.e. calculation of PPS’s), and network propagation using these scores as seeds.**

**we adapted the HotNet2 algorithm to identify protein network modules enriched for signals derived from ATAC-seq and CUT&RUN datasets obtained from Tfh cells at four stages of differentiation within a healthy human subject’s tonsil. By framing the search for drivers as a search for identifying protein subnetworks with concentrated epigenomic signals, our approach recovered known regulators of Tfh cell differentiation, and we have implicated several novel candidate genes involved in Tfh cell differentiation.** 



network propagation represents a flexible and broadly useful approach for studying biological processes that require the integration of disparate epigenomic datasets. This



Calculation of a gene-centric peak proximity scores
For a biological (EBSeq) pattern of interest, peak proximity scores (PPS’s) for a gene (“gene- centric PPS”) were calculated for a biological pattern of interest by calculating a weighted sum of the number of peak clusters from the pattern in the gene’s vicinity using distance-based weights to reflect the varying contributions of promoter (± 10 kb) [3x weight], proximal enhancer (± 25 kb) [2x weight] and distal enhancers (± 250 kb) [1x weight] to gene regulation. Gene-centric PPS scores were first calculated individually for each epigenomic dataset. These were then combined into a final PPS score using equal weights for the different datasets.

Network propagation using random walk with restart
For each pattern, we used network propagation to integrate the gene-centric scores. We used the random walk with restart algorithm (equivalent to an insulated heat diffusion process) as implemented in HotNet2 (Leiserson et al., 2015) to identify high scoring sub-networks. Relevant HotNet2 code is available at https://github.com/raphael-group/hotnet2. We ran network propagation on the union of the high-quality reference human binary and co-complex protein interactomes (Das and Yu, 2012). 



### Duan et al. 2021

Refuting the widely held intuition that incorporating more types of omics data always produces better results, **our analyses showed that there are situations where integrating more omics data negatively impacts the performance of integration methods.** Our analyses also suggested several effective combi- nations for most cancers under our studies, which may be of particular interest to research- ers in omics data analysis.

we observed that the methods NEMO and SNF perform very well in all three criteria

several commonly-used combinations ofomics data types can indeed improve the accuracy of all the ten methods as measured by both clustering and clinical metrics. On the other hand, our analysis indicates that integrating more types ofomics data may negatively impact the performance on cancer subtyping, refuting the widely held intuition that incorporating more types ofomics data always helps produce better results

using the four omics data types for integration analysis were not always better than the results obtained by using three or two-omics data types. Similarly, the results ofthe three-omics integration were not always bet- ter than the results oftwo-omics integration. While this observation is counter-intuitive and deserves further investigation, we believe that it is the consequence ofthree intertwined fac- tors: (1) the negatively-correlated noises in the omics data which may cancel out useful infor- mation; (2) the redundancy in different types ofomics data; and (3) the computational/ statistical challenges that integrating excessive datasets impose on the integration methods pre- venting them from making the best use ofthe information, if any, in the multi-omics datasets, to calculate optimal solutions.

another surprising observation on the significance ofDNA methylation data in the effectiveness ofintegration methods

**Our analysis on the effectiveness ofthe 11 possible combinations, however, shows that only one combination in 2 and 3 omics datasets can be considered to be effective.** This is unexpected as methylation had been proven to play an important role in cancer and has been the most common data type used in previous research on integrated multi-omics data for cancer subtyping. We leave it to future work to understand this unexpected phenomenon, but speculate that it is the result of the following three factors ofour experiments: (1) certain characteristics ofthe methylation data that do not fit the model assumption ofsome ofthe ten integration methods, resulting in a much lower overall effectiveness score for those ofthe methylation-participating combina- tions in our evaluation; (2) the data processing step in which we mapped all the features to genes in methylation and CNV datasets as previous studies did, resulting in a loss ofsome sig- nificant information; and (3) performance criterion that uses clinical-based metrics as well as clustering-based metrics, in comparison to previous studies where only one metric may be used

### Paull et al. 2021 - MOMA

, using a network-based approach, identified 407 master regulator (MR) proteins responsible for canalizing the genetics of individual samples from 20 cohorts in The Cancer Genome Atlas (TCGA) into 112 transcriptionally distinct tumor subtypes

. MR proteins could be further organized into 24 pan-cancer, master regulator block modules (MRBs), each regulating key cancer hallmarks and predictive of patient outcome in multiple cohorts. Of

Of all somatic alterations detected in each individual sample, >50%were predicted to induce aberrant MR activity, yielding insight into mechanisms linking tumor genetics and transcriptional identity and establishing non-oncogene dependencies.

**step1: gene expression profiles from 20 TCGA cohorts (Table S1) were first transformed to protein activity profiles by using the Virtual Proteomics by Enriched Regulon Analysis (VIPER) al- gorithm**

Candidate MR proteins were then identified by Fisher’s integration of p values from (1) their VIPER-measured activity, (2) functional genetic al- terations in their upstream pathways by Driver-Gene Inference by Genetical-Genomic Information Theory (DIGGIT) analysis (Chen et al., 2014), and (3) additional structure and literature- based evidence supporting direct protein-protein interactions between MRs and proteins harboring genetic alterations, via the Predicting Protein-Protein Interactions (PrePPI) algorithm (Zhang et al., 2012) (steps 2 and 3) (Figure



vector of integrated (Log_10 p)2 (MOMA scores) to weigh each MR’s contribution in a tumor subtype clustering step (step 4) (Figure S1D). Finally, genomic saturation analysis upstream of top candidate MRs identified those most likely to control the subtype transcriptional identity (step 5) (Figure S1D). Finally, this was followed by identification and functional characteriza- tion of MR block sub-modules, termed MR-Blocks (MRBs), recurring across multiple subtypes (step 6) (Figure S1E). See STAR methods for a detailed description of each step

**VIPER has been extensively validated as an accurate method-**
**ology to measure a protein’s activity, on the basis of the enrichment of its tissue-specific activated and repressed tran- scriptional targets (regulon) in over and under-expressed genes (Alvarez et al., 2016)—i.e., akin to a highly multiplexed gene- reporter assay. To generate accurate regulons for 2,506 regulatory proteins annotated as transcription factors (TFs) and co-transcription factors (co-TFs) in Gene Ontology (GO) (Ash- burner et al., 2000; The Gene Ontology Consortium, 2019), we used Algorithm for the Reconstruction of Accurate Cellular Net- works (ARACNe) (Basso et al., 2005); see STAR methods for ARACNe and VIPER accuracy metrics**

**For each candidate MR, we first identified candidate upstream modulator proteins by using the Conditional Inference of Network Dynamics (CINDy) algorithm (Giorgi et al., 2014) and then assessed whether the presence of genomic alterations in their encoding genes was associated with differential MR activity (activity Quantitative Trait Locus analysis [aQTL]). These two steps comprise the DIGGIT algorithm, which was highly effective in elucidating key driver mutations missed by prior analyses in GBM**

. Minimum cohort size reflected the need to generate accurate regulatory network models by using the ARACNe algorithm



### Cao et al. 2021

we 5 propose a **computational framework called GLUE (graph-linked unified embedding), which utilizes 6 accessible prior knowledge about regulatory interactions to bridge the gaps between feature spaces**

**we introduce GLUE (graph-linked unified embedding), a modular framework for integrating 28 unpaired single-cell multi-omics data and inferring regulatory interactions simultaneously. By 29 modeling the regulatory interactions across omics layers explicitly, GLUE bridges the gaps between 30 various omics-specific feature spaces in a biologically intuitive manner.**

**Integrating unpaired single-cell multi-omics data via graph-guided embeddings**

we model cell states as low-dimensional cell embeddings learned 7 through variational autoencoders27, 28. Given their intrinsic differences in biological nature and assay 8 technology, each omics layer is equipped with a separate autoencoder that uses a probabilistic 9 generative model tailored to the layer-specific feature space

**Taking advantage of prior biological knowledge, we propose the use of a knowledge-based graph 2 (“guidance graph”) that explicitly models cross-layer regulatory interactions for linking layer- 3 specific feature spaces**; the vertices in the graph correspond to the features of different omics layers, 4 and edges represent signed regulatory interactions.

, adversarial multimodal alignment is performed as an iterative optimization procedure, guided 8 by feature embeddings encoded from the graph2

**GLUE employs omics-specific variational autoencoders to learn low-dimensional cell embeddings from each omics 14 layer.** The data dimensionality and generative distribution can differ across omics layers, but the cell embedding 15 dimensions are shared. A graph variational autoencoder is used to learn feature embeddings from the prior 16 knowledge-based guidance graph; these embeddings are then used as data decoder parameters. The feature 17 embeddings effectively link the omics-specific autoencoders to ensure a consistent embedding orientation. Last, an 18 omics discriminator is employed to align the cell embeddings of different omics layers via adversarial learning.

Combining omics-specific autoencoders with graph-based coupling and adversarial alignment, we 24 designed and implemented the GLUE framework for unpaired single-cell multi-omics data 25 integration with superior accuracy and robustness. By modeling regulatory interactions across omics 26 layers explicitly, GLUE uniquely supports model-based regulatory inference for unpaired multi- 27 omics datasets, exhibiting even higher reliability than regular correlation analysis on paired datasets

The whole package of GLUE, along with tutorials and demo 30 cases, is available online at https://github.com/gao-lab/GLUE for the community.

-------





##### Margolin et al. 2006  - ARACne

This method uses an information theoretic approach to eliminate the majority of indirect interactions inferred by co- expression methods.

ARACNE reconstructs the network exactly (asymptotically) if the effect of loops in the network topology is negligible, and we show that the algorithm works well in practice, even in the presence of numerous loops and complex topologies

ARACNE (Algorithm for the Recon- struction of Accurate Cellular Networks), a novel informa- tion-theoretic algorithm for the reverse engineering of transcriptional networks from microarray data

ARACNE defines an
edge as an irreducible statistical dependency between gene expression profiles that cannot be explained as an artifact of other statistical dependencies in the network



Within the assumption of a two-way network, all statisti- cal dependencies can be inferred from pairwise marginals, and no higher order analysis is needed. While not imply- ing that this is always the case for biological networks, it is important to understand whether this assumption may allow the inference of a subset of the true interactions with fewer false positives. Thus we identify candidate interac- tions by estimating pairwise gene expression profile mutual information



We then filter MIs using an appropriate threshold, I0, computed for a specific p-value, p0, in the null-hypothesis of two independent genes. This step is basically equivalent to the Relevance Networks method [6] and suffers from the same significant limitations; namely, genes separated by one or more intermediaries (indirect relationships) may be highly co-regulated without implying an irreduci- ble interaction, resulting in numerous false positives

Thus in its second step, ARACNE removes the vast major-
ity of indirect candidate interactions (φij = 0) using a well- known information theoretic property, the data process- ing inequality (DPI, discussed in detail later), that has not been previously applied to the reverse engineering of genetic networks



##### Hudson et al. 2009

a new algorithm that correctly identifies the gene containing the causal mutation from microarray data alone. 

a coordinated, simultaneous, weighted integration of three sources of microarray information: transcript abundance, differential expression, and differential wiring. 

Our new approach identifies causal regulatory changes by globally contrasting co-expression network dynamics. The entirely data-driven ‘weighting’ procedure emphasises regulatory movement relative to the phenotypically relevant part of the network.

we examined the difference in the specific behaviour or co-expression of targeted pairs of genes between the two crosses, by subtracting the correlation coefficient in Wagyu from that in Piedmontese. This

n an attempt to assess the importance of each DE gene to the
change in phenotype, we propose a new metric: the ‘‘phenotypic impact factor (PIF).’’ PIF is a mathematical abstraction designed to ‘weight’ for the contribution the various DE genes make to the difference in the molecular anatomy of the two systems, based purely on their numerical properties.

The values were generated by combining the amount ofDE between the crosses, coupled with the average abundance calculated for both crosses at all time points for each ofthe 85 DE genes.

, the PIF metric is not particularly well suited
to regulators, although they were included in the analysis. Regulators are often stably expressed at close to baseline levels making detection of isolated changes in expression level challeng- ing and possibly misleading. To account for this, we ascribed ‘‘regulatory impact factors’’ (RIFs) to each of the 920 regulators based on their cumulative, simultaneous, DW to the DE genes, accounting for the PIF of the DE genes. This metric was intended as a mathematical abstraction to represent the relative importance of the regulators in driving the phenotypically relevant part of the network described above, based on differences in their correla- tions.

PIF=The average expression (state 1 and state 2 combined) multiplied by the DE (see above for definition), computed for all DE genes

RIF=The cumulative DW of each regulator relative to the target DE genes, weighted for PIF

DW=The difference in co-expression between a specified pair of genes in two different states.

##### Madhamshettiwar et al. 2012

comparative evaluation of nine state-of-the art gene regulatory network inference methods encompassing the main algorithmic approaches (mutual information, correlation, partial correlation, random forests, support vector machines) using 38 simulated datasets and empirical serous papillary ovarian adenocarcinoma expression-microarray data. We





##### Brunner et al. 2021

a novel **stochastic dynamical systems model that predicts gene expression levels from methylation data of genes in a given GRN**. 

https://github.com/kordk/stoch_epi_lib.

we developed a dynamic interaction network model [25] that depends on
epigenetic changes in a gene regulatory network (GRN). Dynamical systems integrate a set of simple interactions (i.e., transcription factor (TF) binding to a promoter region and subsequent gene expression) across time to produce a temporal simulation of a physical process (i.e., gene regulation in a given GRN). Therefore, the predictions of a dynamical systems model (e.g., TF binding and unbinding events, gene expression levels) emerge from a mechanistic understanding of a process rather than the associations between data (e.g., predicting an outcome from a set of predictor variables). A dynamical systems model can predict gene expression using epigenetic data and a GRN by simulating hypothesized mechanisms of transcriptional regulation. Such models provide predictions based directly on these biological hypotheses, and provide easy to interpret mechanistic explanations for their predictions. The dynamical systems approach offers a number of unique characteristics. First, a stochastic dynamical system provides us with a distribution of gene expression estimates, representing the possibilities that may occur within the cell. Next, the mechanistic nature of the approach means that the model can provide a biological explanation of its predictions in the form of a predicted activity level of various gene-gene regulatory interactions. Finally, a dynamical systems approach allows for the prediction of the effects of a change to the network. To our knowledge, there are no studies that have taken a dynamical systems approach to predicting gene expression from methylation data and a GRN.







##### Azad et al. 2021

Signalling transduction pathways (STPs) 

methodologies with **a fully Bayesian approach in discovering novel driver bio-markers in aberrant STPs given high-throughput gene expression (GE) data**.

‘PathTurbEr’ (Pathway Perturbation Driver) uses the GE dataset derived from the lapatinib (an EGFR/HER dual inhibitor) sensitive and resistant samples from breast cancer cell lines (SKBR3).

Differential expression analysis revealed 512 differentially expressed genes (DEGs) and their pathway enrichment revealed 13 highly perturbed singalling pathways in lapatinib resistance, including PI3K-AKT, Chemokine, Hippo and TGF-β singalling pathways.

the aberration in TGF-β STP was modelled as a **causal Bayesian network (BN)** using three MCMC sampling methods, i.e. Neighbourhood sampler (NS) and Hit-and-Run (HAR) sampler that potentially yield robust inference with lower chances of getting stuck at local optima and faster convergence compared to other state-of-art methods. Next,

, we examined the structural features of the optimal BN as a statistical process that generates the global structure using p1-model, a special class of Exponential Random Graph Models (ERGMs), and MCMC methods for their hyper-parameter sampling.

This step enabled key drivers identification that drive the aberration within the perturbed BN structure of STP

A BN is defined as a directed acyclic graph (DAG), ‘G’ = (V,E),
where ‘V’ is a set of random variables (here genes/proteins) and ‘E’ represents relationships among those random variables. Each random variable is associated with conditional probability distributions given its parents except the root variable, which are only associated with corresponding prior probability distri- butions. Causal BNs possess causal edges, for example X−→Y
indicates causality from random variable ‘X’ to random variable ‘Y’. An important property of BNs is called ‘Markov condition’, which must be satisfied for the probability distributions of the constituent random variables. ‘Markov condition’ states that each variable in the BN should be conditionally independent of other random variables (i.e. nodes) that are its non-descendants (i.e. not their children/grand-children) given its parents.



Bayesian network (BN) has been applied to model
aberrant STPs by studying case/control gene expression data [8] using MCMC methods such as Metropolis–Hastings or Gibbs sampling algorithm. Like other non-MCMC approaches (i.e. score-based methods such as Greedy, Hill-Climbing, Tabu search algorithm), these MCMC approaches also suffer local optimal problem, which means models claimed as optimal from those methods may not be globally optimal. Recently, we have shown that compared to Metropolis–Hastings algorithm, newly developed Neighbourhood sampler (NS) and Hit-and-Run sampler find BN structures with better accuracy and faster convergence towards true distribution



after modelling data-driven STP (optimal BN), statistical approaches relying on probabilistic models should be adopted to assess the probability of forming each interaction (i.e. aberrant activities among signalling proteins) within that BN structure. Being originally proposed in social network analyses, Exponential Random Graph Models (ERMGs), or p*-models [13], are central of statistical modelling of networks, where each possible edge in a network is considered as a random variable and modelled as combinations structural properties of the constituent nodes, e.g. global density of nodes, attractiveness/expansiveness of nodes (commonly known as sociality parameters). Note, non-statistical methods such as descriptive approaches may also model hub nodes based on their degree distributions which may not address stochastic nature of the underlying data and experimental measurements [14]. In the context of this study, hub nodes should be highly social, which would allow us to analyse the possibilities of detecting those hub nodes as a key bio-marker underlying the aberrant STP structure formation.

The full code base for PathTurbEr is available in here: https://github.co m/Akmazad/PathTurbEr/ 

modelled the signal transduction path- ways (STPs) as Bayesian Networks (BNs),



**a fully Bayesian statistical modelling approach for analysing the statistical aspect of the perturbed STP structure** yielded from the MCMC sampling algorithms of BN structure learning (see previous subsection). We have used p1-model, initially proposed by Holland and Leinhardt [18], is a special class of exponential families of distributions, for this study that offers robust and flexible parametric models, which are used to evaluate the probability that a gene to be hub in the perturbation network inferred in previous subsection

, the exponential family of distribution (i.e. p1-model) is a common choice for Pr(u) as it explicitly allows the distributions parameters to be tied with dif- ferent network statistics (e.g. global density of nodes, in-degree, out-degree, etc.) that control the formation of that network at the first place

we have employed a fully Bayesian approach to infer the posteriori ofparameters of the above p1-model. For that purpose, MCMC sampling methods i.e. Gibbs sampling approach were adopted, which

. Our approach relies on a hierarchical Bayesian model 



##### Raj-Kumar et al. 2019

**Principle Component Analysis-based iterative PAM50 subtyping (PCA-PAM50) to perform intrinsic subtyping in ER status unbalanced cohort**s

This method leverages **PCA and iterative PAM50 calls to derive the gene expression-based ER status and a subsequent ER-balanced subset for gene centering.**

Applying PCA-PAM50 to three different breast cancer study cohorts, we observed improved consistency (by 6–9.3%) between intrinsic and clinical subtyping for all three cohorts. Particularly, a more aggressive subset of luminal A (LA) tumors as evidenced by higher MKI67 gene expression and worse patient survival outcomes, were reclassified as luminal B (LB) increasing the LB subtype consistency with IHC by 25–49%.

ftp://ftp.wriwindber.org

PCA-PAM50 enhances the consistency of breast cancer intrinsic and clinical subtyping by reclassifying an aggressive subset of LA tumors into LB. PCA-PAM50

The originally defined normal-like1 (Normal) breast cancer subtype is now less frequently used6–8. Clinical subtyping of BC is based on immunohistochemistry (IHC) assays for the estrogen receptor (ER), progesterone receptor (PR), and human epidermal growth factor receptor 2 (HER2). More institutions now include Ki67, thus classifying tumors into triple-negative (TN; ER−/PR−/HER2−), HER2+ (ER−/PR−/HER2+), LA (ER+/HER2−/Ki67−), LB1 (ER+/HER2−/Ki67+) and LB2 (ER+/HER2+)9–12. For clarity, we will refer to clinical subtyping as IHC subtyping

IHC subtyping is the only accepted molecular assay for patient treatment decision-making13–17. BC intrinsic and clinical subtypes do not completely match even when comparing intra-subtypes, especially for LB18

The PAM505 classifier, which is also deployed in Genefu Bioconductor package20, makes calls based on the 50 gene centroid correlation distance to LA, LB, Basal, Her2 and normal-like centroids

the application of the PAM50 algorithm has its challenges. The two main challenges are (1) balancing ER status and

(2) the gene centering procedures5,21. The PAM50 classifier works accurately if the original cohort/dataset is ER status-balanced. However this is often not the case with most genome-wide studies. In such cases, a conventional strategy is to choose a subset which is ER status-balanced and use the median derived from that subset to gene center the entire cohort. In practice, an ER-balanced subset is chosen based on IHC-defined ER status. There have been reports that IHC-defined ER status, which is based on protein expression, not being completely consistent with ER status defined by gene expression22,23. This inconsistency may impact the accuracy of the subsequent gene centering procedure which is aimed to minimize the bias of the dynamic range of the expression profile per sequencing technology. As a result, such inconsistency may contribute to the discrepancy between the IHC and PAM50 subtyping results.

Hence, we explored the possibility of using a gene expression-based ER-balanced sub- set for gene centering leveraging principal component analyses (PCA) and iterative PAM50 calls to avoid intro- ducing protein expression-based data into a gene expression-based subtyping method. We validated our method termed PCA-PAM50 usig three different primary breast tumor datasets: an in-house 118 patients cohort, The Cancer Genome Atlas (TCGA)24 Breast cancer RNA-Seq 1097 patients cohort obtained from the Genomic Data Commons [https://gdc.cancer.gov/], and the Molecular Taxonomy of Breast Cancer International Consortium (METABRIC)25 discovery set of 997 patients. Our method resulted in improved consistency between PAM50 calls and IHC subtypes compared to the conventional method for all three cohorts, and this improved consistency is attributable to re-classification of an aggressive subset of LA tumors as LB

The PAM50 calls resulting from the use of the primary ER-balanced subset are called conventional intrinsic subtypes. The secondary ER-balanced set was based on principal com- ponent 1 (PC1) separation that agreed with IHC status: ER-negative

The PAM50 calls result- ing from the use of the secondary ER-balanced subset are called intermediate intrinsic subtype. The tertiary ER-balanced set was based on the intermediate intrinsic subtype’s Basal and LA calls

The PAM50 calls resulting from the use of the tertiary ER-balanced set are called refined intrinsic subtype



##### Tripathi et al. 2017

sgnesR (Stochastic Gene Network Expression Simulator in R) is an R package that provides an interface to simulate gene expression data from a given gene network using the stochastic simulation algorithm (SSA)



##### Zhang et al. 2020

we propose a new framework with a new metric to identify driver modules with low-frequency mutation genes, called **iCDModule. Inspired by the gravity model, we integrate the coverage and mutual exclusivity in mutation information, define a new metric between gene pairs, called mutation impact distance, to help identifying potential driver genes sets**, including those have extremely low mutation rates but play an important role in functional networks. A genetic network is constructed by combining the defined mutation impact distance and then the driver module identification problem is formalized as the maximum clique solution problem, and an improved ant colony optimization algorithm is used to solve it. iCDModule is applied to TCGA breast cancer, glioblastoma, ovarian cancer to test performance.

MDPFinder [9], Multi-dendrix [12], ComMDP and SpeMDP [13] used integer linear programming for settling the problem of maximum cover-exclusive sub-matrix to detect driver pathways.

Hotnet [16], Hotnet2 [17], Hierarchical Hotnet [18], use a thermal diffusion method on PPI network and the diffusion value is used to extract pathways or modules with high con- nectivity. MEMo [19] use interaction networks and functional relationship graphs to derive the largest clique in similar graphs, and combine the mutual exclusivity to process the largest clique. Babur et al. [20] introduce a method based on seed growth on the genetic network. It applies TCGA datasets to detect pan-cancer pathways or modules, and determines the growth strategy through a properly defined by mutual exclusivity score. BeWith [21] put the interaction density and mutual exclusivity as the optimization goal, using an ILP method to solve it. MEMCover [5] combines mutation data and interaction data to identify mutual exclusivity mutation gene in the same or different cancer tissues

we propose a new approach to de novo identify driver modules with low-frequency mutation genes, called iCDModule. The

The mutation score of each gene is calculated based on the mutation damage rate and coverage characteristics to measure the contribution of the mutation to cancer. 

Inspired by the gravity model [22], based on the mutation scoring function and mutation exclusivity, a new metric is defined between gene pairs, namely mutation impact distance. Construct a network of mutation genes with nodes corresponding to the mutation genes. If the mutation impact distance between genes is greater than the average of all non- zero mutation impact distances and they have edge in the PPI network, generate edge between them. The significance of the edges between gene pairs takes into account coverage, exclusivity, and functional correlation.

The problem of driver module identification is transformed into a maximum clique solution problem, and an improved ant colony optimization algorithm is used to settle this problem. Finally,





##### Gosline et al. 2012

**SAMNet, for Simultaneous Analysis of Multiple Networks, that is able to interpret diverse assays over multiple perturbations.** The algorithm uses a constrained optimization approach to integrate mRNA expression data with upstream genes, selecting edges in the protein–protein interaction network that best explain the changes across all perturbations

The result is a putative set of protein interactions that succinctly summarizes the results from all experiments, highlighting the network elements unique to each perturbation.

the human dataset measured cellular changes in four different lung cancer models of Epithelial-Mesenchymal Transition (EMT), a crucial process in tumor metastasis. SAMNet

SAMNet, for Simultaneous Analysis of Multiple Networks, an algorithm that uses a network flow model to integrate two distinct high-throughput experiments across multiple conditions.

SAMNet uses a constrained optimization formulation based
on the multiple commodity flow problem to model multiple experiments simultaneously as ‘‘commodities’’ that must transit from a common source to a common sink through a shared protein interaction network. Each edge in the interaction network has a particular capacity, and therefore must be ‘shared’ by all commodities. This constraint forces the algorithm to select interactions that are unique to each cellular perturba- tion, thus avoiding the selection of common stress pathways, a common pitfall of other optimization approaches.

SAMNet is a powerful tool for modeling diverse sources of high throughput data across multiple experiments







##### Gabr et al. 2013

, we consider **the problem of finding causal orderings of nodes in such protein interaction networks to discover signaling pathways. We adopt color coding technique to address this problem.** Color coding method may fail with some probability. 

. Our key contribution in this paper is elimination of the key conservative assumptions made by the traditional color coding methods while computing its success probability. We do this by carefully establishing the relationship between node colors, network topology and success probability. 





##### Erdogdu et al. 2017

Formulate the induction and control of gene regulatory networks (GRNs) from gene expression data using **Partially Observable Markov Decision Processes (POMDPs)**

partial observability would be a more natural and realistic method for handling the control of GRNs

We propose a method for the construction of POMDP model of GRN from only raw gene expression data which is original and novel. Then, we introduce a novel approach to decompose/factor the POMDP model into sub-POMDP’s in order to solve it efficiently with the help of divide-and-conquer strategy.

we focus on the GRN control problem. The problem
requires maintaining certain expression level for a single gene or certain expression levels for a group of genes

One of the important aspects of GRN control problem is that
it is not possible to obtain complete state information. 

n this paper, we are proposing to model the GRN control prob-
lem in a more natural and realistic way, mainly as a Partially Observable Markov Decision Process (POMDP). In fact, there are some aspects of the problem that are not fully observable, and unfortunately all the above mentioned models assume full observ- ability and hence simplify the problem. We argue that it is only possible to solve the GRN control problem in a realistic setting if partially observability is properly accounted in the model







##### Li et al. 2020

To identify driver modules with rarely mutated genes, we propose **a functional similarity index to quantify the functional relationship between rarely mutated genes and other ones in the same module**. Then, we develop a method to detect **Driver Modules with Rarely mutated Genes (DMRG)** by incorporating the functional similarity, coverage and mutual exclusivity

If a rarely mutated gene is in- cluded in a driver module, the functional relationships between this rare gene and other ones in the same module can be high. Thus, we propose a functional similarity index to quantify these relationships with network propagation, which can amplify weak similarities between different genes. This functional similarity can successfully capture functional relationship between driver genes with any mutation frequencies in the same driver module, especially for rarely mutated genes. Then, we develop a method to detect Driver Modules with Rarely mutated Genes (aka, DMRG) by incorporating the functional similarity, coverage and mutual exclusivity. 



##### Hu et al. 2012

network analysis to two independent human breast cancer data-sets and three different mouse populations developed for quanti-tative analysis of metastasis. Analysis of these datasets revealed
that the gene membership of the networks is highly conserved
within and between species, and that these networks predicted
distant metastasis free survival.

cross-species network analysis
of metastatic breast cancer. Using two publicly available human breast cancer gene-expression datasets that represent the natural progression of disease, as well as three experimental mouse populations used for identification of inherited metastasis sus- ceptibility genes, we have identified two gene networks that are independent predictors of metastatic disease in a meta-analysis of 1,881 human tumors (7). Unlike previously described human prognostic signatures, the networks significantly overlap between the two human datasets. In addition, significant overlap was also observed for the networks generated from the mouse samples

To determine which of the networks were associated with dis-
ease progression, Kaplan-Meier analysis was performed using the gene expression-based outcome (GOBO) algorithm

number of genes of a given network shared across species

##### Ihmels et al. 2005

the ‘‘**differential clustering algorithm’**’ for revealing conserved and diverged co-expression patterns. 

Existing approaches for comparative gene expression analyses emphasize mostly conserved co-regulation patterns, rather than differences in expression patterns [8,9,11].

To better capture differential expression patterns, we developed a novel approach, termed the **differential clustering algorithm (DCA), for systematically characterizing both similarities and differences in the fine structure of co-regulation patterns**

The DCA is applied to a set of orthologous genes that are
present in both organisms. As a first step, the pair-wise correlations between these genes are measured in each organism separately, defining two pair-wise correlation matrices (PCMs) of the same dimension (i.e., the number of orthologous genes) (Figure 2A). Next, the PCM of the primary (‘‘reference’’) organism is clustered, assigning genes into subsets that are co-expressed in this organism, but not necessarily in the second (‘‘target’’) organism. Finally, the genes within each co-expressed subgroup are re-ordered, by clustering according to the PCM of the target organism. This procedure is performed twice, reciprocally, such that each PCM is used once for the primary and once for secondary clustering, yielding two distinct orderings of the genes.

The results of the DCA are presented in terms of the rearranged PCMs. Since these matrices are symmetric and refer to the same set of orthologous genes, they can be combined into a single matrix without losing information. Specifically, we join the two PCMs into one composite matrix such that the lower-left triangle depicts the pair-wise correla- tions in the reference organism, while the upper-right triangle depicts the correlations in the target organism (

An automatic scoring method is then applied to classify clusters into one ofthe four conservation categories: full, partial, split, or no conservation of co-expression

##### Kao et al. 2004

For complex transcriptional networks, more sophisticated tools are required to **deconvolute the contribution of each regulator**. Here, we demonstrate the utility of **network component analysis in determining multiple transcription factor activities based on transcriptome profiles and available connectivity information regarding network connectiv- ity**

NCA takes advantage of the connectivity information to decompose DNA microarray data to determine both TFA and the control strength (CS) of each regulatory pair

everal linear decompositions of the matrix log [Er] have been used in the study of gene expression array, such as singular value decomposition (2) and independent component analysis (3). Although these decomposition tech- niques have strong statistical foundations, their molecular basis is difficult to pinpoint. The solution obtained by NCA is based not on any hypothesis
of relationship between the TFAs, but on the structure of [CS], namely, the connectivity structure of the network linking TFs and genes. Specifically, such constraints involve, for example, setting to zero the elements CSij when gene i is not regulated by TFAj, but can also include constraints on the polarity of the regulation (induction or repression). We demonstrated (1) that if the underlying transcriptional network satisfies the following properties, such decomposition becomes unique up to some normalization factors: (i) The connectivity matrix [CS] must have full-column rank. (ii) When a node in the regulatory layer is removed along with
all of the output nodes connected to it, the resulting network must be characterized by a connectivity matrix that still has full-column rank.
(iii) The log [TFAr] matrix must have full row rank. In other
words, each regulatory signal cannot be expressed as a linear combination of the other regulatory signals. This criterion requires M ⬎ L as a necessary but not sufficient condition







##### Ma et al. 2016

several methods for **active learning of causal networks** has been developed recently14–18. The active learning methods utilize both observational and experimental data to discover causal networks. These methods typically first construct a draft of the causal network, generally represented as an unori- ented or partially oriented graph, from observational data. Then, the methods select a variable for experimen- tation/manipulation to further refine the graph. The experimental data obtained from the targeted experiment is used to update the draft of the causal network. The process of variable selection, experimentation, and causal network update is repeated until some termination criterion is satisfied, e.g. all edges in the causal network are oriented. Since randomized controlled experiments are costly, active learning methods employ various heuris- tics when selecting variables for experimentation in order to minimize the required number of experiments. It

Among the active learning methods examined in this study, ODLP variants achieved the best local pathway
reconstruction quality with low cost on the 5 transcription factors examined. In



##### Maere et al.  2008

the distance measures used in traditional clustering algorithms have difficulties in detecting one of the most prominent features of perturbational data, namely partial correlations between expression profiles. Biclustering methods on the other hand are specifically designed to capture such partial correlations. However, most biclustering algorithms do not provide measures for pair-wise expression correlation between genes, but rely on emergent properties of groups of genes and conditions (modules) in order to identify statistically significant subpatterns in the data. This reliance complicates the elucidation of less modular regions in the underlying transcriptional network.

gene expression profiles are discretized into three categories (upregulated, downregulated, unchanged) based on p-values for differential ex- pression. For each pair of profiles, we then assess the probability that the observed overlap of upregulated and downregulated fields is generated by chance. The resulting correlation p-values are corrected for multiple testing and translated to edges in a coexpression net- work, which is then clustered into (overlapping) expression modules using a graph clustering procedure that identifies densely connected components in the network. Relevant condition sets are then determined for all modules and the modules are screened for enrichment of Gene Ontology categories and transcription factor binding sites. Finally, a regulation pro- gram is learned for each module in an attempt to explain the expression behavior of the module’s genes as a function of the expression of a limited set of regulators (transcription factors and signal transducers). 

ENIGMA, that addresses some of these issues. **ENIGMA leverages differential expression analysis results to extract expression modules from perturbational gene expression data.** 



Our goal was to build a method that: (i) leverages differential expression analy- sis results to extract co-differential expression networks and expression modules from perturbational gene expres- sion data, (ii) is able to detect significant partial coexpression relationships between genes and overlap between modules, (iii) depends on parameters that can be auto- matically optimized or set on reasonably objective grounds. (iv) produces a realistic amount of modules, and (v) visually integrates the expression modules with other data types such as Gene Ontology (GO) information [28], transcription factor (TF) binding data, protein and genetic interactions, in order to facilitate the biological interpreta- tion of the results. 

ENIGMA takes as input a set of perturbational expression data, externally calculated p-values for differ- ential expression (e.g. using the limma package in Biocon- ductor [29]) and other data types if available. ENIGMA uses a novel combinatorial statistic to assess which pairs
of genes are significantly co-differentially expressed (henceforth abbreviated as coexpressed for the purpose of readability). The resulting coexpression p-values are cor- rected for multiple testing and translated to edges in a coexpression network, which is clustered into expression modules (i.e. groups of significantly co-differentially expressed genes) using a graph-based clustering algorithm inspired on the MCODE algorithm [30]. The clustering procedure depends on two parameters that control the density of individual modules and the overlap between modules. The main reason why we chose a two-tier clus- tering approach (data → coexpression network → cluster- ing) is that it allows simulated annealing-based optimization of the clustering parameters to obtain opti- mal coverage of the coexpression network, in terms of module overlap and redundancy. The

In the post- processing phase, ENIGMA determines relevant condition sets for each module, visualizes their substructure and overlap with other modules, screens the modules for enriched GO categories, suggests potential regulators for the modules based on regulator-module coexpression links and enrichment of TF binding sites, and overlays protein and genetic interaction data.





##### Mall et al. 2017

The ability to detect statistical relevant changes in the interaction patterns induced by the progression of the disease can lead to the discovery of novel relevant signatures. Several procedures have been recently proposed to detect sub-network differences in pairwise labeled weighted networks. Methods: In this paper, we propose an **improvement over the state-of-the-art based on the Generalized Hamming Distance adopted for evaluating the topological difference between two networks and estimating its statistical significance.** The proposed procedure exploits a more effective model selection criteria to generate p-values for statistical significance and is more efficient in terms of computational time and prediction accuracy than literature methods. Moreover, the structure of the proposed algorithm allows for a faster parallelized implementation.







##### Wouters et al. 2019

cancer cells from each sample form a distinct cluster per patient, whereas the corresponding normal host cells from various patients cluster together according to their cell type (Tirosh et al. 2016; Puram et al. 2017; Lambrechts et al. 2018). This observation is somewhat counterintuitive because cells with similar gene expression profiles are known to occur in multiple tumors, for instance cells in specific cell cycle stages. 

gene regulatory network inference using SCENIC has been shown to normalize away part of these tumor- specific differences, resulting in one pan-tumor cluster of cycling cells (Aibar et al. 2017). Nonetheless, the unsupervised discovery of common transcriptional states remains a challenge

**SCENIC network inference to the single-cell expression matrix**

**. A transcription factor with its candidate targets is called a regulon. SCENIC yields a regulon-cell matrix with regulon activities across all single cells, and provides therefore an alternative dimensionality reduction**. A UMAP visualization based on the regulon-cell matrix reveals three candidate cell states in an unsupervised manner,



##### Xie et al. 2020

methods are urgently needed which can separate the impact of true regulatory elements from stochastic changes and downstream effects. We propose **the differential network flow (DNF) method to identify key regulators of progression in development or disease.** **Given the network representation of consecutive biological states, DNF quantifies the essentiality of each node by differ- ences in the distribution of network flow, which are capable of capturing comprehensive topological dif- ferences from local to global feature domains.** 

s a new approach for quantifying the essentiality of genes across networks of different biological states

The existing differential network analy- sis methods mainly fall into two categories. The first category is focused on capturing linear or nonlinear correlation differences in gene expression between two gene regulatory networks (GRNs). For instance, DDN [10] is the first algorithm to detect topological differences by lasso regression in network inference. DISCERN computes a novel perturbation score to capture how likely a given gene has a distinct set of regulators between different condi- tions, which is shown to be robust to errors in network structure estimation. pDNA [12] incorporates prior information into differ- ential network analysis using non-paranormal graphical models, which relaxes the assumption of normality of omics data to find more cancer-related genes



The existing differential network analy- sis methods mainly fall into two categories. The first category is focused on capturing linear or nonlinear correlation differences in gene expression between two gene regulatory networks (GRNs). For instance, DDN [10] is the first algorithm to detect topological differences by lasso regression in network inference. DISCERN computes a novel perturbation score to capture how likely a given gene has a distinct set of regulators between different condi- tions, which is shown to be robust to errors in network structure estimation. pDNA [12] incorporates prior information into differ- ential network analysis using non-paranormal graphical models, which relaxes the assumption of normality of omics data to find more cancer-related genes

The second category is focused on topological differences
between constructed GRNs. For instance, DEC captures the global differential eigenvector connectivity to prioritize nodes in net- works [7]. DiffRank [13] computes the linear combination of differ- ential connectivity and differential betweenness centrality to order genes. DCloc [14] computes the average proportion of changes of each node’s neighborhood as a significance score by iteratively removing edges with different thresholds. DiffNet [15] evaluates topological differences between two networks based on general- ized hamming distance and its statistical significance. TKDS [8] measures the importance of genes by calculating the graphlet vec- tor distance

Two issues remain which the above methods fail to address.
First, all the above methods utilize networks which assume the existence of edges based on co-expression 

The second issue is more subtle. While techniques such as spec-
tral analysis provide a global perspective on connectivity, these approaches fail to encapsulate the flow of information inherent in all biological networks. Network

we propose **the differential network flow (DNF)**
**method to identify key regulators between two networks under different biological conditions. This algorithm is built upon the idea of network flow and information theory.** Rewiring of a GRN can be characterized as a dynamic pattern of network flow [19], such a flow-based model captures multiple (from local to global) features of network structure. Information theory is able to quantify the uncertainty in networks, making networks built upon information-theoretic measurements a more acceptable representation of biological systems at the molecular scale [20]. Therefore, DNF is capable of capturing comprehensive topo- logical differences by quantifying the flow in a network.

DNF is built upon the ideas of network flow and information
theory. The novelty of DNF lies in quantifying node-to-node infor- mation entropy according to the network flow in a gene regulatory network, and in characterizing each node as a distribution of net- work flow, which is equal to the distribution of information entropy. The distribution differences of one gene in different net- works represents its essentiality in the biological process responsi- ble for the network’s evolution. Genes are ordered by the magnitude of this difference to establish a ranking

we
**employ a three-step process to integrate both transcriptomics and proteomics datasets in GRN-construction.** First, a network skeleton is built by differential expression analysis using the tran- scriptomics dataset. Specifically, the skeleton gene sets are selected based on a given criterion, such asjlog2FoldChangej > u, p-value < v, where u represents the fold change of gene expression and v rep- resents the statistical significance of differential expression. Sec- ond, the known corresponding protein–protein interactions in the STRING database (http://string-db.org) are used to establish the gene-gene network for the selected genes. For example, sup- pose p skeleton genes are selected in the first step, then a network skeleton with p nodes is described by an adjacency matrix Ai;j, such that Ai;j >0, i, j =1, ..., p, if protein i and protein j are functionally associated. Finally, the absolute value of the spearman correlation coefficient (scc) of expression is adopted to estimate the strength of connections between adjacent genes, and edges Ai;jfor which scc(i, j) < 0.1 are discarded (see

we construct a pair of GRNs based on a specific transcriptomics dataset (e.g. cancer and control samples) and a generic proteomics dataset







##### Van de Sande et al. 2020

SCENIC reconstructs regulons (i.e., transcription factors and their target genes) assesses the activity of these discovered regulons in individual cells and uses these cellular activity patterns to find meaningful clusters of cells. Here

. First, coexpression modules are inferred using a regression per-target approach (GRNBoost2). Next, the indirect targets are pruned from these modules using cis- regulatory motif discovery (cisTarget). Lastly, the activity of these regulons is quantified via an enrichment score for the regulon’s target genes (AUCell). 

t**he SCENIC pipeline consists of three steps. First, candidate regulatory modules are inferred from coexpression patterns between genes (Steps 5 and 6). Next, coexpression modules are refined by the elimination of indirect targets using TF motif information (Step 6). Finally, the activity of these discovered regulons is measured in each individual cell and used for clustering (**Steps 7 and 8; Fig. 1)

<u>Network inference</u> (Step 5) In a first step, given a predefined list of TFs, regulatory interactions between these factors and putative target genes are inferred via regression-based network inference10 from the expression or count matrix. More specifically, for every expressed gene, a regression model is built that predicts the expression of this gene across cells (the response) from the expression of these predefined TFs (the independent variables or regressors). A regulatory interaction between gene and TF is inferred by assessing the importance of these factors in the regression model

Our workflow relies on a tree-based regression model and, by default, uses an efficient and distributed implementation based on gradient boosting machine regression (i.e., GRNBoost2 (ref. 11)). This algorithm is a reimplementation of the GENIE3 algorithm, which was based on a random forest1

The output of this step is a list of adjacencies connecting a TF with a target gene. A weight or importance is associated with these connections to distinguish strong from weak regulatory interactions

<u>Module generation</u> (Step 6) From these regulatory interactions, inferred from coexpression patterns, modules (i.e., a TF and its predicted target genes) are generated

Multiple strategies are combined so that, for every factor, multiple modules of different sizes are created

Subsequently, modules are differentiated between transcriptional activation and repression based
on the correlation patterns between the expression of the regulator and its targets: t

, the resulting modules with <20 genes are not retained, as these modules tend to be less
stable and more sensitive to target gene dropout.

<u>Motif enrichment and TF-regulon prediction</u> (cisTarget step; Step 6) Modules contain direct and indirect targets of a regulator because these regulatory interactions are only inferred based on coexpression patterns. Therefore, in this step, the putative regulatory regions of these target genes are searched for enriched motifs by comparing scores of cis-regulatory modules (CRMs) near the genes in the module, with the remaining genes in the genome. CRM scoring is done using hidden Markov models, with a large collection of position weight matrices (PWMs) following the procedures described in iRegulon13 and i-cisTarget14,15. If the motif of the regulator TF is significantly enriched in one of its modules, this regulator and its predicted targets are retained for further analysis.

<u>Cellular enrichment</u> (AUCell step; Step 7) A similar ranking and recovery framework is used to quantify the activity of the predicted regulons in the individual cells that make up the scRNA-seq experiment. In detail, each individual cell’s tran- scriptome is modeled as a whole-genome ranking based on the expression of its genes. The enrichment of a regulon is subsequently assessed via recovery of its targetome on the cell’s whole- genome ranking. The AUC metric measures the relative biological activity of a regulon in a given cell

(Optional) <u>Binarization of cellular regulon activity</u>

<u>Clustering of cells based on regulon activity</u> (Step 8) Quantification of regulon activity via AUCell is a biological dimensional reduction: the number of discovered regulons (k) is typically much lower than the number of genes (n), and therefore each cell can be represented by a k-dimensional vector instead of being represented by a point in n-dimensional space.

This eliminates the need to reduce the dimensionality via principal component analysis (PCA)
before applying a nonlinear projection technique such as t-distributed stochastic neighbor embedding (t-SNE) or uniform approximation and projection (UMAP) for rendering visual groupings16.In addition, this regulon-based clustering also reveals the GRNs underlying the gene expression profiles



##### Wu et al. 2021

In this study, we propose a new recognition method of driver modules, named **ECSWalk to solve the issue of mutated gene heterogeneity and improve the accuracy of driver modules detection, based on human protein–protein interaction networks and pan-cancer somatic mutation data**. This study first utilizes high mutual exclusivity and high coverage between mutation genes and topological structure similarity of the nodes in complex networks to calculate interaction weights between genes. Second, the method of random walk with restart is utilized to construct a weighted directed network, and the strong connectivity principle of the directed graph is utilized to create the initial candidate modules with a certain number of genes. Finally, the large modules in the candidate modules are split using induced subgraph method, and the small modules are expanded using a greedy strategy to obtain the optimal driver modules. 

although the previous methods can
detect the gene sets with high mutual exclusivity and high coverage, they just focus on the mutual exclusivity and cov- erage between genes, instead of the topological structure of complex networks. To effectively solve the problem of mutated gene heterogeneity and improve the accuracy of driver modules, this study proposes a driver module detec- tion algorithm (ECSWalk) based on gene mutation and human protein–protein interaction network. The algorithm takes into account aspects, such as high mutual exclusivity and high coverage between genes, and high similarity of topological structure. First,

First, the complex network topology analysis method is used in human protein–protein interaction network data to calculate the topological similarity between network nodes, and then the two characteristics of high cov- erage and high mutual exclusivity of the mutated genes are combined to obtain the weight of the vertices and edges in the human protein–protein network. The weights of vertices in the human protein–protein network are obtained accord- ing to the coverage of mutated gene, and the random walk with restart strategy is utilized to calculate the weights of edges in the network by the three characteristics, namely, the coverage, the mutual exclusivity, and the similarity of the topological structure between the nodes

Second, based on the weighted network constructed in the previous step, the large modules are split into several candidate gene sets using the method of the induced subgraph. In addition, the greedy strategy is utilized to add the nodes in the leaf module to the seed module to achieve the optimal gene sets. These mutated gene sets with high mutual exclusivity, high coverage and high similarity of the topological structure are likely to work as driver modules in cancer



##### Wang et al. 2019

we explore **evolutionary algorithms**, and their applications with sparse matrix representations. Our approach can speed up the optimization process and find good solutions, uncovering the underlying GRNs

##### Wu et al. 2021

**ShareNet, a Bayesian framework for boosting the accuracy of cell type-specific gene regulatory networks by propagating information across related cell** types via an information sharing structure that is adaptively optimized for a given single-cell dataset.

The techniques we introduce can be used with a range of general network inference algorithms to enhance the output for each cell type. We

We introduce ShareNet, a Bayesian information sharing frame-
work for increasing the accuracy of predicting cell type-specific regulatory associations from single-cell transcriptomic data (Fig. 1b). Our framework draws upon the intuition that many of the regu- latory interactions (and non-interactions) are shared across different cell types, due to shared developmental lineages, regulatory pro- grams, or biophysical constraints. Thus, by propagating information across related cell types, we hope to reduce noise and boost the ac- curacy of inferred GRNs in all cell types. Since we do not have full knowledge of the sharing patterns underlying a given dataset, we designed our framework to adaptively learn a multifactorial, infor- mation sharing structure that best explains the data in all the study’s cell types. Importantly, our framework is widely applicable, as it can serve as an additional layer on top of existing state-of-the-art network inference algorithms to enhance their accuracy in estimat- ing the GRNs of all cell types in a dataset. 





##### Parikh et al. 2010 

we describe **a Bayesian network approach that addresses a specific network within a large dataset to discover new components**. Our algorithm draws individual genes from a large gene-expression repository, and ranks them as potential members of a known pathway. We

Information theory approaches, such as ARACNE, compare expres- sion profiles between all genes using mutual information as a generalized measure of correlation

Bayesian networks are useful because they can model higher than pairwise orders of dependences between genes and can incorporate existing knowledge 

We used Bayesian networks [20] to model the core PKA pathway (Figure



##### Suo et al. 2015

An analysis pipeline is built for **integrating genomic and transcriptomic alterations from whole-exome and RNA sequence data and functional data** from protein function prediction and gene interaction networks.

The method accumulates evidence for the functional implications of mutated potential driver genes found within and across patients. A driver-gene score (DGscore) is developed to capture the cumulative effect of such genes.

To contribute to the score, a gene has to be frequently mutated, with high or moderate mutational impact at protein level, exhibiting an extreme expression and functionally linked to many differentially expressed neighbors in the func- tional gene network. The

these methods do not utilize iso- form-level information and the potential drivers are generally not validated in terms of patients’ clinical outcomes such as survival

we summarize the effects of potential driver genes into a single value DGscore and assess its clinical value as prognostic biomarker. In



##### Pillai et al. 2021

Various markers or regulators associated with distinct phenotypes in melanoma have been identified, but, how does a network of interactions among these regulators give rise to multiple “attractor” states and phenotypic switching remains elusive. 

we inferred a network of transcription factors (TFs) that act as master regulators for gene signatures of diverse cell-states in melanoma. Dynamical simulations of this network predicted how this network can settle into different “attractors” (TF expression patterns), suggesting that TF network dynamics drives the emergence of phenotypic heterogeneity.

**To identify the master regulators for the differentially expressed genes obtained from WGCNA, we used geWorkbench (Floratos et al., 2010). At first, we identified a baseline transcriptional interaction network for the dataset, using ARACNE (Algorithm for the Reconstruction of Accurate Cellular Networks) (Margolin et al., 2006).** A p-value of 10-7 was set to determine the mutual information threshold and the software was run for 100 bootstraps with data processing inequality set to 0. Fisher’s exact test was used to identify master regulators from a list of candidate master regulators (Lambert et al., 2018). Only those transcription factors (TFs) enriched for in the WGCNA modules with p-value < 0.05 were considered for further analysis. This list was cross validated against CHEA, ENCODE and ARCHS4 databases by using EnrichR (Chen et al., 2013) to identify potential TFs regulating each module. Only those TFs identified by both analyses (ARACNE and EnrichR) were considered as master regulators (Table

RAndom CIrcuit PErturbation (RACIPE) (Huang et al., 2017) was used to generate an ensemble of ordinary differential equation (ODE) models. Each model represents a collection of modified Hills equations for each gene, with randomized kinetics parameters sampled from user-defined ranges.

##### Ronellenfitsch et al. 2017

Complimentary to perturbation approaches, we extract functionally related groups of genes by analyzing the standing variation within a sampled population. To distinguish bi- ologically meaningful patterns from uninterpretable noise, we focus on correlated variation and develop a novel **density-based clustering approach that takes advantage of a percolation transition generically arising in random, uncorrelated data.**

we address the complementary challenge of **identifying the underlying regulatory relationships among genes from the standing variation in expression across sampled individuals.** **Rather than seeking to fully infer the underlying gene regulatory network topology from this inherently (and often prohibitively) noisy class of data (11, 14), we focus on identifying functional modules — sets of genes that demon- strate significant evidence for co-regulation.** Extracting gene modules from standing variation can be addressed by clustering expression patterns across samples, and has been attempted in the past with varying degrees of success (13, 15–20). Yet a primary challenge remains to distinguish true regulatory relationships from noise, and these efforts have depended on expert insights about the specific biological systems to appro- priately pre-filter genes, tune analysis parameters, and filter results. We have developed a novel, **data-driven approach motivated**
**by the theory of percolation on random graphs** (21–24). The method is conceptually simple yet robustly applicable, reliably yielding interpretable gene clusters across diverse data sets without fine-tuned optimization and filtering steps. We exploit the generic behavior of random geometric networks close to the percolation critical point, from which we devise a null model for the noise. This noise model in turn provides a basis for identifying statistically significant branches within the cluster hierarchy

.We leverage the standing variation across unperturbed samples to reveal functional modules in gene regulatory net- works (Fig. 1A). Groups of functionally related genes are expected to share a common pattern of expression variation across samples, the similarity of which can be quantified by a correlation-based distance measure



##### McClure et al. 2019

 six additional mutual information methods in the MINET R package (ARACNE, CLR, MIM, MINET, MRNET, MRNETB) [[43](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007241#pcbi.1007241.ref043)







##### Tapia-Carrillo et al. 2019

**an extension of the original Master Regulator Inference Algorithm (MARINa) analysis. This modified version of MARINa utilizes a restricted molecular signature containing genes from the 25 human pathways in KEGG's signal transduction category.**

TMRs were inferred using the MARINa (Lefebvre et al., 2010). **MARINa identifies TMRs through an enrichment of TF regulons (a TF with its targets) with differentially expressed genes between the two phenotypes (breast cancer vs. adjacent healthy mammary tissue). TMR inference with MARINa requires as input a network of regulons, a gene expression, molecular signature, and a null model** (Lefebvre et al., 2010) (Figure

To obtain a regulon set from the data, we used the expression
matrix of the tumor samples and a list of transcription factors in the TFCheckpoint curated database (Tripathi et al., 2013

As a first step, transcription factors are associated with other
genes expressed in the tissue. We used the mutual information- based algorithm ARACNe
(Margolin et al., 2006) which
calculates the pairwise mutual information for a pair of genes using the empirical probability distributions of their expression levels. For this network all possible interactions between TFs and genes in the expression matrix were calculated and kept if itspvalue was below 0.005

Mutual information can detect both indirect and direct
relationships. ARACNe constrains the number of indirect interactions applying the data processing inequality theorem (DPI), which considers that, in a triangle of interactions, the weakest one has a greater probability of being indirect if its difference is large with respect to the other two interactions (Hernández-Lemus

The
type of association (activation or repression) of the transcription factors is determined from the Spearman correlation of the TF with the levels of expression of all its targets (Lefebvre et al., 2010). This calculation was performed by the **aracne2regulon function in the viper R package (Alvarez et al., 2016). This function transforms the undirected MI network from ARACNE into a regulons network that is directed.**

**in the standard MARINa workflow, the molecular signature is built by comparing the expression level distributions of all genes between two conditions (e.g., healthy and diseased). For this work we built a molecular signature using only those genes annotated within the signal transduction pathways category in the KEGG database** (Kanehisa



To estimate the probability that a gene set enrichment score depends on the biological context and thus is not merely random, a null model was generated by random permutation between cases and control samples and recalculating differential expression values

**With the molecular signature, the regulon network and the null model, MARINa estimated the top regulons that enrich the most differentially expressed genes in the molecular signature through a gene set enrichment analysis**

The difference of this work with respect to MARINa lies in
the use of a specific set of genes (signal transduction signature) which produces a ranking of the regulons for this specific subset. It is important to note that the regulons network used as input is the same as in regular MARINa, but the ranking is focused on the specific gene signature. The set of genes that constitute each regulon may include genes that are not in the molecular signature and can be part of a more diverse range of biological functions than signal transduction. This is the reason why we performed a subsequent enrichment analysis of the regulons with all KEGG human pathways



##### Grechkin et al. 2016

DISCERN takes two expression datasets as input: an expression dataset of diseased tis- sues from patients with a disease of interest and another expression dataset from matching normal tissues. **DISCERN estimates the extent to which each gene is perturbed—having distinct regulator connectivity in the inferred gene-regulator dependencies between the dis- ease and normal conditions**. This approach has distinct advantages over existing methods. First, DISCERN infers conditional dependencies between candidate regulators and genes, where conditional dependence relationships discriminate the evidence for direct interac- tions from indirect interactions more precisely than pairwise correlation. Second, DISCERN uses a new likelihood-based scoring function to alleviate concerns about accuracy of the specific edges inferred in a particular network.

Most analysis methods that compare gene expression datasets from two conditions address the question ofwhich genes are significantly differentially expressed between conditions. **The DISCERN method addresses a distinct question concerning which genes are significantly rewired in the inferred gene-regulator network in disease tissues**








##### Jung 2019

We propose a method of Knowledge-based Evaluation of Dependency DifferentialitY (KEDDY), which is a statistical test for differential functional protein networks of a set of genes be- tween two conditions with utilizing known functional protein–protein interaction information. Unlike other approaches focused on differential expressions of individual genes or differentiality of individual interactions, **KEDDY compares two conditions by evaluating the probability distributions of functional protein networks based on known functional protein–protein interactions**



##### Gundogdu et al. 

the deep neural networks constrained by several types of prior biological information, e.g. signaling pathway information, as a way to reduce the dimensionality of the scRNA-seq data.

**including in the DNN architecture pathway knowledge allows obtaining a smaller architecture** (less nodes and hence faster inference), which is easier to interpret [22] and that performs as well as other methodologies in a set of cell type identi?cation benchmarks

**Prior Biological Information Integration**
**In order to incorporate the biological priors, the ?rst hidden layer was adjusted in two ways: 1) each neuron/node corresponds to one biological unit, in this case there are as many neurons as pathways and 2) the weights that arrive to a neuron are ?xed to zero when no input gene participates in the pathway associated to the node. In this way, biological priors were incorporated using known gene clusters with de?ned functions (the pathways) at the same time that the size of the model is reduced, which can help with over-?tting as well as training and inference time**



##### Giorgi et al. 2014

**CINDy (Conditional Inference of Network Dynamics), a novel algorithm for the genome-wide, context specific inference of regulatory dependencies between signaling protein and transcription factor activity, from gene expression data. T**he algorithm uses a novel adaptive partitioning methodology to accurately estimate the full Condition Mutual Information (CMI) between a transcription factor and its targets, given the expression of a signaling protein



##### Morone et al. 2015

we map the problem onto **optimal percolation in random networks to identify the minimal set of influencers**, which arises by minimizing the energy ofa many-body system, where the form of the interactions is fixed by the non- backtracking matrix15
ofthe network.



##### Zhang et al. 2017

We propose a **new differential network analysis method to address the above challenges. Instead of using Gaussian graphical models, we employ a non-paranormal graphical model that can relax the normality assumption. We develop a principled model to take into account the following prior information: (i) a differential edge less likely exists between two genes that do not participate together in the same pathway; (ii) changes in the networks are driven by certain regula- tor genes that are perturbed across different cellular states and (iii) the differential networks esti- mated from multi-view gene expression data likely share common structures**.



##### Ben Guebila et al. 2021

 GRAND (https://grand.networkmedicine.org) as a database for computationally-inferred, context-specific gene  regulatory network models that can be compared between biological  states, or used to predict which drugs produce changes in regulatory  network structure. The database includes 12 468 genome-scale networks  covering 36 human tissues, 28 cancers, 1378 unperturbed cell lines, as  well as 173 013 TF and gene targeting scores for 2858 small  molecule-induced cell line perturbation paired with phenotypic  information. GRAND allows the networks to be queried using phenotypic  information and visualized using a variety of interactive tools.



##### Gamez et al. 2015

Statistical analyses were conducted to associate miRNAs and
proteins. As a first approach to describe associations present in our database, we choose **probabilistic graphical models** com- patible with high dimensionality. The result is an undirected graphical model with local minimum Bayesian Information Criterion (BIC; ref. 14) obtained after executing the next steps: first the spanning tree with maximum likelihood is found and then (I), a forward search is performed by successively adding edges that reduce the BIC and still preserve the decomposability (15) of the initial graph (II). In the first stage, in order to learn a Markov tree structure from a random sample of a supposed multidimensional normal population, we used the extension of the Chow–Liu solution (16), according to which, for categorical data, the maximum likelihood structure is given by the max- imum weight spanning tree (17) with empirical mutual infor- mation quantities (18) as edge weights. In the Gaussian case, a similar reasoning applies, but now the mutual information value is ?(1/2)log(1 ? r2), where r is the empirical correlation coefficient between the two variables (nodes) joined by the edge. Given that the algorithm is invariant under monotone transformations of the variables, r2 can be used as a weight. In the second phase, we introduce the BIC criterion that penalizes more complex models and then, simpler graphs are generated. This is a fundamental objective in high-dimensional problems. Both methods are implemented in the open-source statistical programming language R (19). In particular, the functions minForest and stepw,inthe gRapHD package (20), are used for phases I and II, respectively

 Functional node identification and activity measurement To identify functional nodes within the probabilistic graphical
models, we split it in several branches. Then, we used **Gene Ontology analyses to investigate which function or functions were over-represented in each branch.** To measure the functional activity ofeach node, we calculated the mean expression ofall the proteins included in one branch related with a concrete function. Differences in node activity between ERþ and TNBC samples were assessed by class comparison analyses

One obvious limitation of proteomics when compared with genomics is that genomics can measure the expression of all known genes in the same experiment, whereas proteomics only provides a measurement of peptides that are both detected and identified. This

##### Zhan et al. 2019

we identified 73 breast cancer pa- tients from the TCGA and CPTAC projects matched whole slide images, RNA-seq, and proteomic data. By calculating 100 different morphological features and correlating them with the transcriptomic and proteomic data, we inferred four major biological processes asso- ciated with various interpretable morphological fea- tures. These processes include metabolism, cell cycle, immune response, and extracellular matrix development, which are all hallmarks of cancers and the associated mor- phological features are related to area, density, and shapes of epithelial cells, fibroblasts, and lymphocytes. In addition, protein specific biological processes were inferred solely from proteomic data, suggesting the importance of pro- teomic data in obtaining a holistic understanding of the molecular basis for tumor tissue morphology. 



TCGA aggregates an extensive col- lection of omics and clinical datasets from large cohorts of patients for more than 30 types of cancers (24). It also ar- chives histopathology images for solid tumor samples from which omics data were sampled. Currently, more than 24,000 histopathology images are available and can be visualized at the Cancer Digital Slide Archive (CDSA, http://cancer. digitalslidearchive.net/). In addition, The NCI Clinical Proteomic Tumor Analysis Consortium (CPTAC) (https://proteomics. cancer.gov/programs/cptac) program also provides high- throughput proteomic data for some of the TCGA tumor spec- imens such as breast cancer, ovarian cancer, and colorectal cancer based on mass-spectrometry technology. These



Based on this pipeline we performed a series of analysis correlating and integrating molecular data, morphological features, and clinical outcome using data from TCGA and CPTAC Breast invasive carcinoma (BRCA) project.

First, we performed a correlative analysis between multi-omics data (including pro- teomic, transcriptomics data) and morphological features ex- tracted from histopathology images. We observed that pro- teomic and transcriptomic data shared consistent correlation pattern with various morphological features at genome scale. However, comparing to transcriptomics data, proteomic data can identify specific protein-related biological processes as- sociated with morphological features that otherwise cannot be inferred from transcriptomic data. More comprehensive analysis revealed that four major categories of biology pro- cesses related to the hallmarks of cancer (6) are associated with different morphological feature based on the correlated proteomic data. Furthermore, we examined the relationship between nuclear morphology and patient outcome (i.e. sur- vival time). Both prognostically favorable and unfavorable morphological features have been identified. The biological processes associated with these prognostic morphological features were also identified based on proteomic data. The biological processes such as immune responses, cell cycle, and extracellular matrix development have been previously associated with cancer patient outcome. In summary, our work linked molecular data, morphology, and clinical out- come, which led to new insights and hypotheses into the relationships between cancer tissue development and molec- ular events, thus contributing to a more comprehensive un- derstanding of breast cancer. The entire process and work- flow can be applied to other cancers.



we performed correlation analysis between imaging features and mRNA or protein (MS- based global proteomic data) profiles by calculating Spear- man’s rank correlation coefficients (?). 

we compared the distribution of correlation coefficients
for image-mRNA and image-protein pairs for individual morpho- logical features

automated image analysis was carried out and ten types of cell-level features from tissue images were extracted following the three main steps: 1) nuclei segmentation, 2) cell-level feature measurement, and 3) aggregation of cell-level measurements into patient-level statistics.

In Step 1, the nuclei of all cells in the image are automatically segmented based on our previous workflow (31). In Step 2, ten types of cell-level features were ex- tracted, including seven types of morphological and spatial traits and three types of pixel traits in the RGB color space. The seven types of morphological and spatial features of cell nuclei were: major axis length (Major_Axis), minor axis length (Minor_Axis), the ratio of major to minor axis length (Ratio), nuclear area (Area), mean distance to neighboring cells (Mean_Distance), maximum distance to neighboring cells (Max_Distance), and minimum distance to neighboring cells 
(Min_Distance).

 ). The seven types of morphological and spatial fea- tures of cell nuclei can be summarized as nucleic area (Area), nucleic shape (Major_Axis, Minor_Axis, and Ratio), and cell density (Mean_ Distance, Max_Distance and Min_Distance). In Step 3, 5-bin histo- gram and five distribution statistics (i.e. mean, standard deviation or S.D., skewness, kurtosis, and entropy) were calculated for each of the ten types of morphological features to aggregate the measurements over the whole slide image. Thus for each type of feature, ten meas- urements (i.e. five histogram bins and five distribution statistics) were generated and 100 image features were generated in total for the ten types of morphological features. The centers of the five bins were determined by clustering each type of cell-level features from all patients instead of a single patient, which ensured that the histogram features are comparable and consistent across the entire patient cohort. The value of each feature based on the five bins of the histogram represented the relative percentage of corresponding im- age feature over the entire slide for a patient.



##### Yanovich et al. 2018

utilized mass spectrometry–based proteomic anal-ysis on more than 130 clinical breast samples to demon-strate intertumor heterogeneity across three breast cancer
subtypes and healthy tissue.

classification ambigu- ities still exist, and although multiple studies have been focusing onmolecularcharacterization ofbreast cancerprofiles, the clinical routine still relies primarily on immunostaining of ER, PR, and HER2 as a basis for classification (10).



Previously, the TCGA performed reversed-phase protein arrays (RPPA) in order to evaluate how genetic alterations affect the protein level. However, this method is limited to hundreds of pre-selected proteins and several phos- phorylation sites, and does not provide a system-wide view ofthe proteome. Analyzing



Analyzing the cancer proteome in an unbiased man- ner, two recent mass spectrometry (MS)-based proteomics studies analyzed dozens of tumors and evaluated the relevance of gene- expression levels to the proteome level. 



In the first study, we analyzed a cohort of 40 breast tumors derived from different subtypes, quantified the proteins using the super-SILAC technol- ogy (11), and reached the depth of more than 10,000 proteins



The second study, which was conducted by the Clinical Proteomic Tumor Analysis Consortium (CPTAC; ref. 13), ana- lyzed MS-based proteomics and phosphoproteomics data from a cohort containing 77 samples, used iTRAQ quantification and reached the depth of 12,553 proteins and 33,239 phosphoryla- tion sites. They further integrated the proteomic with the genomic profiles of the same cohort from the TCGA, and associated the mutational landscape with the corresponding proteomic and phosphoproteomic profiles. 

Both studies showed limited con- cordance between the overall protein and mRNA profiles. This was evident both in the 19-protein signature identified by Tya- nova and colleagues, in which a third was found to be regulated exclusively on the proteomic level, and also in the global protein– mRNA correlation analysis of the CPTAC, which resulted in a median Pearson correlation of 0.39

The larger scope of the CPTAC study enabled them to perform unsupervised classification of the tumors and identify novel proteomic and phosphoproteomic subgroups. Compared with the four expected mRNA-based subtypes, they identified only three proteomic groups, out of which two corresponded to the known classification (a basal cluster and a luminal cluster), whereas the third one was a new "stromal-enriched" cluster. This cluster was characterized by samples originating from all sub- types, and thus maybe associated with stroma-rich regions within the tumors

In this study, in order to explore intertumor heterogeneity and examine proteomics-based classification, we integrated three datasets from both published and unpublished resources, all generated using the same approach. The entire cohort consisted of 212 formalin-fixed paraffin-embedded (FFPE) samples origi- nating from different breast cancer subtypes, as well as healthy tissue. Out ofthese, 131 samples were selected for analysis, based onthe data qualityandsample annotation. Unsupervised analysis of the integrated cohort enabled us to identify a novel luminal breast cancer subtype. Comparing our findings to independent published data, we found that this classification is reproduced only on the proteomic level. This proteomic cohort enabled us to both support and challenge the known breast cancer classifica- tion, and thus to facilitate the analysis of similarities and differ- ences between distinct omics levels.



Batch effects originating from the dataset integration-related variability were corrected for using principal component analysis (PCA) and component subtraction. We found that the second component of a PCA separates between datasets (9% of the variation), and therefore subtracted this component before inte- grated downstream analysis.



For the **consensus clustering analysis**, ratios across samples and
across proteins were z-score normalized followed by their clus- tering using the consensus cluster algorithm (26). The algorithm performs multiple iterations ofhierarchical clustering, whereas in each iteration only a subset ofthe data is considered and the range of the number of output clusters is determined by the user. This results in a separate consensus matrix for each number ofclusters, describing the frequency oftwo samples clustering together out of the total number of iterations. We used R environment (version 3.2.3) and the ConsensusClusterPlus package (27)



Kyoto Encyclopedia ofGenes and Genomes (KEGG) pathway proteomaps



classified the 109 tumor samples using the consensus clustering algorithm (26).



Exami- nation of the distribution of the classical subtypes between the consensus clusters showed partial concordance with the known classification; most TN samples belonged to consensus cluster (CC) 4; HER2 samples spread across clusters; and ERþ samples (LumA or LumB) were divided into three clusters (Fig. 2B). These analyses, as well as t-SNE analysis (Fig. 2C), revealed that although CC1,CC2,andCC3are all associated withERþsamples, CC2 and CC3 are closer to the TN-enriched cluster, CC4; and CC1 is relatively segregated from them. Notably, running the consen- sus clustering algorithm with k¼ 3, the ERþCC2 was co-clustered (Supplementary Fig. S2C). To examine whether contaminating with the ER? CC4, thus indicating their higher global proximity we repeated the analysis and omitted all extracellular matrix (ECM) proteins from the dataset using Gene Ontology Cellular Component annotation. Re-analyzing with consensus clustering resulted in three clusters, reproducing CC1 and CC3, and creating a combined cluster of CC2 and CC4. This suggests that ECM proteins do not play a key role in the separation of the luminal samples (



Previous studies highlighted major differences in unsupervised
classification based on proteomic compared to the genomic approaches (1–3, 7). We therefore sought to examine whether the proteomic signature of differentially expressed proteins in CC2 and CC3 also separates the LumA samples based on their RNA profiles. In contrast to the separation based on the protein levels (Fig. 5C), hierarchical clustering oftheRNAdata ofthe same genes showed no clear differentiation between the intrinsic sub- types, even when sorted according to subtype (Fig. 5D). We applied the same analysis to an independent gene expression dataset, from the Metabric study (2012). Sorting the samples according to the 10 integrative clusters identified in the Metabric analysis and plotting the gene expression of the same signature yielded results similar to the TCGA RNA analysis, and did not manage to differentiate clearly between the clusters (Fig.





Unsupervised tumor classification using the consensus clustering algorithm separated the data into four clusters, consisting of TN-enriched group, and three ER+ enriched groups,  while HER2 samples were distributed across groups



The inability to identify a distinct group of HER2 samples in proteomics data might stem from the diverse phenotype of these tumors beyond HER2 expression. This result was previously shown by the CPTAC consortium (13) using similar techniques and a smaller cohort, and also by the Metabric study (2), which shows that HER2 is high in several integrative clusters based on copy number aberrations and mRNA expression. Additionally, a mass cytometry-based proteomic analysis of breast cancer also suggested that HER2 expressing cells do not represent a single phenotype, and are instead characterized byvariable expression of other proteins (44). The RPPA results did separate the HER2 subtype, presumably due to the pre-selection of a large number of HER2 related proteins in the arrays. In the current work, we consistently assessed HER2 and its associated proteins in an unbiased manner. They were indeed high in the HER2 clinical subtype, but showed variable expression across consensus clus- ters. While the importance ofHER2 is not questionable given the efficacy of HER2-directed therapies, our results suggest that this pathway affects only a subset of the proteome and does not represent a distinct subtype when examining the global proteome profiles. Thus, it is useful as a therapeutic marker rather than a representative of a defined tumor subtype. Additional

major discrepancies have been identified between all 'omics'- based classifications (including ours), which are based on aver- aging of large tumor regions, as compared to IHC approaches. Apart from potentially overlooking the effect of tumor microen- vironment components (e.g., ECM,immune infiltrates, and fibro- blasts), "averaging" the cancer cells masks the intratumor hetero- geneity that is well documented in breast cancer (45, 46). This intratumor heterogeneity is likely to underlie some of the
observed inconsistencies between classification using the clinical IHC approaches and the "omics" approaches





Comprehensive investigation of the three ERþ clusters identified one LumB cluster and two LumA clusters. One of the two 
LumA clusters was characterized with increased activity of TN- associated signaling proteins, thus suggesting a novel, sub- subtype of LumA tumors.

Our classification was validated using an independent dataset obtained from TCGA, which also assisted in the characterization of this novel subtype as potentially more aggressive and therapy resistant. The mechanism underlying this resistance relies on the coordinated interaction between ER itself and growth factor signaling protein

For example, KRT18 TN-enriched CC4; and the TN markers EGFR and CD44 were was significantly higher in all ERþ clusters compared to the significantly higher in CC4 when compared to CC3 and CC1, respectively. Examination of the LumB marker MKI67 suggested thatCC1 represents this subtype. This was further supported byan enrichment analysis that identified gene sets ofRB1-regulated and MYC-related targets, with which LumB subtype is known to be associated (1), as significantly enriched withinCC1samples when compared to either CC2 or CC3 (Supplementary Fig. S3A). This analysis also showed that the majority of significantly enriched processes are higher in CC1 in both comparisons. Very few processes are higher in CC2 and CC3, and the common ones include mostly interferon related pathways



Having associated CC4 with TN subtype and CC1 with LumB subtype, CC2 and CC3 remained as two possible LumA clusters

### Cao et al. 2021

In this study, an integrated strategy with a combination of tran- scriptomics, proteomics, glycomics and glycoproteomics was applied to explore the dysregulation of glycogenes, glycan structures and glycoproteins in chemoresistance of breast cancer cells. In paclitaxel (PTX) resistant MCF7 cells, 19 differentially expressed N-glycan-related proteins were identified, of which MGAT4A was the most significantly down-regulated, consistent with decrease in MGAT4A expression at mRNA level in PTX treated BC cells. Glycomic analysis consistently revealed suppressed levels of multi-antennary branching structures using MALDI-TOF/TOF-MS and lectin microarray. Several target glycoproteins bearing suppressed levels of multi- antennary branching structures were identified, and ERK signaling pathway was strongly suppressed in PTX resistant MCF7 cells. Our findings demonstrated the aberrant levels of multi-antennary branching structures and their target glycoproteins on PTX resistance.



N.B we are aware that not splitting the datasets is bias and will produce over-fitted models; however, here it is important to note that we do not use the results of **LASSO regression as a stand-alone method for selecting proteins**, but as a way of filtering results from differential abundance analysis

We performed LASSO regression with tenfold cross- validation, as implemented in the R-package, GLMNET [50]. We used contrasts from the differential expression analysis which yielded significant hits, including (a) BC subtype, (b) TIL status, and (c) hormone receptor sta- tus of ER and PgR as a response

The full set of pro- teins retained after filtering were used as input for the LASSO regression. We

We ran each LASSO model 10 times with 10 different random seeds and extracted the overlap of selected proteins across runs to a consensus set

We performed RF with the contrasts from the DAA, which yielded significant hits, including (a) breast can- cer subtype; (b) TIL status, and (c) hormone receptor status of ER and PgR as a response. 





### Liu et al. 2019

**independent component analysis (ICA)** to human breast cancer proteogenomics data to retrieve mechanistic information. We show that as an unsupervised feature extraction method, ICA was able to construct signatures with known biological relevance on both transcriptome and proteome levels. Moreover, pro- teome and transcriptome signatures can be associated by their respective correlation with patient clinical features, providing an integrated description of phenotype-related biological processes. 

**independent component analysis to proteomic and transcriptomic data of 77 breast cancer samples to extract pathway-level molecular signa- tures**



Independent component analysis (ICA) is an unsuper- vised learning method widely used in signal processing and has been applied to cancer genomics with notable success

This approach decomposes the molecular profiles into linear combinations of non-Gaussian independent sources or components, each of which is comprised of weighted contri- butions from individual genes. Therefore, ICA reduces the dimensionality of original data by representing the molecular profile of each sample as weighted sum of several “meta- genes” or “meta-proteins,” and the weight of specific meta- gene/protein (mixing scores) in one sample reflects the “ac- tivity” of that component in the sample.

Different from the more conventional dimension reduction method of principal component analysis (PCA), which seeks to find uncorrelated factors that explain the variance among the data, and works the best when the underlying components are normally dis- tributed, ICA are able to discover more informative represen- tations of high-dimensional biological signals, which are usually super-Gaussian and contain more close-to-zero values than a normally-distributed sample



As **clinical features are also available for the CPTAC samples, molecular signa- tures can be constructed from clusters of meta-genes/pro- teins that show activity patterns correlated with these clinical features**.



taking advantage of a specific clinical fea- ture as an “anchor,” this method may help extract patterns at different biological levels and across different cohorts, which may originate from the same cellular functionality



signatures extracted from different data sets were filtered based on their intrinsic stability and association with known clinical features (see Experimental Procedures below) and grouped into modules that showed similar correlation pat- terns to clinical features. 

Subsequent gene set enrichment analysis revealed the biological relevance of these modules to



Independent Component Analysis—Data were presented in a p ⫻ n
matrix X with genes in rows and samples in columns. The goal of ICA is to decompose the p ⫻ n data matrix into the product of a source matrix S (p ⫻ k) and a mixing matrix A (k ⫻ n). The ith column of the source matrix represents coefficients of each of the p genes for the ith inde- pendent component. The coefficient vector of each component could be considered as p random samples that revealed probability distribu- tion of a specific random variable. The mutual information between any pair of such variables is minimized, and the components are therefore statistically independent. Genes with coefficients of positive or negative values in the same component indicated that their levels may be en- hanced or suppressed in the same biological process. The ith row of the mixing matrix represents contributions of the ith component in the source matrix to profiles of each of the n samples. Rows of the mixing matrix therefore record the activity of each components (or meta-genes/ meta-proteins) across n samples.



All computations were carried out on the R platform. Package “fastICA” which implements the iterative FastICA algorithm (15) was used to extract non-Gaussian independent components with logcosh contrast function. Components were subsequently assigned to clusters using the “cluster” package. Clusters were visualized with 2d t-SNE using the R package “tsne.” The number of clusters was determined as equal to number of components extracted at each run of ICA. When the number of samples is small comparing to the number of features, which is usually the case for biological data, it is convenient to retrieve as many as independent signal sources as possible, and the number of compo- nents extracted is equal to sample size

The package “pcaMethods” was used to calculate principal com-
ponents for comparison with independent components. 

Each cluster was also associated with clinical features as following: First, 22 clinical features were recoded into ordinal variables (supplemental Table S1). Second, ordinary linear regression models were built with corre- sponding mixing scores for members in a component cluster to predict each of the ordinal responses. Counts of significant associa- tions between components and clinical features (p value for the slope coefficient ⬍ 5.9 ⫻ 10⫺7, Bonferroni correction at ␣⫽ 0.01 level) were tabulated. Hierarchical clustering with complete linkage was applied to the clinical associations of independent components clusters ex- tracted from both transcriptome and proteome data.

The stability of these signatures could be inspected by visu- alizing all meta-genes and meta-proteins with t-distributed stochastic neighbor embedding (t-SNE), a widely used dimen- sionality reduction technique (19) 

We recoded 22 clinical features into ordinal factors and used linear regres- sion to assess their correlation with activity scores of meta- genes or meta-proteins in each signature cluster. To select the most significant associations, Bonferroni correction was applied to control for multiple comparisons at ␣ ⫽ 0.01 level, such that a nominal p value of 5.9 ⫻ 10⫺7
was set as the
significance threshold (corrected ␣ ⫽ 0.01/(77 ⫻ 22)). Within

The biological relevance of the clusters identified with ICA
could be further validated by inspecting the average gene coefficients (the centroids of the clusters).

The corre- lation between mixing scores and clinical features provided a valuable opportunity to find the link between independently extracted components from different data types



**To integrate**
**the proteomic and transcriptomic information with guidance from clinical relevance, we applied the hierarchical clustering algorithm to the vectors of clinical association counts for the most clinically relevant meta-protein and meta-gene clusters**



Proteome and transcriptome signatures could therefore be grouped based on their similarity in functional indications, as the metric of direct correlation between gene coefficients of different signatures was limited by noise introduced by high dimensionality (supplemental

Combined network vi- sualization of GO terms enriched in proteome and transcrip- tome signatures allowed us to examine complementary func- tional modules on different biological levels.

### Krug et al. 2020



Proteogenomic analyses of 122 primary breast cancers provide insights into clinically relevant biology, including cell cycle dysregulation, tumor immunogenicity, aberrant metabolism, and heterogeneity in therapeutic target expression.

‘‘proteogenomics’’ approachwas applied to 122 treatment-naive primary breast cancers accrued to preserve post-translational modifications, including protein phosphoryla- tion and acetylation.

Proteogenomics challenged standard breast cancer diagnoses, provided detailed anal- ysis of the ERBB2 amplicon, defined tumor subsets that could benefit from immune checkpoint therapy, and allowed more accurate assessment of Rb status for prediction of CDK4/6 inhibitor responsiveness

. Phospho- proteomicsprofilesuncoverednovelassociationsbetweentumorsuppressorlossandtargetable kinases.Ace- tylproteomeanalysis highlightedacetylationonkeynuclear proteins involved in theDNAdamage responseand revealed cross-talk between cytoplasmic and mitochondrial acetylation and metabolism. Our

Proteogenomics is an approach to tumor profiling that com-
bines next-generation DNA and RNA sequencing with mass spectrometry-based proteomics to provide deep, unbiased quantification of proteins and post-translational modifications such as phosphorylation

The Clinical Pro- teomic Tumor Analysis Consortium (CPTAC) seeks to perform deep-scale proteogenomics profiling across multiple cancer types. Our initial proteogenomics analysis of BRCA using resid- ual samples from The Cancer Genome Atlas (TCGA) provided proof of principle that proteogenomics represented an advance in BRCA profiling (Mertins

proteogenomics characterization of the largest cohort to date of BRCA samples that were acquired to minimize ischemic time, maximizing fidelity and reducing pre- analytical variability. We offer the first comprehensive report of the BRCA acetylome;

The PAM50 model was applied to RNA sequencing (RNA-seq)
data to determine representation of intrinsic subtypes (Parker

Somatic copy number alteration (SCNA) data were analyzed to detect focal and arm-level events (Mermel et al., 2011; Figures S2C and S2D) with confirmation of anticipated effects on mRNA and protein abundance 

To explore intrinsic cohort structure using the full complement of proteogenomics data, single-omic and multi-omics clustering were performed for SCNA, mRNA, protein, and individual phos- phosite and acetylation site abundance using non-negative ma- trix factorization (NMF) (Lee

Although NMF yielded between two and six clusters in single-omic ana- lyses (Figure S3A), integrative multi-omics analysis converged on four NMF clusters, with	

Clusters designated luminal A-inclusive (NMF LumA-I) and basal-inclusive (NMF Basal-I) were almost entirely composed of tumors with the cor- responding PAM50 assignments



Two clusters showed sample compositions that were discor-
dant with PAM50 subtypes. The luminal B-inclusive cluster (NMF LumB-I) comprised all but one LumB case but also included a subset of PAM50 LumA samples. A

The two luminal clusters also showed remarkable dichotomies in pathway space, supporting the concept that, although heterogeneous, these are biologically separate tumor types. For example, cancer hallmark gene set enrichment scores for LumA-I versus LumB-I were significantly anti-correlated even though estrogen response-related terms were positively enriched in both (Figures

Notably, a mixed PAM50 LumA/B cluster was also observed when clus- tering the global RNA data in isolation, indicating that PAM50 classification, a method simplified for clinical purposes, does not capture all biological distinctions between LumA and LumB

random
forest classifiers were trained on protein or mRNA data to distin- guish PAM50 LumA samples assigned to the NMF LumB-I clus- ter from PAM50 LumA samples assigned to the NMF LumA-I cluster. When these classifiers were applied to METABRIC data (Curtis et al., 2012), samples from patients with NMF fea- tures that drove PAM50 LumA samples into the NMF LumB-I cluster had outcomes that were intermediate between the re- maining PAM50 LumA samples and the PAM50 LumB samples (Figure 1C; Figure S3K). This finding supports the NMF assign- ment of some PAM50 LumA samples to the higher-risk LumB-I cluster.

HER2-inclusive cluster (NMF HER2-I) was remarkably het- erogeneous. Although predominantly composed of HER2-en- riched PAM50 subtype samples and samples with centrally confirmed, clinically positive ERBB2 status, NMF HER2-I also included tumors from all four other PAM50 subtypes, suggesting the presence of unifying biological features in NMF informatic space that are absent in the PAM50-based classification

To identify putative therapeutic targets specific for each NMF subtype, phosphoproteomic data were used as kinase activation surrogates

Phosphorylated kinases enriched in each NMF subtype were identified using outlier enrichment analysis (Black- Sheep Python package) (Blumenberg

The BlackSheep approach also associated phosphorylated ki-
nase outliers with recurrent somatic mutations (Figure

Tumor metabolic characteristics were profiled at the level of the proteome, and unsupervised clustering of differentially expressed (DE) metabolism-related proteins

Protein acetylation (Ac) has been implicated in cellular meta-
bolism in addition to roles in epigenetic regulation



RNA-based immune cell deconvolution signatures and protein- level signatures for immune modulators (Thorsson et al., 2019) revealed a range of immune-related features across all four intrinsic subtypes (Table S6), including the immune checkpoint proteins PD1 and PD-L1 at the RNA and phosphosite levels

. RNA level profiles inferred for individual acquired and innate immune cell types generally tracked with CIBERSORT absolute scores in each sub- type (

RNA-based estimates of overall I-TME provided by CI- BERSORT absolute scores



Finally, stromal, fibroblast, mast cell, endothelial cell, and neutrophil signatures were elevated in PAM50 LumA tumors with higher CIBERSORT scores but lower overall in LumB and Basal tumors



To identify potential drivers of immunogenicity across common BRCA subtypes, PD-L1 mRNA levels were correlated separately with proteomics data from PAM50 luminal and basal cases 



Proliferation rate is a critical prognostic feature in BRCA, and the cell cycle is a target for endocrine therapy (Ellis et al., 2017) and CDK4/6 inhibition in ER+, ERBB2? advanced BRCA (Pernas et al., 2018).



To compare PG fea- tures with cell cycle control in hormone receptor (HR)+/ERBB2 PG? and triple-negative BRCA (TNBC) tumors, the multi-gene proliferation score (MGPS; Figure 6A; Table S6) was generated for each sample (Ellis

e high-quality, multi-omics resource we created allows inves- tigators to explore correlations between the genomic landscape and the downstream effects in the BRCA proteome, phospho- proteome, and acetylproteome, extending and refining analytical opportunities provided by prior studies (Bouchal et al., 2019; Jo- hansson et al., 2019; Mertins et al., 2016; Tyanova et al., 2016).
Figure



he entire workflow described under ‘Multi-omics clustering’ has been implemented as a module for Broad’s cloud platform Terra
(https://app.terra.bio/). The docker containers encapsulating the source code and required R-packages for NMF clustering and ssGSEA have been submitted to Dockerhub



**non-negative matrix factorization algorithm (NMF) was used to decipher de novo mutation signa- tures in cancer somatic mutations stratified by 96 base substitutions in tri-nucleotide sequence contexts**



SignatureA- nalyzer exploited the **Bayesian variant of the NMF algorithm** and enabled an inference for the optimal number of signatures from the data itself at a balance between data fidelity (likelihood) and model complexity (regularization) (Kasar

the inferred signatures were compared against known signatures derived from COSMIC (Tate et al., 2019) and cosine similarity was calculated to identify the best match. Mutational



APOBEC enrich- ment was also assessed using TrinucleotideMatrix and plotApobecDiff functions of the maftool package (Mayakonda



Non-negative matrix factorization (NMF) implemented in the NMF R-package (Gaujoux and Seoighe, 2010) was used to perform un- supervised clustering of tumor samples and to identify proteogenomic features (proteins, phosphosites, acetylsites, RNA transcripts and somatic copy number alterations) that showed characteristic abundance patterns for each cluster. Briefly, given a factorization rank k (where k is the number of clusters), NMF decomposes a px n data matrix Vinto two matrices Wand H such that multiplication of Wand H approximates V. Matrix H is a kx n matrix whose entries represent weights for each sample (1 to N) to contribute to each cluster (1 to k), whereas matrix Wis a px k matrix representing weights for each feature (1 to p) to contribute to each cluster (1 to k). Matrix H was used to assign samples to clusters by choosing the k with maximum score in each column of H. For each sample, we calculated a cluster membership score as the maximal fractional score of the corresponding column in matrix H. We defined a ’’clus- ter core’’ as the set of samples with cluster membership score > 0.5. Matrix W containing the weights of each feature in a certain cluster was used to derive a list of representative features separating the clusters using the method proposed in Kim and Park (2007). Cluster-specific features were further subjected to a 2-sample moderated t test (Ritchie et al., 2015) comparing the feature abundance between the respective cluster and all other clusters. Derived p values were adjusted for multiple hypothesis testing using the methods proposed in Benjamini and Hochberg (1995)



To enable integrative multi-omics clustering, we required all data types (and converted if necessary) to represent ratios to either a common reference measured in each TMT plex (proteome, phosphoproteome, acetylproteome) or an in-silico common reference calculated as the median abundance across all samples (mRNA, see ‘‘RNA quantification’’). All data tables were then concatenated and only features quantified in all tumors were used for subsequent analysis. Features with the lowest standard deviation (bottom 5th percentile) across all samples were deemed uninformative and were removed from the dataset. Each row in the data matrix was further scaled and standardized such that all features from different data types were represented as z-scores.



Since NMF requires a non-negative input matrix, the data matrix of z-scores was further converted into a non-negative matrix as follows:
1) Create one data matrix with all negative numbers zeroed. 
2) Create another data matrix with all positive numbers zeroed and the signs of all negative numbers removed.
3) Concatenate both matrices resulting in a data matrix twice as large as the original, but with positive values only and zeros and hence appropriate for NMF.
    The resulting matrix was then subjected to NMF analysis leveraging the NMF R-package (Gaujoux and Seoighe, 2010) and using
    the factorization method described in Brunet et al. (2004).

To determine the optimal factorization rank k (number of clusters) for the multi-omic data matrix, a range of clusters between k = 2 and 8 was tested. For each k we factorized matrix Vusing 50 iterations with random initializations of Wand H. To determine the optimal factorization rank we calculated two metrics for each k: 1) cophenetic correlation coefficient measuring how well the intrinsic structure of the data was recapitulated after clustering and 2) the dispersion coefficient of the consensus matrix as defined in Kim and Park (2007) measuring the reproducibility of the clustering across 50 iter- ations. The optimal k was defined as the maximum of the product of both metrics for cluster numbers between k = 3 and 8 



Having determined the optimal factorization rank k, and in order to achieve robust factorization of the multi-omics data matrix V, the
NMF analysis was repeated using 1000 iterations with random initializations of Wand H and partitioning of samples into clusters as described above. Due to the non-negative transformation applied to the z-scored data matrix as described above, matrixWof feature weights contained two separate weights for positive and negative z-scores of each feature, respectively.



In order to reverse the non- negative transformation and to derive a single signed weight for each feature, each row in matrixWwas first normalized by dividing by the sum of feature weights in each row. Weights per feature and cluster were then aggregated by keeping the maximal normalized weight and multiplying with the sign of the z-score from the initial data matrix. Thus, the resulting transformed version of matrixWsigned contained signed cluster weights for each feature present in the input matrix



In order to functionally characterize the clustering results, normalized enrichment scores (NES) of cancer-relevant gene sets were calculated by projecting the matrix of signed multi-omic feature weights (Wsigned) onto Hallmark pathway gene sets (Liberzon et al., 2015) using ssGSEA (Barbie et al., 2009). To derive a single weight for each gene measured across multiple omics data types (protein, RNA, phosphorylation site, acetylation site) we retained the weight with maximal absolute amplitude



### Bouchal et al. 2019

Proteotyping of 96 breast tumors by SWATH mass spectrometry

we present a quanti- tative proteotyping approach based on sequential windowed acquisition of all theoretical fragment ion spectra (SWATH) mass spectrometry and establish key proteins for breast tumor classification. The study is based on 96 tissue samples representing five conventional breast cancer subtypes. SWATH proteotype patterns largely recapitulate these sub- types; however, they also reveal varying heterogene- ity within the conventional subtypes, with triple nega- tive tumors being the most heterogeneous. Proteins that contribute most strongly to the proteotype- based classification include INPP4B, CDK1, and ERBB2 and are associated with estrogen receptor (ER) status, tumor grade status, and HER2 status

SWATH-MS Data Matrix Consisting of 2,842 Consistently Quantified Proteins across 96 Patient Samples

We performed **unsupervised hierarchical clustering on the proteotypes of the pooled samples**. Figure 1A shows that pools of lymph-node-positive and negative samples of each subtype clustered closely together, indicating high reproduc- ibility of our measurements. Moreover, clustering revealed proteotype similarity between less aggressive luminal A and luminal B subtypes, whereas the more aggressive HER2 and triple-negative subtypes formed a separate cluster. The luminal B HER2+ group was more similar to the cluster with high aggressiveness, in agreement with its worse therapy response

we systematically correlated the quantitative proteo-
types of the 96 individually measured breast cancer samples and ordered the resulting correlation coefficients according to the classical tumor subtypes

we found a high correlation of some samples of the HER2- enriched subtype with some (mostly lymph-node-positive) luminal B HER2+ samples, indicating that a higher degree of similarity in HER2+ tumors of luminal B and HER2-enriched subtypes were apparent from the proteotype

we found that clus- tering by proteotypes closely recapitulates conventional tumor subtyping, but we also found that some of these subtypes are more heterogeneous (triple-negative tumors) than others

we next constructed **a decision tree** to classify the 96 tumors into the five conventional subtypes based on their proteotypes

we compared our comprehensive SWATH-MS da- taset against the five microarray datasets of 883 patients mentioned above

6% of protein-level observations and 7%–15% of transcript-level observations (depending on the set of patients) exhibited statistically significant changes

13%–28% of differentially abundant proteins also showed a statistically significant change with the same di- rection on the transcript level. From the reverse perspective, 9%–18% of significantly regulated transcripts showed a signifi- cant change with the same trend also on protein level

A decision tree constructed from the five independent transcriptomics datasets using expression data for 1,036 genes resulted in a tree with three nodes and similar structure



##### Yasar et al. 2020

The approach used in variable selection is the **Penalized Logistic Regression model**. Penalized Logistic Regression adds a penalty parameter (also known as regularization) to logistic regression, which has a lot of variables so that the coefficients of variables that contribute less to the model are zero. The most frequently used punished regression model is **Lasso regression**. The coefficients of some variables that contribute less to the Lasso regression are forced to be exactly zero. Only the most important variables are kept in the final model. The lasso logistic regression model established in the study has been applied 100 times, and variables that included more than 15 in the model have been included in the study

The deep learning model constructed in this study is based on a multi-layer feed-forward artificial neural network trained with stochastic gradient descent using the back- propagation technique.

**The aim of this study is to classify the molecular subclasses of breast cancer using a deep learning algorithm based on proteomic data.**

In the deep learning model based on proteomics data, the protein code and relative importances of the top 10 proteins, which are also important in classifying subtypes of breast cancer



##### Osz et al. 2021

integrate available proteome-level breast cancer datasets to identify and validate new prognostic biomarker candidates.

We established a database integrating protein expression data and survival information from four independent cohorts for 1229 breast cancer patients.

(n proteins = 7'432)



### Lawrence et al. 2015

integrate these results with data generated in-house and from publicly accessible genomics and drug sensitivity resources



quantitative proteomics anal- ysis of 20 human-derived breast cell lines and 4 primary breast tumors to a depth of more than 12,000 distinct proteins

We integrated proteomics data with exome sequence resources to identify genomic aberrations that affect protein expression. We per- formed a high-throughput drug screen to identify protein markers of drug sensitivity and understand the mechanisms of drug resistance

We assembled a panel of 20 human breast cell lines and four clinical tumors to analyze the proteomic landscape of TNBC (Fig- ure 1A). These included 16 triple-negative cell lines covering mesenchymal-, luminal-, and basal-like subtypes, as well as three receptor-positive and one non-tumorigenic cell line to serve as a basis for comparison (Lehmann

To facilitate use and dissemination of the data, we have developed a web resource (https://zucchini.gs. washington.edu/BreastCancerProteome/) in which protein abundances can be queried and correlated to genomic and drug sensitivity data, as presented below. 

We used hierarchical clus- tering to identify patterns based on correlation of protein expres- sion profiles. This approach classified the panel of cell lines into two overarching groups containing four clusters (Figure 3A). To illustrate the relationship between driver gene alterations and proteome profiles, we show the most frequent census mutations and copy-number aberrations for each cell line

As has been observed previously (Cancer Genome Atlas Network, 2012), PIK3CA mutations were associated with luminal breast cancer subtypes (80% of the cell lines in cluster 1), whereas TP53 mutations were characteristic of TNBC (100% of the cell lines in clusters 3 and 4). Mutations in the tumor suppressor NF1 were exclusive to the mesen- chymal-like subtype (cluster 4) and BCR mutations were exclu- sive to luminal cells (cluster 1).

To better illustrate this, we used principal component analysis (PCA) to project the distances between each proteome onto a two-dimensional coordinate system. Some of the sample proteomes formed tight clusters, while others were more distantly related to those in the same group

The composition of each cluster showed striking similarity to subtypes defined by mRNA expression arrays and morpholog- ical studies

To better understand the biology of each subtype, we compared the distribution of protein abundance within gene ontology categories. Interest- ingly, luminal-like cells expressed higher levels of pathways associated with proliferation, such as cell cycle, growth factor signaling, metabolism, and DNA damage repair mechanisms (Figure 3E; Figure S3B). TNBC cell types, particularly the tumors and more invasive cells, expressed higher levels of pathways associated with metastasis, such as ECM-receptor interaction, cell adhesion, and angiogenesis (Figure 3E; Figure S3B). The ex- pressions of proliferation and metastasis pathways were mutu- ally exclusive, an observation also made in an analysis of mRNA expression profiles from claudin-low tumors (Prat et al., 2010). Thus, therapies targeting immune and metastatic sig- naling are an exciting avenue for TNBC treatment

despite overall concor- dance of whole proteome profiles with various cellular pheno- types, in most cases the expression of particular cancer proteins did not uniformly belong to one subtype or another



The identification and quantification of protein isoforms resulting from alternative splicing is a significant challenge in proteomics, arising from the reduced number of isoform-specific peptides that are amenable to analysis by mass spectrometry. For this da- taset, we first relied on isoform-specific peptides to unambigu- ously identify proteins mapping to the same gene in the Uniprot sequence database. This led to the identification of 1,860 protein isoforms that corresponded to 844 genes, 52 of which were members of the COSMIC census. Next, we examined the relative quantification of protein isoforms. Protein isoforms share long segments of identical sequence but are missing certain protein domains, resulting in altered signal intensity from those parts of the protein.



We integrated publicly avail- able exome sequence and gene CN data from COSMIC (Forbes et al., 2011) with proteome profiles from 18 cell lines. Protein abundance trended positively with gene CN.



We reasoned that mutations in certain driver genes, such as
those in the same signaling pathway, would likely converge to regulate common effectors. To determine the global effects of driver gene mutations on protein expression, we systematically evaluated gene-protein associations for frequently mutated census genes (n R 3 cell lines) by comparing the abundance of each protein in cell lines with versus without a mutation, and plotted this information as a network.

it demonstrates that dysregulation of cell-cycle protein abundance may be a common effect of diverse genetic mutations.

We reasoned that mutations in certain driver genes, such as
those in the same signaling pathway, would likely converge to regulate common effectors. To determine the global effects of driver gene mutations on protein expression, we systematically evaluated gene-protein associations for frequently mutated census genes (n R 3 cell lines) by comparing the abundance of each protein in cell lines with versus without a mutation, and plotted this information as a network.

Proteomics of Drug Sensitivity To generate a resource for drug sensitivity prediction, we screened the 16 TNBC cell lines from our panel against a library of 160 compounds at eight different concentrations spanning four orders of magnitude. We used this data to determine the IC50, defined as the dose required to reach a 50% reduction in cell viability, for each drug in each cell line



Next, we combined our pharmacological dataset with publicly accessible data from the Genomics of Drug Sensitivity in Cancer (CRx) resource (Yang et al., 2013) and performed regression analysis against mass spectrometry-derived protein abun- dances to discover proteomic markers of drug sensitivity or resistance. We used hierarchical clustering to analyze global patterns among drug sensitivity-protein expression relation- ships, revealing many distinct clusters (Figure 7B). Drugs target- ing proteins in the same pathway (e.g., BRAF andMEK inhibitors) showed similar correlation profiles. Interestingly, proteins that were part of the same pathways or complexes also clustered together, which did not occur using protein expression data alone



**the integration of proteomics and drug sensitivity data using regression analysis** provides a rich resource to iden- tify unexpected modes of action and to discover new features of target pathways.

to demonstrate the potential clinical utility of these re-
sults, we asked how many proteins from the drug association analysis could be identified in primary tumors. We found that 73% (6,798/9,292) were quantifiable in the four clinical speci- mens we analyzed (Figure



##### Kalocsay et al. 2020

We performed quantitative proteomics on 61 human-derived breast cancer cell lines to a depth of ~13,000 proteins. The

. All datasets are freely available as public resources on the LINCS portal.

e https://lincs.hms.harvard.edu/db/datasets/20343



### **Wu et al. 2021**

**SCSubtype: intrinsic subtyping for scRNA-seq data**. Since unsu- pervised clustering could not be used to find recurring neoplastic cell gene expression features between tumors, we asked whether we could classify cells using the established PAM50 method. Due to the inherent sparsity of single-cell data, we developed a scRNA-seq-compatible method for intrinsic molecular subtyp- ing. We constructed ‘pseudobulk’ profiles from scRNA-seq for each tumor and applied the PAM50 centroid predictor. To identify a robust training set, we used hierarchical clustering of the pseu- dobulk samples with The Cancer Genome Atlas (TCGA) dataset of 1,100 breast tumors using an approximate 2,000-gene intrinsic breast cancer gene list3
(Extended Data Fig. 2a,b). Training samples
were selected from those with concordance between pseudobulk PAM50 subtype calls and TCGA clusters



To further validate SCSubtype, we calculated the degree of epi- thelial cell differentiation (DScore)33
and proliferation34 , both of
which are independently associated with the molecular subtype of each cell. Basal_SC cells tended to have low DScores and high proliferation scores whereas LumA_SC cells showed high DScores and low proliferation scores (Fig. 2d and Extended Data Fig. 2e), as observed across PAM50 subtypes in TCGA (Extended Data Fig. 2f)

we inves- tigated the biological pathways driving intratumor transcriptional heterogeneity (ITTH) in an unsupervised manner, using integrative clustering of tumors with at least 50 neoplastic cells, to generate 574 gene signatures of ITTH. These gene signatures identified seven robust groups, ‘gene modules’ (GMs), based on their Jaccard simi- larity (Extended Data Fig. 3a). Each GM was defined with 200 genes that had the highest frequency of occurrence across the ITTH gene signatures and individual tumors (Supplementary Table 5), mini- mizing the contribution of a single tumor to any particular module

Gene-set enrichment identified shared and distinct functional
features of these GMs (Fig. 2e). GM4 was uniquely enriched for hall- marks of cell cycle and proliferation (for example, E2F_TARGETS)  driven by genes including MKI67, PCNA and CDK1

Spatially mapping breast cancer heterogeneity. To gain insights into the spatial organization of cell types, we performed spatially resolved transcriptomics on six samples



We
earlier showed that GMs were enriched for distinct
microenvironment-associated pathways and factors; thus, we hypothesized that GMs would be spatially organized in breast tumors. We

we define the cellular architecture of breast tumors at three levels. First, a detailed cellular taxonomy that includes new cell types and states and new methods for characterizing cellular heterogeneity (Fig. 8g). Second, a spatial map of cellular locations and interactions within tumors that reveals coor- dination of tumor and host cell phenotypes within tissue and reveals spatial relationships between cells. Third, using deconvolution, we observed groups of tumors with similar cell type proportions and prognostic associations, named ecotypes, often driven by spe- cific clusters of co-segregating cells



We **deconvoluted all primary breast tumor** datasets in the
METABRIC cohort40
. Supporting the validity of the predictions,
and SCSubtype, we observed significant enrichment (Wilcoxon test, P < 2.2 × 10−16
) of the four SCSubtypes in tumors with matching
bulk-PAM50 classifications. Significant enrichment (Wilcoxon test, P < 2.2 × 10−16
) of cycling cells in basal, LumB and HER2E tumors was
also shown (Extended Data Fig. 10c). Consensus clustering revealed nine tumor clusters with similar estimated cellular composition (‘ecotypes’) (Fig. 8a–c). These ecotypes displayed correlation with tumor subtype, SCSubtype cell distributions and a diversity of major cell types (Fig. 8a). Ecotype 3 (E3) was enriched for tumors contain- ing Basal_SC, Cycling and Luminal_Progenitor cells (the presump- tive cell of origin for basal breast cancers28
) and a basal bulk PAM50
subtype (Fig. 8a,b). In contrast, E1, E5, E6, E8 and E9 consisted pre- dominantly of luminal cells. Ecotypes also possessed unique patterns of stromal and immune cell enrichment. E4 was highly enriched for immune cells associated with antitumor immunity (Fig. 8a), includ- ing exhausted CD8 T cells (LAG3/c8), along with TH
1 (IL7R/c1) and
central memory CD4 T cells (CCR7/c0). E2 primarily consisted of LumA and normal-like tumors (Fig. 8b) and was defined by a clus- ter of mesenchymal cell types, including endothelial CXCL12+
and
ACKR1+ cells, s1 MSC iCAFs and depletion of cycling cells (Fig. 8a)



we investigated the association between ecotypes and the
integrative genomic clusters identified in the METABRIC cohort40 (Extended Data Fig. 10l). E3 had a high proportion of cancers from integrative genomic cluster 10, which also predominantly con- sisted of basal-like tumors with similarly poor 5-year survival. E7 had a high proportion of ERBB2-amplified and HER2E integrative genomic cluster 5 tumors. These were the worst prognosis groups in both the METABRIC and ecotype analyses. However, most ecotypes did not clearly associate with a specific integrative genomic cluster or PAM50 subtype, which is reflected by the role of stromal and immune cells in defining ecotypes. This lack of unique association suggests that ecotypes are not a simple surrogate for molecular or genomic subtypes



**Pseudotemporal ordering to infer cell trajectories.** Cell differentiation was inferred for mesenchymal cells (CAFs, PVLs and endothelial cells) using the Monocle 2 (ref. 46
) v.2.10.1 with default parameters as recommended by the
developers. Integrated gene expression matrices from each cell type were first exported from Seurat v.3 into Monocle to construct a CellDataSet. All variable genes defined by the differentialGeneTest function (cutoff of q < 0.001) were used for cell ordering with the setOrderingFilter function. Dimensionality reduction was performed with no normalization methods and the DDRTree reduction method in the reduceDimension step



Differential gene expression, module and pathway enrichment. Differential gene expression was performed using MAST67
v.1.8.2. All DEGs from each
cluster (log fold change >0.5, P threshold of 0.05 and adjusted P threshold of 0.05; Supplementary Tables 9 and 10) were used as input into the ClusterProfiler package69
v.3.14.0 for gene ontology functional enrichment. Results were clustered,
scaled and visualized using the pheatmap package v.1.0.12. Cytotoxic, TAM and dysfunctional T cell gene expression signatures were assigned using the AddModuleScore function in Seurat v.3.0.0 (ref. 36
). The list of genes used for
dysfunctional T cells were adopted from Li et al.38. The TAM gene list was adopted from Cassetta et al.10
. The cytotoxic gene list consists of 12 genes that translate to
effector cytotoxic proteins (GZMA, GZMB, GZMH, GZMK, GZMM, GNLY, PRF1 and FASLG) and well-described cytotoxic T cell activation markers (IFNG, TNF, IL2R and IL2).



**Tumor ecotype analysis using deconvolution.** CIBERSORTx57 v.1.0 and DWLS58 (accessed from https://github.com/dtsoucas/DWLS on 30/11/2020) were
used to deconvolute predicted cell fractions from a number of bulk transcript profiling datasets (see Supplementary Note for specific parameters). To prevent confounding of cycling cell types, we first assigned all neoplastic epithelial cells with a proliferation score greater than 0 as cycling and then combined these with cycling cell states from all other cell types to generate a single cycling cell state. Normalized METABRIC expression matrices, clinical information and PAM50 subtype classifications were obtained from METABRIC (https:// www.cbioportal.org/study/summary?id=brca_metabric). Tumor ecotypes in the METABRIC cohort were identified using SKmeans-based consensus clustering (as implemented in cola v.1.2.0) of the predicted cell fraction from either CIBERSORTx or DWLS in each bulk METABRIC patient tumor. When comparing ecotypes between methods (that is, consensus clustering results from using the cell abundances of all cell types or just the 32 significantly correlated cell types from CIBERSORTx deconvolution and the consensus clustering results from CIBERSORTx or DWLS cell abundances), the number of tumor ecotypes was fixed as 9 and the tumor overlaps between all ecotype pairs was calculated (Supplementary Tables 7 and 8). Common ecotypes were then identified by identifying the ecotype pairs with the largest average METABRIC tumor overlap. Differences in survival between ecotypes were assessed using Kaplan–Meier analysis and log-rank test statistics using the survival v.2.44-1.1 and survminer v.0.4.7 R packages



All processed scRNA-seq data are available for in-browser exploration and download through the Broad Institute Single Cell portal at https://singlecell. broadinstitute.org/single_cell/study/SCP1039. Processed scRNA-seq data from this study are also available through the Gene Expression Omnibus under accession number GSE176078



Code related to the analyses in this study can be found on GitHub at https://github. com/Swarbricklab-code/BrCa_cell_atlas





JacLy: a Jacobian-based method for the inference of metabolic interactions from the covariance of steady-state metabolome data



Onecommon approach that utilizes inherent variability in steady-state data is correlation
based inference methods, especially the Gaussian Graphical Model (GGM) approach (Çakır et al., 2009; Montastier et al., 2015; Wang et al., 2016). Correlation based approaches are capable ofdetecting strong interactions in the metabolic network to some extent. However, they infer interactions only in undirected manner, and they have limited power in the detection of weak interactions (Çakır et al., 2009). A directed network inference approach from steady-state metabolome data is also available in the literature (Steuer et al., 2003; Öksüz, Sadıkoğlu & Çakır, 2013; Çakır & Khatibipour, 2014). The approach is based on the prediction of interaction strengths from the covariance of the data. The network structure information stored in the inherent variability of the data is reflected on the covariance of the data, and later used in the prediction of interaction strengths in the form of a Jacobian matrix. The Jacobian matrix of a cellular interaction system contains a significantly high amount of valuable information both on the structure and dynamic characteristics of the system. Numbers in this matrix easily provide us with detailed information on the underlying interactions in the network, such as direction of interactions, nature of interactions (positive or negative effects), and strengths of interactions (Steuer et al., 2003; Öksüz, Sadıkoğlu & Çakır, 2013).



### Cao et al. 2021 [pancreas]

To investigate the association between methylation and proteomics expression, for each gene, we first calculated Z scores for its mRNA expression, protein, and phosphorylation levels and beta values for DNA methylation. We then calculated Pearson correlation scores with its associated significance between methylation and gene expression, protein, and phosphorylation levels for all pairs of genes, respectively.
e10

Non-negative matrix factorization (NMF)-based multi-omics clustering was performed similar to as previously described (Gillette

Non-negative matrix factoriza-tion (NMF)-based proteogenomics subtyping using multiomics
data from all 140 tumors



### Valle et al. 2020

Topic modeling was introduced to classify texts of natural language by inferring their topic structure from the frequency of words. This paper assumes that analogously the cancer subtype identity, which is crucial for the correct diagnosis and treatment plan, can be extracted from gene expression patterns with similar techniques. Focusing on breast and lung cancer, we show that state-of-the-art topic modeling techniques can successfully classify known subtypes and identify cohorts of patients with different survival probabilities. The topic structure hidden in expression data can be looked at as a biologically relevant low-dimensional data representation that can be used to build efficient classifiers of expression patterns.

The problem of finding a topic structure in a dataset was recently recognized to be analogous to the community detection problem in network theory

algorithms which try to identify the “topics” associated to a given document from the word usage have to face the same type of challenges we are facing here. In this analogy, the cancer samples play the role of the documents, the words are the genes, the number of times a particular word is used in a given document is the analogous of the expression level of a particular gene in a given sample. The topics are the gene sets (the “signatures”) we use to cluster samples into subtypes. The goal of topic modeling is to identify the “topic” of a given document from the word usage exactly as our goal is to identify the cancer subtype from the gene expression pattern.

The major advantage of topic modeling methods with respect to standard clustering approaches is that they allow a “fuzzy” type of clustering [7]. The output of a topic modeling algorithm is a probability distribution of membership, i.e., the probability of a given document to be composed by a given set of topics and, at the same time, the probability of a word to characterize a given topic. In our context, this means that we have as output a set of values that quantify the probability of a given sample to belong to a particular cancer subtype and the relevance of a given gene in driving this identification

The most popular tool to perform this kind of analysis is the so-called Latent Dirichlet Allocation
(LDA) algorithm, which basically assumes a Dirichlet prior distribution for the topic distribution. There is no particular motivation in the natural language context as well as in our biological context for the Dirichlet prior. Its motivation comes from the fact that it allows a simple solution of the allocation problem and thus the algorithm can be efficiently applied to databases with a large number of documents and words. However, the lack of biological motivation for the prior and the large number of free parameters represent a possible limitation of LDA in our context.

Another common approach, often used in addressing gene expression data, is the Nonnegative
Matrix Factorization [8]. The major drawback of this approach is that it facilitates the detection of sharp boundaries among subtypes and this could be a limitation in very heterogeneous settings

In recent years, some important advances in the field have laid the foundations for overcoming
some of the limits of LDA. First, it has been realized [9,10] that there is a strong connection between topic modeling analysis of complex databases and the community detection problem in bipartite graphs, which is a well know and much studied problem in network theory [11]. Second, a very effective class of community detection algorithms, based on hierarchical stochastic block modeling (hSBM), has been adapted to the topic modeling task [9].

we apply for the first time **hSBM-based topic modeling** to the study of cancer gene
expression data.





### Wang et al. 2021

We performed Lasso linear regression between histone acety-

lation sites and the protein and acetylation abundances of histone acetyltransferases (HATs), bromodomain-containing proteins (BRDs), and deacetylases (HDACs). It revealed potential connec- tions between HATs and BRDs and H2B acetylation sites, for example CREBBP and EP300, whose protein and acetylation levels showed substantial respective associations with H2B- K12, K13, K16, K17, and K21 and H2B-K21, and K24 sites (Fig- ure 5B).



f CausalPath (Babur et al., 2018) to the protein and phosphoprotein expression data (Figure S7A, Table S5)showed upregulation of the hypoxia pathway in mesenchymal tumors, evi- denced by significant activation of multiple HIF-1 downstream targets (networkpermutationp= 0.0012)



Using the Library of Integrated Network-Based Cellular Signa-

tures (LINCS) (Keenan et al., 2018; Stathias et al., 2019), we calculated the similarity between alteration-specific RNA or phosphoprotein signatures from our study with corresponding transcriptional (L1000 assay) (Subramanian et al., 2017) and phosphoproteomic LINCS signatures (P100 assay) (Litichevskiy et al., 2018) to identify compounds predicted to reverse tumor signatures of the cohort.



Multi-omics subtyping using non-negative matrix factorization (NMF) We selected the following proteogenomic features to the sample availability: CNV, bulk RNA, protein, and phosphoprotein expres- sion. Due to limited sample amounts, not all tumors were analyzed for DNA methylome, metabolome, and miRNA. We used non- negative matrix factorization (NMF) implemented in the NMF R-package (Gaujoux and Seoighe, 2010) to perform unsupervised clustering of tumor samples and to identify proteogenomic features that show characteristic expression patterns for each cluster. Briefly, given a factorization rank k (where k is the number of clusters), NMF decomposes a p 3 n data matrix V into two matricesWand H such that multiplication ofWand H approximates V. Matrix H is a k3n matrix whose entries represent weights for each sample (1 to N) to contribute to each cluster (1 to k), whereas matrixWis a p3 k matrix representing weights for each feature (1 to p) to contribute to each cluster (1 to k). Matrix Hwas used to assign samples to clusters by choosing the k with maximum score in each column of H. For each sample, we calculated a cluster membership score as the maximal fractional score of the corresponding column in matrix H. We defined a ‘‘cluster core’’ as the set of samples with cluster membership score > 0.5. Matrix Wcontaining the weights of each feature to a certain cluster was used to derive a list of representative features separating the clusters using the method proposed in (Kim and Park, 2007). To enable integrative multi-omics clustering we enforced all data types (and converted if necessary) to represent log2-ratios to
either a common reference measured in each TMT plex (proteome, phosphoproteome), an in silico common reference calculated as the median abundance across all samples (RNA gene expression) or to gene copy numbers relative to matching normal blood sample (CNV). All data tables were then concatenated and all rows containing missing values were removed. To remove uninforma- tive features from the dataset prior to NMF clustering, we removed features with the lowest standard deviation (bottom 5th percentile) across all samples. Each row in the data matrix was further scaled and standardized such that all features from different data types were represented as z-score



Unsupervised clustering of DNA methylation Methylation subtypes were segregated based on the top 8,000 most variable probes using k-means consensus clustering as previ- ously described (Sturm et al., 2012). We first removed underperforming probes (Zhou et al., 2017), and then the samples with more than 30% missing values. Remaining missing values were imputed using the mean of the corresponding probe value. We then per- formed clustering 1000 times using the ConsensusClusterPlus Rpackage



Unsupervised clustering of miRNA expression Unsupervised miRNA expression subtype identification was performed on mature miRNAs expression (log2 TPM) from 98 tumors with miRNA-seq available using Louvain clustering (Blondel et al., 2008) implemented in louvain-igraph v0.6.1. Top 50 differentially expressed miRNAs from each miRNA-based subtype were selected



Stemness scores were calculated as previously described (Malta et al., 2018). Firstly, we used MoonlightR (Colaprico et al., 2020)to query, download, and preprocess the pluripotent stem cell samples (ESC and iPSC) from the Progenitor Cell Biology Consortium (PCBC) dataset (D



We integrated somatic mutation, CNV, DNA methylation, RNA, protein, phosphorylation (phospho) and acetylation (acetyl) levels via iProFun (Song et al., 2019) to investigate the functional impacts of DNA alterations in GBM.



Weaggregated a set of interacting proteins (e.g. kinase/phosphatase-substrate or complex partners) from OmniPath (downloaded on 2018-03-29) (T€urei et al., 2016), DEPOD(downloadedon2018-03-29) (Duanet al., 2015),CORUM(downloadedon2018-06-29) (Ruepp et al., 2010), Signor2 (downloaded on 2018-10-29) (Perfetto et al., 2016), and Reactome (downloaded on 2018-11-01) (F

For each kinase-substrate protein pair supported by previous experimental evidence (OmniPath, NetworKIN, DEPOD, and SIGNOR), we tested the associations between all sufficiently detected phosphosites on the substrate and the kinase. For a kinase-substrate pair to be tested, we required both kinase protein/phosphoprotein expression and phosphosite phosphorylation to be observed in at least 20 samples in the respective datasets and the overlapped dataset. We then applied the linear regression model using lm function in R to test for the relation between kinase and substrate phosphosite.



Immune subtypes of the GBM tumors were generated on the consensus clustering of the cell type enrichment scores by



We investigated pathways from Hallmark, KEGG, WIkipathways, and REACTOME, positively or negatively aligned with averaged H2B and H3/H4 acetylation level. H2B acetylation was calculated by averaging acetylation of all H2B peptides detected. Since H3 and H4 histones are strongly correlated with each other, we averaged acetylation of histones H3 and H4 peptides together to obtain H3/H4 acetylation value



To test the association between HATs/HDACs protein and acetylation levels of histone sites, we fitted Lasso regression model with HATs/HDACs and histone protein expression as independent variables and a histone acetylation site as a dependent variable. Lasso



pathways from Hallmark, KEGG, WIkipathways, and REACTOME



We assumed that true biological activity of a pathway is regulated by collective changes of expression levels of majority of proteins involved in this pathway; then a difference in a pathway activity between tumors can be assessed by a difference in positioning of expression levels of proteins involved in this pathway in ranked list of expression levels of all proteins in each of tumors



Following
this idea, we assessed relative positioning of pathway proteins between tumor by determining two probabilities: (1) probability of
pathway proteins to occupy by random chance the observed positions in a list of tumor proteins ranked by expression level from
the top to the bottom and, similarly, (2) probability to occupy by random the observed positions in a list of expression levels ranked
from the bottom to the top. Then, the inferred relative activation of a given pathway across tumors was assessed as a negative log-arithm of the ratio of the above ‘‘top’’ and ‘‘bottom’’ probabilities. Thus, for a pathway of a single protein, its relative activity across
tumors was assessed as a negative log of ratio of two numbers: a number of proteins with expression level bigger than an expression
level of given protein, and a number of proteins with expression levels less than an expression level of given protein. For pathways of
multiple proteins, the ‘‘top’’ and ‘‘bottom’’ probabilities were computed as geometrically averaged P values computed for each of
proteins using Fisher’s exact test, given protein’s ranks in a list of pathway proteins and in a list of ranked proteins of a tumor, a num-ber of proteins in a pathway, and the total number of proteins with the assessed expression level in a given tumor. The thermody-namic interpretation of the inferred pathway activity scoring function is a free energy associated with deviation of the system from
equilibrium either as a result of activation or suppression. Thus, the scoring function is positive, when expression levels of pathway’s
proteins are overrepresented among top expressed proteins of a tumor, and it is negative, when pathway’s proteins are at the bottom
of expressed proteins of a tumor; the scoring function is close to zero, when expression levels are distributed by random. Given any
biological axis, e.g. histone acetylation levels in each of tumors, one can determine pathways which are significantly correlated or anti
correlated with the axis.



Causative pathway interaction discovery using CausalPath To discover the causative pathway interactions in our proteomic and phosphoproteomic data, we took the normalized expression of protein with < 10% missing values and phosphoprotein with < 25% missing values across all tumor and normal samples as the input to CausalPath (commit 7c5b934). We ran CausalPath in the mode that tests the mean values between test and control groups (value-transformation = significant-change-of-mean), where the test group being the tumors of one subtype and control group being the rest of the tumors. The pathway interaction discovery data source was Pathway Commons v9 (built-in- network-resource-selection = PC). Additionally, we enabled the causal reasoning if all the downstream targets of a gene were active or inactive (calculate-network-significance = true, use-network-significance-for-causal-reasoning = true, permutations- for-significance = 10000). 



### Satpathy et al. 2021

We performed non-negative matrix factorization (NMF)-based single- and multi-omic unsupervised clustering on CNA, RNA, protein, phosphoprotein, and acetylprotein datasets from 108 tumors, excluding ubiquitylprotein data as it was not available for the entire cohort (Figure 1A). The five resulting multi-omic subtypes (Figures 2A and 2B; Figure S2A) were named based on their predominant pathway associations and similarities to previously defined RNA clusters (Wilkerson

We used non-negative matrix factorization (NMF)-based multi-omic clustering using protein, phosphosite, acetylsite, RNA transcript and gene copy number variants (CNV) as previously described (Gillette

To enable integrative multi-omic clustering we enforced all data types (and converted if necessary) to represent ratios to: (i) a com-
mon reference measured in each TMT plex (proteome, phosphoproteome, acetylproteome), (ii) an in silico common reference calcu- lated as the median abundance across all samples (RNA expression) or (iii) to matching blood normal for CNA data. The CNA data was further filtered to only retain genes significantly altered (GISTIC2 thresholded of +2 or ?2) in at least 5% of all tumors. All data tables were then concatenated and all rows with missing values were removed. To remove uninformative features from the dataset prior to NMF clustering we removed features with the lowest standard deviation (bottom 5th percentile) across all samples. Each col- umn in the data matrix was further scaled and standardized such that all features from different data types were represented as z-scores. The resulting data matrix of z-scores into was converted to a non-negative input matrix required by NMF

For functional characterization of clustering results by single sample Gene Set Enrichment Analysis (ssGSEA), we calculated normalized enrichment scores (NES) of cancer-relevant gene sets by projecting the matrix of signed multi-omic feature weights (Wsigned) onto Hallmark pathway gene sets (Liberzon et al., 2015) using ssGSEA (Barbie et al., 2009). To derive a single weight for each gene measured across multiple omics data types (protein, RNA, phosphorylation site, acetylation site, CNA) we retained the weight with maximal absolute amplitude. We used the ssGSEA implementation available on https://github.com/broadinstitute/ ssGSEA2.0 using the following parameters:

The entire workflow described above has been implemented as a module for PANOPLY (Mani et al., 2020) which runs on Broad’s
Cloud platform Terra (https://app.terra.bio/). T

The CMAP (Lamb et al., 2006; Subramanian et al., 2017) is a collection of about 1.3 million gene expression profiles from cell lines treated with bioactive small molecules (?20,000 drug perturbagens), shRNA gene knockdowns (?4,300) and ectopic expression of genes. The CMAP dataset is available on GEO (Series GSE92742). For this analysis, we use the Level 5 (signatures from aggregating replicates) TouchStone dataset with 473,647 total profiles, containing 36,720 gene knock-down profiles, with measurements for 12,328 genes. See https://clue.io/GEO-guide for more information. 

Protein abundance comparisons were performed between all 5NMFsubtypes using the Wilcoxon rank-sum test and p values were adjusted using the Benjamini & Hochberg method. Subtype-specific differentially expressed proteins were identified based on their differential expression (adjusted p value < 0.05) in at least 2 out of the 4 comparisons, and by having a concordant fold change among all comparisons. The identified differentially expressed genes and proteins were then filtered for gene symbols measured in the L1000 assay (978 landmark genes). These NMF-specific signatures were used as input to calculate normalized weighted connectivity scores (WTCS) against the Library of Integrated Network-Based Cellular Signatures (LINCS) L1000 perturbation-response signa- tures. The scores were computed using the sig_queryl1k_tool pipeline (https://hub.docker.com/u/cmap) and the LINCS L1000 Level 5 compound (trt_cp) signatures from CLUE (https://clue.io, ‘‘Expanded CMap LINCS Resource 2020 Release’’). 



a Library of Integrated Network-Based Cellular Signatures- based (LINCS) query for compounds that reversed the EMT-E signature showed enrichment for TGFb inhibitors (

Normalized connectivity scores (NCS) between the Library of Integrated Network-based Cellular Signatures (LINCS) compound signatures and the five NMF multi-omic subtype signatures at the RNA and protein level are represented in the heatmap. Each NMF signature was used to query the CMAP LINCS L1000 resource (2020 release), a collection of more than 1M response signatures of cancer cell lines to different types of chemical (e.g., drugs, tool compounds) and genetic perturbations (e.g shRNAs, sgRNAs, cDNAs). The 20 most negatively connected compounds to each NMF signature (with the most negative NCS) were then selected for visualization. 



Immunohistochemistry-based antibody-specific staining scores in lung tumors were obtained from the Human Protein Atlas (HPA, https://www.proteinatlas.org), in which tumor-specific staining is reported in four levels, i.e., high, medium, low, and not detected. The protein-specific annotations such as enzymes, transcription factor, transporters, secreted, transmembrane and FDA- approved drugs targeting the protein or reported in drugbank were designated



Continuous smoking score Non-negative matrix factorization (NMF) was used in deciphering mutation signatures in cancer somatic mutations stratified by 96 base substitutions in tri-nucleotide sequence contexts. To obtain a reliable signature profile, we used SomaticWrapper to call mu- tations from WGS data. SignatureAnalyzer exploited the Bayesian variant of the NMF algorithm and enabled an inference for the optimal number of signatures from the data itself at a balance between the data fidelity (likelihood) and the model complexity (reg- ularization). After decomposing into three signatures, signatures were compared against known signatures derived from COSMIC (Tate et al., 2019) and cosine similarity was calculated to identify the best match. We also sought to integrate count of total mutations, t, percentage that are signature mutations, c, and count of DNPs, n, into a continuous score, 0 < S < 1, to quantify the degree of confidence that a sample was associated with smoking signature. We referred to these quantities as the data, namely D=CC¸TC¸ N, and used Aand A’ to indicate smoking signature or lack thereof, respectively. In a Bayesian framework, it is readily shown that a suitable form is S = 1 / (1 + R), where R is the ratio of the joint probability of A’ and D to the joint probability of A and D. For example, the latter can be written P(A)・P(CjA)・P(TjA)・P(NjA) and the former similarly, where each term of the former is the complement of its respective term in this expression. Common risk statistics are invoked as priors, i.e., P(A) = 0.9 (W



Ranking tumors by inferred activity of IFN-g pathway We assumed that true biological activity of a pathway is regulated by collective changes of expression levels of majority of proteins involved in this pathway. Then, a difference of a pathway activity between tumors can be assessed by a difference in positioning of expression levels of proteins involved in this pathway in ranked list of expression levels of all proteins in each of tumors. Following



ssGSEA (Barbie et al., 2009) was utilized to obtain pathway scores based on RNA-seq and global prote- omics data using the R package GSVA (Ha¨nzelmann et al., 2013). For this analysis, pathways from the Reactome (Fabregat et al., 2018), KEGG (Kanehisa et al., 2017) and Hallmark (Liberzon et al., 2015) databases were

Independent component analysis ICA was performed with a workflow modified from previously described (Liu et al., 2019). Decomposition was run for 100 times on the matrix of protein abundance difference between tumor/NAT pairs (n = 99). Independent components were in the form of vectors comprised with weights of all genes in the original data. Components extracted from each run were clustered using HDBSCAN al- gorithm (McInnes et al., 2017) with cosine distance as dissimilarity metric, m

Correlation between the extracted signatures and known clinical characteristics were examined by regressing the corresponding mixing scores for all members of a component cluster against 64 sample annotations to obtain within-cluster average of log10 p values





known interactome DBs including Omnipath, Phos- phositeplus, DEPOD, Signor, and CORUM. A



To identify germline genetic variants that explain variation in tumor gene (eQTL) and protein (pQTL) expression, we utilized the gold- standard mapping pipeline at https://github.com/molgenis/systemsgenetics/wiki/eQTL-mapping-analysis-cookbook-(eQTLGen)





We performed phosphosite-specific signature enrichment analysis (PTM-SEA) (Krug et al., 2019) to identify dysregulated phosphor-ylation-driven pathways.



We queried the PTM signatures database (PTMsigDB) v1.9.0 downloaded from https://github.com/broadinstitute/ssGSEA2.0/ tree/master/db/ptmsigdb using the flanking amino acid sequence (+/? 7 aa) as primary identifier. We used the implementation of PTM-SEA available on GitHub (https://github.com/broadinstitute/ssGSEA2.0) using the command interface R-script (ssgsea- cli.R). The



Association analysis between KGG-site abundances and E3 ligases and DUBs A list of known human E3 ubiquitin and ubiquitin-like ligases and DUBs was compiled from (Medvar et al., 2016; Nijman et al., 2005). We then fit a linear model using limma in R with the formula kgg_site_abundance ?protein_abundance, followed by empirical bayes shrinkage



Genes in these site-wise clusters were used for pathway enrichment analysis against the KEGG, Reactome, and WikiPathways databases using g:profiler (Raudvere et al., 2019). Pathway enrichment was performed using a gene background containing all observable genes in the K-GG dataset. Pathway enrichment results were imported into Cytoscape (Shannon et al., 2003) using the Enrichment Map app (Merico et al., 2010) for network analysis of pathways. Pathways were connected using the gene overlap and clustered (pathway cutoff q-val < 0.1; jaccard overlap > 0.375). Each cluster was manually annotated from the pathways con- tained in it to facilitate interpretation



DepMap genetic dependency dataset (CRISPR Avana Public 20Q3) that contained 18119 genes and 789 cell lines (https://depmap.org/portal/download/ file: Achilles_gene_effect.csv). O

CausalPath (Babur et al., 2018) searches for known biological mechanisms that can explain correlated proteomic changes in terms of causal hypotheses

OncoKB database for oncoprotein and tumor suppressor classification (excluded proteins that have both annotations), and used PhosphoSitePlus and Signor databases for the activating/inhibiting classification of phosphorylation sites. 



NeoFlow (https://github.com/bzhanglab/neoflow) for neoantigen prediction (Wen et al., 2020). Specifically, Optitype (Szolek et al., 2014) was used to find human leukocyte antigens (HLA) for each sample based on WES data. Then we used netMHCpan (Jurtz et al., 2017) to predict HLA peptide binding affinity for somatic mutation–derived variant peptides with a length between 8-11 amino acids. The



We used Customprodbj (Wen et al., 2020)(https://github.com/bzhanglab/customprodbj) for customized database construction.



Remaining variant peptides were further filtered using PepQuery (http://www.pepquery.org)(Wen et al., 2019) with the p value cutoff of 0.01. Competitive filtering based on unrestricted posttranslational modification searching was enabled in PepQuery validation. The spectra of variant peptides were annotated using PDV (http://pdv.zhang-lab.org)(Li et al., 2019b



PROGENy (Schubert et al., 2018) was used to generate activity scores for EGFR based on RNA expression data.



##### Babur et al. 2010

GEM (Gene Expression Modulation) is a probabilistic framework that predicts modulators, their affected targets and mode of action by combining gene ex- pression profiles, protein–protein interactions and transcription factor–target relationships.

Wang et al. (5) propose MINDy, an information- theoretic algorithm for detecting modulators. They test the conditional mutual information (CMI) between the transcription factor and the target gene, and its depend- ency on the modulator candidate. This is, in essence, the aforementioned non-linearity principle. Building upon the same principle, we present GEM (Gene Expression Modulation), a probabilistic method for detecting modulators of transcription factors using a priori know- ledge and gene expression profiles. For a modulator/tran- scription factor/target triplet, GEM predicts how a modulator–factor interaction will affect the expression of the target gene. GEM improves over MINDy by detecting two new classes of interaction that would result in strong correlation but low ?CMI, can filter out logical-or cases and offers a more precise classification scheme. A detailed comparison of GEM and MINDy is provided in the discussion



##### Babur et al. 2015

a novel method for the identification of sets of mutually exclusive gene alterations in a given set of genomic profiles. We scan the groups of genes with a common downstream effect on the signaling network, using a mutual exclusivity criterion that ensures that each gene in the group significantly contributes to the mutual exclusivity pattern. We test the method on all available TCGA cancer genomics datasets, and detect multiple previously unreported alterations that show significant mutual exclusivity and are likely to be driver events

We are expanding on these approaches by combining
detailed prior pathway information with a novel statistical metric to improve both accuracy and biological inter- pretation and to validate the results. Specifically, we are using a large aggregated pathway model of human sig- naling processes to search groups of mutually exclusively altered genes that have a common downstream event. To



##### Robertson et al. 2017

Groups of samples with similar abundance profiles were identified by unsupervised consensus clustering with ConsensusClusterPlus (CCP) 1.20.0. Calculations were performed using Spearman correlations, partitioning around medoids (PAM) and 10,000 iterations



PARADIGM Integrated Pathway Analysis Integrated Pathway Levels (IPLs) mRNA expression, SCNA, and pathway interaction data for 80 UM samples were integrated using the PARADIGM software (Sedge- wick et al., 2013). Briefly, this procedure infers integrated pathway levels (IPLs) for genes, complexes, and processes, using pathway interactions, and genomic and functional genomic data from each patient sample. 

Pathways were obtained in BioPax Level 3 format, and included the NCIPID and BioCarta databases from http://pid.nci.nih.gov
and the Reactome database from http://reactome.org.

Differential pathway regulators of each PARADIGM clusters were identified by comparing one cluster vs. all others using the t-test
and Wilcoxon Rank sum test with a BH FDR correction. All

Hierarchical MARINa (hMARINa) Estimate of Kinase Activity MARINa is well suited for the analysis of TF activity, because TF proteins are directly involved in changes in expression of their tar- gets



statistically significant findings from the PARADIGM and (h)MARINa differential pathway regulator analyses were examined for consistency. For each cluster, pathway regulators with similar findings across the two methods were identified as ‘‘consistent pathway features.’’

 MARINa was run via the VIPER R package (http://www.bioconductor.org/packages/release/bioc/html/viper.html)(Alvarez et al., 2016); and hMARINa was per- formed by extending the functionality of the package



PhosphositePlus and the SuperPathway (see Curated TF Regulome and Curated Kinase Regulome sections)



To identify genes that were associated with time to metastasis in M3 cases, we censored time and status at 5 years for the 32 Lau-
rent and 33 TCGA records.

Creating a Curated Transcription Factor (TF) Regulome A compendium of TFs and their targets (TF regulons) were created by combining information from four databases:
(i) SuperPathway (Sedgewick et al., 2013): This is the same interaction network used in the PARADIGM analysis (above). Only links that correspond to regulation at the transcriptional level were retained for MARINa and hMARINa use.
(ii) Literome (Poon et al., 2014): The network was filtered to include only transcription links in which the regulator is a known TF. (iii) Multinet (Khurana et al., 2013): The network was reduced to links that correspond to regulation on transcriptional level. (iv) ChEA (Lachmann et al., 2010): Data from the Gene Expression Atlas (Petryszak et al., 2014) was used to filter the inferred links in the ChEA database. Specifically, the context likelihood of relatedness (CLR) method (Faith et al., 2007) was used to compute a measure of association between every pair of genes. The top 10% of gene pairs based on the CLR score were intersected with the ChEA network and the overlapping pairs were added to the final combined network.
The combined network includes 72,915 transcriptional regulatory links between 6,735 regulators and their targets. Only regulators
with at least 15 targets were considered in the final analysis, which resulted in a final network consisting of 419 TFs with 58,363 total targets (covering a set of 12,754 unique targets). Creating a Curated Kinase Regulome Proteins identified as kinases in Manning (Manning et al., 2002) or Uniprot (UniProt Consortium, 2014) were aggregated into a list of 546 kinases. Protein substrates were extracted from PhosphositePlus (Hornbeck et al., 2014) on March 7, 2015. Kinase-substrate interactions were retained if the kinase appeared in the Manning-Uniprot kinase list and the kinase was identified as a human protein in the PhosphositePlus database. The final compendium consisted of 5,388 links between 342 kinases and 2,260 unique substrates.



Brett et al. 2021

Omic and Multidimensional Spatial (OMS) Atlas generated from four serial biopsies of a metastatic breast cancer patient during 3.5 years of therapy. This

the HTAN Data Coordinating Center (https://humantumoratlas.org/) 

aggregation of publicly available molecular interactions and biological pathway databases provided by the Pathway Commons (PC) resource.23 The aggregated data is represented in the standard Biological Pathway Exchange (BioPAX) language and provides the most complete and rich representation of the biological network models stored in PC. These complex biochemical reactions were reduced to pairwise relationships using rules to generate a Simple Interaction Format (SIF) representation of BioPAX interactions. The reduction of BioPAX interactions to the SIF allows for the representation of pairwise molecular interactions in the context of specific binary relationships. The



Regulon enrichment: This method leverages pathway information and gene expression data to produce regulon-based protein activity scores. Our method tests for positional shifts in experimental-evidence supported networks consisting of transcription factors and their downstream signaling pathways when projected onto a rank-sorted gene-expression signature. The gene-expression signature is derived by comparing all features to the median expression level of all samples considered within the data-set. After weights have been assigned to the regulatory network, the positive and negative edges of each regulator are rank ordered. The first component of the enrichment signature, the local delta concordance signature, is derived by capturing the concordance between the magnitude of the weight assigned to a particular edge and the ranked position of that edge. The features associated with activation, positive edges within the regulatory network, are monotonically ranked from most lowly to highly expressed in the restricted feature space, where the features that are repressed are ranked by a monotonically decreasing function. This component of the signature considers positive and negative edges independently, which captures support for an enrichment signature even if one of the edge groups is underrepresented in the network graph. The second component of the enrichment signature, the local enrichment signature, captures positional shifts in the local gene ranked signature and integrates those shifts with the weights assigned to overlapping features for a given regulon and the expression data set. The last component of the enrichment signature considers the entire feature space and projects the rank-sorted local signature onto this global ranked feature space. We derive our global enrichment signature from this projection for each regulator we consider. We use the median of robust quantile-transformed ranked positions as the enrichment scores for both the local enrichment and global enrichment signatures. We then integrate these three individual signatures together, which allows us to capture differences between individual regulator signatures within the context of an individual patient as well as at a cohort level. Reverse



CausalPath (https://github.com/PathwayAndDataAnalysis/causalpath; commit 9f8d6f8) was used for integrated pathway analysis of protein, phosphoprotein, gene abundance, and transcriptional regulator activity.24 Briefly, CausalPath is a hypothesis generating tool that uses literature- grounded interactions from Pathway Commons to produce a graphical representation of causal relationships that are consistent with patterns in a multi-omic datasets.24 This integrative approach allows for holistic evaluation of signaling networks and pathway activity across longitudinal biopsies.



### Boizard et al. 2021

a new prioritization strategy called PRYNT (PRioritization bY protein NeTwork) that employs a combination of two closeness‑based algorithms, shortest‑path and random walk, and a contextualized protein–protein interaction (PPI) network, mainly based on clique consolidation of STRING network.



Other approaches use biological networks in order to prioritize candidates (e.g. MaxLink21, ToppNet18). One of the network-based software most commonly used by biologists in order to interpret high-throughput expression data is Ingenu- ity Pathway Analysis (IPA)22. This suite is based on a PPI network containing millions of structured, manually curated experimental observations. In IPA, the “Upstream Regulator Analysis” (URA) algorithm prioritizes disease candidates using in-house causal network approach to elucidate upstream biological causes that can explain the observed molecular changes23,2



One of the main limitations hampering the use of IPA is that the software is proprietary and therefore its use cannot be broadly generalized to the biology community

PRYNT is based on the integration of Search Tool for the Retrieval of Interacting (STRING, version 10.5)29 PPI network and a combination of shortest-path and random walk, two closeness-based algorithms as it has been previously shown in the literature that this method outperformed other computational methods1





##### Rodchenkov et al. 2020

Pathway Commons (https://www.pathwaycommons. org) is an integrated resource of publicly available information about biological pathways including bio- chemical reactions, assembly of biomolecular com- plexes, transport and catalysis events and physical interactions involving proteins, DNA, RNA, and small molecules 

Data is collected from multiple providers in stan- dard formats, including the Biological Pathway Ex- change (BioPAX) language and the Proteomics Stan- dards Initiative Molecular Interactions format, and then integrated. Pathway

Pathway Commons provides biolo- gists with (i) tools to search this comprehensive re- source, (ii) a download site offering integrated bulk sets of pathway data (e.g. tables of interactions and gene sets), (iii) reusable software libraries for work- ing with pathway information in several program- ming languages (Java, R, Python and Javascript) and (iv) a web service for programmatically query- ing the entire dataset

Pathway Commons currently con- tains data from 22 databases with 4794 detailed hu- man biochemical processes (i.e. pathways) and ∼2.3 million interactions

www.pathguide.org)



##### Babür et al. 2020

We used CausalPath to map site-specific changes in protein phosphorylation in our data to known, causal relations curated through Pathway Commons.25 We recently developed CausalPath to mimic how a scientist might search literature sources to understand changes in protein phosphorylation within a quantitative dataset – but, on the scale of hundreds of thousands of events – and, with the computational capacity to logically infer how thousands of phosphorylation changes may have occurred. CausalPath

**Babür et al. 2008**

We have built a microarray data analysis tool, named PATIKAmad, which can be used to associate microarray data with the pathway models in mechanistic detail, and provides facilities for visualization, clustering, querying, and navigation ofbiological graphs related with loaded microarray experiments. PATIKAmad is freely available to noncommercial users as a new module ofPATIKAweb at http://web.patika.org

d PATIKAmad, within PATIKAweb [5], which is a Web interface to the PATIKA database for querying, visu- alizing, and analyzing biological networks. 



##### Touré et al. 2020

A large variety of molecular interactions occurs between biomolecular components in cells. When a mo- lecular interaction results in a regulatory effect, exerted by one component onto a downstream component, a so- called ‘causal interaction’ takes place. Causal interactions constitute the building blocks in our understanding of larger regulatory networks in cells. These

we propose a checklist that accommodates current representations, called the Minimum Information about a Molecular Interaction CAusal STatement (MI2CAST). 



##### Kim et al. 2021

a comprehensive analysis using multi-omics data has been conducted. In addition, a pathway activity inference method has been developed to facilitate the inte- grative effects of multiple genes.

a novel integrative pathway activity in- ference approach, iDRW and demonstrated the effectiveness of the method with respect to dichotomizing two sur- vival groups.

The reason why we inferred pathway activities is that our approach aims at not only integrating multi-omics data based on the graph, but also a pathway-level representation of multiple genomic profiles to better analyze the prioritized pathways considering interactions between genes and to improve outcome prediction performance using ma- chine learning models

we propose a general framework
for integrative pathway activity inference on the multi-omics net- work and investigate multiple network scenarios

To reflect the interaction effects of genes, we designed a directed gene–gene graph in multiple layers by assigning within-layer interactions and
between-layer interactions considering multiple scenarios. We inferred pathway activities by performing a random walk with re- start (RWR) on the multi-layered network. 

As a result, iDRW trans- forms the multiple genomic profiles into a single pathway profile on the graph. 

The inferred pathway profile is validated with the out- come prediction models. We prioritized pathways, visualized the multi-omics network and extensively analyzed the pathway activity patterns.

iDRW jointly prioritizes potential driver pathways and genes on multi-omics data

##### ràv

A breast carcinoma context-specific network model of transcriptional  regulation was assembled with the ARACNe, based on 851 RNA-seq  expression profiles obtained from TCGA. ARACNe was run with 100  bootstraps, a p value threshold of 10−8, and 0 data  processing inequality (DPI) tolerance, generating a network of 1,748 TFs associated with 18,783 target genes by 459,569 interactions. The  regulatory models were generated from the ARACNe results using the VIPER package from Bioconductor (http://bioconductor.org/packages/release/bioc/html/viper.html).

The gene expression signatures for 20 ER+ and 11 TNBC metastases (MET-GES)  were computed with paired Student’s t test by comparing their profiles  against the matching primary tumor ones. Then, the enrichment of each [regulatory protein](https://www.sciencedirect.com/topics/neuroscience/regulator-protein) [regulon](https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/regulon) on the MET-GESs was inferred by the VIPER algorithm ([Alvarez et al., 2016](https://www.sciencedirect.com/science/article/pii/S2211124717310318#bib2), [Aytes et al., 2014](https://www.sciencedirect.com/science/article/pii/S2211124717310318#bib5)), as implemented in the VIPER package for R available from Bioconductor (https://www.bioconductor.org/packages/release/bioc/html/viper.html). Statistical significance was estimated by permuting the samples uniformly at random 1,000 times.

For single patient-based analysis, gene expression signatures were computed by comparing each MET expression profile with the matching primary  tumor expression profile. A null model for statistical testing was  generated by permuting the samples uniformly at random 1,000 times. The  most consistent MRs across all tumors were prioritized according to the  average dysregulation (normalized enrichment score [NES]) divided by its SD.



##### Ietswaart et al. 2021

GeneWalk (github.com/churchmanlab/genewalk) that identifies individual genes and their relevant functions critical for the experimental setting under examination. After the automatic assembly of an experiment-specific gene regulatory network, GeneWalk uses representation learning to quantify the similarity between vector representations of each gene and its GO annotations, yielding annotation significance scores that reflect the experimental context.

new methods are required to generate functional information about individual genes under particular conditions of interest or biological contexts. To address this need, we developed GeneWalk, a knowledge-based machine learning and statistical modeling method that highlights the gene functions that are relevant for a specific biological context

deep learning to condense information [33–36], and generation of gene networks de- rived from database aggregation efforts

GeneWalk is developed to generate functional relevance information about individual
genes in a biological context under study. 

GeneWalk first automatically assembles a biological network from a knowledge base and the GO ontology starting with a list of genes of interest (e.g., differentially expressed genes or hits from a genetic screen) as in- put

The network structure is learned through random walks using an un- supervised network representation learning algorithm

The resultant vector representations enable a quantitative comparison between genes and GO terms, highlighting the GO terms most relevant for the biological context under study

, GeneWalk provides for each input gene its direct GO annotations ranked by their statistical relevance. 



##### Broyde et al. 2021

Oncoprotein-specific molecular interaction maps (SigMaps) for cancer network analyses

Although net- works derived from pairwise interaction assays or computational inference might mitigate the excessive simplicity and linearity of cancer pathways, they generally do not account for nor discrimi- nate between cellular contexts1

we developed an integrative machine
learning (ML) framework (**OncoSig**) for the systematic, de novo reconstruction of tumor-specific molecular interaction signaling maps (**SigMaps**), anchored on any oncoprotein of interest. Specifically,

an oncoprotein-specific SigMap recapitulates the molecular architecture necessary to func- tionally modulate and mediate its activity within a specific cellular context, including its physical, cognate-binding partners. To

OncoSig infers context-specific SigMaps by train- ing an ML algorithm to integrate complementary evidence from transcriptional and post-translational interactions inferred from three-dimensional (3D) structural data, as well as from gene expres- sion and mutational profiles in large-scale repositories

four complementary evidence sources integrated by OncoSig. Additional evi- dence can be easily incorporated in the framework. We

**example: KRAS context**

1. KRAS-specific, structure-based protein–
   protein interactions (PPIs), as inferred by the Predicting Protein– Protein Interactions (PrePPI) algorithm6,7
   , by combining struc-
   tural homology to protein complexes in structural databases and non-structure-related data. PrePPI scores represent the likelihood ratio (LR) of any predicted PPI based on a random interaction null model
2. include transcriptional interactions inferred by the
   Algorithm for the Reconstruction of Gene Regulatory Networks (ARACNe)
3. Third, we use Virtual Inference of Protein activity by Enriched Regulon analysis (VIPER)10
   to associate the mutational state of a
   candidate upstream modulator protein with differential activity of the anchor protein or the mutational state of the latter with differential activity of its candidate downstream effectors. VIPER
   measures a protein’s activity based on the expression of its transcrip- tional targets—akin to a tissue-specific, highly multiplexed gene reporter assay.
4. Finally, we infer upstream modulators of the anchor protein using Conditional Inference of Network Dynamics (CINDy)11,12
   —a
   refinement of the Modulator Inference by Network Dynamics (MINDy) algorithm11
   . CINDy uses the conditional mutual infor-
   mation to assess changes in the mutual information between the anchor protein and its transcriptional targets as a function of the candidate modulator expression or mutational state

OncoSig accounts for tumor context specificity by lever-
aging evidence from the ARACNe, VIPER and CINDy algo- rithms, whose predictions are based on the analysis of large-scale, tumor-specific molecular profile data, whereas PrePPI provides context-independent, structure-based evidence.

To integrate the evidence from these algorithms, we tested two established ML algorithms: Naive Bayes14
(NB) and Random Forest1

An advantage of the former is that the inference of specific protein–protein relationships can be easily traced back to their sup- porting evidence. In contrast, the latter is better suited to integrate non-statistically independent evidence sources. However, tracing back predictions to the specific supporting evidence is challenging

input to the algorithm: matrix
(Fig. 1c), with ~20,000 rows—one for each protein in the human proteome.

The gold standard set (GSS) vector in the second col- umn describes proteins known to be functionally related (PGSS) or unrelated (NGSS) to the anchor protein

Remaining columns represent ~36,000 features corresponding to PrePPI, ARACNe, CINDy and VIPER confidence scores, respectively, supporting physical or functional interactions between row-specific and column-specific proteins

. For each protein, statistically sig- nificant algorithm scores are reported in different columns, hence accounting for the greater number of columns than rows

The entire matrix, except for the gold standard column, is identi-
cal for any anchor protein of choice

In addition to the 61 true-positive predictions, SigMapKRAS
LUAD
includes 201 novel predictions at false discovery rate (FDR) = 1%
I
(red circles) (Fig. 2b(ii)), including 30 proteins predicted as KRAS physical interactors (11%, blue and black circles) in BioGRID21
,
134 druggable proteins (51%, red and black circles) in the Drug Repurposing Hub22
and 33 proteins meeting both criteria

SigMaps effectively recapitulate known biology and can help priori- tize novel functional or physical interactors, including many drug- gable ones, for validation

1)PrePPI6
predicts interactions between a protein (red) and its physical and/or functional interactors (gray). 2) The
ARACNe algorithm9
predicts transcription factors or signaling molecules (red) that transcriptionally regulate target genes (blue). 3) CINDy12
predicts
signaling molecules (orange/red) that post-translationally modify transcription factors (blue boxes), which, in turn, leads to differential expression of a transcription factor’s targets (blue diamonds). 4) The VIPeR algorithm10
infers downstream effectors (blue) and upstream regulators (orange) for a given
protein (red). VIPeR associates i) the protein (red) with a mis-sense mutation (black dot) with the activity change of transcription factors (blue) and ii) signaling molecules (orange) with mis-sense mutations (black dots) with activity of the protein (red). c,



we developed an unsupervised version of the algorithm (OncoSigUN
) to extend the analysis to arbi-
trary proteins of interest, without protein-specific training sets.



The term ‘pathway,’ although widely used, is a loosely defined bio- logical concept. Here we propose a fundamentally different rep- resentation (SigMap) of the signaling and regulatory machinery necessary to modulate and affect the function of a specific protein of interest in a specific tissue context, which is equivalent to a protein’s mechanism of action. 



Our data suggest that SigMaps provide a more unbiased, com-
pact and realistic representation of a protein’s mechanism of action, compared to available network representations and algorithms

OncoSig generates a single integrated score repre-
senting the probability that a protein belongs to a specific SigMap. Use of PrePPI is instrumental for identifying PPIs, whereas ARACNe, VIPER and CINDy provide critical tissue specificity and additional evidence supporting both physical and functional interactions

the feature matrix (Fig. 1c) was reduced to contain only interactions with the specific protein of interest, leaving only four of the 36,000 columns, one for each of the algorithms

Proteins were then scored based on aggregate voting across the ten OncoSigRF
classifiers described

The rationale is that, once a
sufficient number of diverse training sets is available, they can be used to assess the generic contribution of each evidence source (that is, its weight) toward classification of a bona fide interaction.

The features derived from the networks were as follows: mutual information
for ARACNe, number of statistically significant triplets for CINDy, negative log P value for VIPER and LR for PrePPI. We coded each feature symmetrically, so that interactions between protein A and protein B were input into the matrix twice, once in the feature vector for A and once in the feature vector for B; all other elements in the RF feature matrix were set to zero

For each of the ten oncogene-centric interactomes, proteins that are part of the PGSS were assigned a ‘1’ within the PGSS vector, whereas all other proteins were assigned a ‘0’ to represent membership in the negative gold standard set (NGSS). OncoSigRF

OncoSigRF was trained and tested with each pathway’s PGSS and NGSS using Monte Carlo cross-validation77
, creating 50 forests each with 50 trees

General procedure for SigMaps and OncoSig. The feature matrix for a protein of interest consists only of those features that correspond to interactions with the protein of interest and, thus, has only four columns, one for each of the four interactomes: PrePPI, ARACNe, CINDy and VIPER. SigMap membership of each protein in the human proteome is determined by aggregate voting after querying the ten OncoSigRF
classifiers described above (KRAS, PI3KCA, TP53, EGFR, BRAF,
STK11, CDKN2A, NTRK3, YAP1 and CTNNB1).

##### Mall et al. 2021 Network-based identification of key master regulators associated with an immune-silent cancer phenotype

The availability of large genomic datasets offers an opportunity to ascertain key determinants of differential intratumoral immune response. 

We follow a network-based protocol to identify transcription regulators (TRs) associated with poor immunologic antitumor activity. 

We use a consensus of four different pipelines to identify TRs affecting immunologic antitumor activity

* 2 state-of-the-art gene regulatory network inference techniques to determine TR regulons
  1. regularized gradient boosting machines 
  2. ARACNE ,
* 3 separate enrichment techniques, 
  1. fast gene set enrichment analysis 
  2. gene set variation analysis 
  3. virtual inference of protein activity by enriched regulon analysis to identify the most important TRs affecting immunologic antitumor activity. 

**through an open-science competition (DREAM Chal- lenge), the authors compared various GRN inference methods on several synthetic and real datasets. In [26], the authors illus- trated the superior performance of RGBM for the DREAM Chal- lenge networks (see Supplementary Figure S1b). Hence, RGBM is the primary GRN inference technique focused on in this work**

Another key component of MRA is to estimate enrichmen-
t/activity scores for TRs in a given sample, taking into consid- eration its regulon.

While techniques such as RGBM utilize a simplistic difference in average expression of positively and negatively regulated targets to estimate the activity of a TR, methods such as virtual inference of protein-activity by enriched regulon analysis (VIPER) [8] and MARINA [6] utilize a dedicated algorithm formulated to estimate TR activity tak- ing into account the TR mode of action, the TR-target gene interaction confidence and the pleiotropic nature of each target gene regulation. Moreover,

**, there exists single sample gene set enrichment analysis [57] techniques such as gene set varia- tion analysis (GSVA [27]) and fast gene set enrichment analysis (FGSEA [53]) to estimate enrichment score for each TR in a given sample.** This

Finally,we perform downstream analysis of the MRs specific to ICR-L using ConsensusPathDB [35] to discover corresponding enriched pathways

Previously tools such as ARACNE, VIPER and NetFactor focused on tran- scription factors (TFs). Moreover, the original RGBM algorithm also exploited an active binding network based on binding sites of TFs. Recently, there have been studies [14, 47, 48]thatextend the hubs of GRNs to regulatory proteins beyond the TFs. For example, in [48], the authors considered a set of 2506 regulatory proteins annotated in GO with TF activity and transcription cofactor activity**. Our set of 3674 TRs was a superset of their set including receptors, kinases, growth factors, signal transduction proteins, transcription co-activators and cofactors as candidate regulators**

Other interesting examples where hubs of the networks were
focused on signal molecules (and not just TFs) include approach such as SigMaps [14] or surface receptors i.e. the receptors inter- actome, to identify active ligand-receptors pairs [47]. These stud- ies used ARACNE + VIPER and generalized the concept of MRA to generic signal molecules (not just TFs) as originally intended

To determine the enrichment score with statistical significance for specific TR regulons, we use the ‘fgsea’ function in the ‘fgsea’ package in R

### Mechanism-Centric Approaches for Biomarker Detection and Precision Therapeutics in Cancer

Yu and Mitrofanova 2021

e computational approaches that identify **mechanism-centric** biomarkers elucidated from gene co-expression networks, regulatory networks (e.g., transcriptional regulation), protein–protein interaction (PPI) networks, and molecular pathways.

the majority of biomarkers are identified from **gene-centric approaches** (we will refer to gene/protein/metabolite etc.,-centric approaches as gene-centric approaches for simplicity), where either a specific gene is investigated (based on previous biological assumptions) or a gene(s) is selected based on differential behavior without connection to the upstream and downstream molecular mechanisms. Gene-centric

Gene-centric findings are often limited in mechanistic interpretability and connectivity to other molecular processes, positioning such biomarkers as passengers, rather than drivers, of the biological process and thus are **often dataset specific**

in **white-box** models (e.g., linear regression and decision trees) the relationship between input variables (i.e., genes) and output variables (i.e., disease outcomes) is understandable/explainable as they often identify linear or monotonic relationships (Zhang

**black-box** models (e.g., neural networks, gradient boosting, or ensemble models such as random forest) are able to capture non-linear/non-monotonic relationships, yet often suffer from model interpretability and subsequent limited clinical adoption (Wang

they mostly capture associative relationships when applied as gene-centric approaches and often miss the complexity of mechanisms inherent in biological systems, especially in the context of cancer.

**mechanism-centric approaches**, which are not focused on single genes and take into account complex mechanisms implicated in cancer initiation, progression, and treatment response. In

<u>Gene Co-expression Network Analysis</u>

* Network Construction: 
  * WGCNA and 
  * lmQCM
* Network Mining: 
  * Centered Concordance Index,
  * Eigengenes,
  * Hubs

**WGCNA** calculates correlation between pairs of genes and transforms the correlation measure into a topological overlap measure in order to minimize effects of noise and spurious associations. The

The resulting matrix is subjected to hierarchical clustering to determine groups of co-expressed genes,

genes cannot be assigned to multiple modules, exposing WGCNA’s limitation since many genes participate in multiple biological processes and often perform multiple functions. An

lmQCM algorithm identifies densely connected subnetworks (i.e., quasi-cliques) using a greedy search algorithm which allows module overlaps

can also identify smaller modules, which can highlight more specific and interpretable biological connections as compared to much larger modules of WGCNA that frequently contain over a thousand genes

Co-expression networks can be mined to determine the functional significance of their modules or identify functionally relevant genes. Here,

Centered Concordance Index has been developed to identify
modules specific to each condition/phenotype.

**CCI** evaluates the concordance of gene expression profiles within a module based on singular value decomposition and is used to identify modules that are highly co-expressed in one condition over another (Han

The CCI is useful in identifying modules specific to phenotype conditions but has yet to be used to associate modules with continuous outcomes.

The **eigengene** approach transforms modules into weighted
vectors, which mathematically correspond to their contribution to the first principal component in principal component analysis

Eigengenes are then able to be associated with clinical features (including continuous outcomes) using correlation/association measures

The translational applicability ofmodules can be hampered by
their relatively large size and might benefit from **identification of hub genes** within modules:

**intramodular connectivity** for gene i is defined as the sum of edge weights between gene i and the other genes in the module (Zhang

**betweenness centrality**, which is a network topology metric used to identify
central nodes in a graph based on a shortest paths algorithm

The betweenness centrality of gene i is a measure of the number of shortest paths connecting any two genes which pass through i.

<u>Regulatory Network Analysis</u>

* Transcriptional Regulatory Networks

  * network construction
    * ARACNe
  * network mining
    * MARINa
    * VIPER

* Multi-Omic Regulatory Network

  * network construction
    * RegNetDriver (step 1)
  * network mining
    * RegNetDriver (step 2)

  

One of the most known and widely experimentally validated methods for transcriptional network reconstruction is **ARACNe**

This information-theoretic algorithm utilizes tissue-specific gene expression profiles to estimate pairwise mutual information between expression levels of TFs/co-TFs and expression levels of their potential (activated or repressed) targets. The

The advantage of using mutual information to measure such relationships lies in its ability to measure not only linear (which would be captured for example by the Pearson correlation) or monotonic (which would be captured for example by Spearman correlation) relationships, but also non- linear associations

Data processing inequality results in a regulatory network that includes primarily direct TF/co-TF-target interactions

The ARACNe network can be effectively interrogated (i.e., mined) using MARINa (Lefebvre et al., 2010) and VIPER two algorithms that identify TFs/co- TFs as driver biomarkers associated with specific phenotypes

Specifically, **MARINa** (Lim et al., 2009; Lefebvre et al., 2010) requires a differentially expressed signature, defined as a ranked list of genes between any two phenotypes of interest. Then, the activated and repressed targets for each TF/co-TF (as inferred by ARACNe) are assessed for their enrichment in the over- and under-expressed parts of this signature (Lefebvre et al., 2010; Figure 3). 

Such enrichment is referred to as TF/co- TF transcriptional activity, and if it is statistically significant, the TF/co-TF is referred to as a Master Regulator (MR). As a result of this analysis, a TF/co-TF is considered an “activated” MR if its activated targets are significantly enriched in the over-expressed part of the signature and/or its repressed targets are significantly enriched in the under-expressed part of the signature. Conversely, a “repressed” MR exhibits the opposite behavior. 

It is important to note that TF/co-TF transcriptional activity is not defined based on **the differential expression of TFs/co-TFs themselves but instead on the differential expression of their transcriptional targets**. This allows the identification of TFs/co-TFs that are not necessarily differentially expressed but are modified on the post-translational level and would otherwise be missed by traditional association methods.

At the same time, **VIPER** estimates TF/co-TF transcriptional
activity **on an individual sample-based level**, as opposed to a two-phenotype signature-based level required by MARINa (Alvarez et al., 2016; Figure 3). In fact, while MARINa requires carefully selected multiple samples of the same phenotype to construct a differential expression signature, VIPER is able to utilize **single-sample analysis by scaling the overall patient cohort (to its average expression for each gene).** Furthermore, several advantages of VIPER include estimation of TF/co- TF activity through a so-called **mode of regulation** (taking into account whether targets are activated, repressed, or their direction cannot be determined), inference of regulator-target **interaction confidence**, and accounting for target overlap between different regulators (Alvarez et al., 2016). 

**RegNetDriver** is an algorithm for multi-omic tissue-specific regulatory network construction and analysis (Dhingra

The regulatory network reconstructed by RegNetDriver represents a **two-layered relationship**: (i) connecting TFs to promoter/enhancer regions; and (ii) further connecting promoter/enhancer regions to their corresponding target genes.

To reconstruct relationships between TFs and promoters/enhancers of potential targets, Dhingra et al. utilize tissue-specific (i.e., prostate epithelium) DNase I hypersensitive sites to define accessible regulatory DNA regions and integrate this information with promoter/enhancer annotations from ENCODE (Encode Project Consortium, 2012) and GENCODE

TFs are then connected to promoters/enhancers based on the enrichment of their binding motifs. Promoters/enhancers

Promoters/enhancers are further connected to their target genes through significant correlation of promoter/enhancer region activity signals (ChIP-seq signals) with target gene expression profiles 

This network is then utilized to **identify TF hubs** with genomic and epigenomic alterations that can potentially cause large perturbations in this tissue-specific network. Specifically, TFs are first mined on degree centrality, such that the top 25% of TFs with the greatest number of outgoing edges are defined as hubs. Next, to identify TF hubs significantly affected on genomic and epigenomic levels in prostate cancer, they are evaluated for the presence of prostate-cancer specific genomic alterations (single nucleotide variants and structural variants) and DNA methylation changes in their coding and non-coding regulatory regions. 

<u>Protein–Protein Interaction Network-Based Analysis</u>

* Network construction
  * Chuang et al., Step 1
* network mining
  * Chuang et al., Step 2

a hybrid approach to combine a PPI network with tissue-specific gene expression profiles across patient samples. The PPI network is comprised of nodes representing proteins and edges representing a characterized PPI, utilizing subnetworks from CellCircuits. Tissue-specific gene expression data are then overlaid onto all PPI subnetworks. For each subnetwork, its **activity in each sample/patient is defined as a combination of z-scores for the subnetwork genes**. This defines patient-specific vectors of subnetwork activities, which are then mined for phenotype associations.

Activities of subnetworks are evaluated for their **association with specific phenotypes** (e.g., metastatic and non-metastatic), where associations can be calculated by mutual information, t-score, or Wilcoxon score and is referred to as the **subnetwork discriminative potential/score**. Next, the method selects subnetworks with a locally maximal discriminative score and performs significance testing to ensure subnetworks are non- random and robust. In

<u>Pathway-Based Analysis</u> 

* pathCHEMO 
* pathER

**Pathways** represent a group of biochemical entities (e.g.,
genes, proteins, etc.), connected by interactions, relations, and reactions (including physical interactions, complex formation, transcriptional regulation, etc.), that lead to a certain product or changes in a cell. Molecular

**pathCHEMO** was specifically developed to compare poor versus good therapeutic response (as categorical outcomes) in cancer. 

 it evaluates differential behavior of biological pathways on both transcriptomic (RNA expression) and epigenomic (DNA methylation) levels between any two phenotypes of interest (Epsi et al., 2019). 

First, an RNA expression treatment response signature is defined as a list of genes ranked by their differential expression between poor and good treatment response. 

Then, genes in each pathway are evaluated for their enrichment in either over-expressed, under-expressed, or differentially expressed (which includes both over- and under- expressed) part of this signature. 

Enrichment in the over- and under-expressed parts separately allows identification of pathways where the majority of genes exhibit a similar behavior (i.e., are either over- or under-expressed), while enrichment in the differentially expressed part of the signature allows identification of pathways where some genes are over-expressed and some are under-expressed (which depicts a complex interplay ofactivation and repression relationships inside a molecular pathway). This enrichment is referred to as the **RNA expression-based activity level** of a molecular pathway.

**DNA methylation-based activity** for each pathway is estimated in the same manner using a DNA methylation treatment response signature

Pathways that are enriched in the RNA expression treatment response signature and the DNA methylation treatment response signature are then integrated to select those that are significantly affected on both expression and methylation levels (Figure

Activity levels of the candidate pathways are further evaluated as biomarkers of therapeutic response in independent patient cohorts

**pathER** applies a pathway-
based approach on a **single-patient level**, which allows the association of pathway activity across a patient cohort to a wide range of therapeutic responses

a **multivariable regression Cox proportional hazards** model to associate pathway activity levels with time-to-therapeutic failure, thus capturing poor, good, and medium therapeutic responses.

### DGCA: A comprehensive R package for Differential Gene Correlation Analysis - McKenzie et al. 2016

offers a suite of tools for **computing and analyzing differential correlations between gene pairs across multiple conditions**

minimize parametric assumptions, DGCA computes **empirical p-values via permutation testing.** To

To understand differential correlations at a systems level, DGCA performs **higher-order analyses** such as measuring the average difference in correlation and multiscale clustering analysis of differential correlation networks

Distinct from differential expression, **differential correlation** operates on the level of gene pairs rather than individual genes (Fig.

While differential coexpression analysis has proven useful
in identifying significantly different modular connectivity patterns, **differential correlation analysis of individual gene pairs** is far more granular. As

Like DiffCorr, DGCA **transforms correlation coefficients to z-scores** and uses **differences in z-scores to calculate p-values** of differential correlation between genes. Like

Like Discordant, DGCA **classifies differentially correlated gene pairs into the nine possible categories**. However,

DGCA differs from the existing differential correlation approaches 

1. DGCA calculates false discovery rate of differential correlations through non-parametric sample permutation.
2.  DGCA can calculate the average difference in correlation between one gene and a gene set across two conditions
3.  DCGA integrates with MEGENA to perform multiscale clustering analysis of differential correlation networks to identify gene mod- ules (clusters) and hub genes
4. DGCA provides comprehensive downstream functional analysis of differen- tial correlation structures including visualization, gene ontology (GO) enrichment, and network tool

DGCA has 3 main inputs 

1. a matrix of **gene expression values**, 
2. a **design matrix** specifying conditions associated with samples, 
3. a **specification of the condi- tions** for comparison (Fig.

**central tendency** re- fers to measures of centrality in a distribution, including the arithmetic mean or median,

**dispersion** refers to measures of spread in a distribution, including the standard deviation and the dispersion index

To **stabilize the variance** of sample correlation coefficients in each condition, the **Fisher z-transformation**

**Fisher z-transformation function serves as a normalizing transformation. The**

Using the **difference in z-scores** dz, a **two-sided p-value** can be calculated using the standard normal distribution. Gene pairs can then be ranked on the basis of their relative strength of differential correlation.
Multiple

several options for adjusting p-values for multiple hypothesis tests, including the conservative Benjamini-Hochberg p-value adjustment method [25, 26] and the local false discovery rate method

difficult to make intuitive sense of the p-values returned because the p-values are originally derived from the difference of z-scores method, which depends on specifying the correct form for the variance of the sample correlation coefficients, andinturnonthebivariatedistributionofthegene expres- sion values. Therefore,

DGCA also offers to generate **per- mutation samples by randomly shuffling the sample labels** across the input conditions and then re-computing the dif- ferential correlation calls

The z-scores from the original and permuted data sets are used to calculate **empirical p-values**, using a reference pool distribution approach

These empirical p-values are used to estimate the proportion of null hypoth- eses in empirical p-values by extrapolating a linear trend from a cubic spline fitted over candidate ranges of the tun- ing parameter lambda [29].

q-values are calculated based on the empirical p-values and the estimated propor- tion of null hypotheses.

At the most basic level, gene pairs can be classified as hav- ing **gain of correlation (GOC) or loss of correlation (LOC)** between one condition compared to another. For

go beyond this binary classification, we also determine if two genes are **significantly correlated in each condition or not**. By

Based upon a **threshold for correlation significance and the sign of correlation in each condition** (i.e., positive or negative), gene-gene correlations in each condition can be **categorized into 3 classes, i.e. significant positive correlation, no significant correlation, and signifi- cant negative correlation**. 

there are **9 classes for differential correlations** between two conditions (Fig.

difference in average correlations of a gene and a set of genes between two conditions.

DGCA **quantifies the median difference in z-transformed correlation coefficients of a gene and a gene set** (henceforth, median difference in z-score) between two conditions. In

a median difference in z-score above 0 indicates a tendency towards a gain of correlation between the given gene and the gene set in the first condi- tion with respect to the second condition, while a median difference in z-score below 0 indicates a tendency towards a loss of correlation. To

To measure the significance of the median change in correlation permu- tation samples to calculate empirical p-values,

DGCA also offers to calcu- late the median difference in z-score between all gene pairs in two conditions.

This approach is similar to the **modular** **differential correlation** calculation that

### CEMiTool: a Bioconductor package for performing comprehensive modular co- expression analyses - Russo et al. 2018

CEMiTool that unifies the **discovery and the analysis of co-expression modules**. 

outperforms existing tools, and provides unique results in a user-friendly html report with high quality graphs

evaluates **whether modules contain genes that are over-represented by specific pathways** or that are **altered in a specific sample group**,

**integrates transcriptomic data with interactome information**, identifying the potential hubs on each network.

WGCNA: 

* following its workflow verbatim is time-consuming and tiresome
* users are often required to manually select parameters and to filter the input genes prior running WGCNA
* limited in terms of the func- tional analyses available for the package users

integrating co-expression infor- mation with protein-protein interaction data can be useful to identify main regulators or hubs.

(CEMiTool), an R package that allows users to easily identify and analyze co-expression modules in a fully automated manner

provides users with 

* a novel unsupervised gene filtering method, 
* automated parameter selection for identifying modules, 
* enrichment and module func- tional analyses, 
*  integration with interactome data. 
*  reports everything in HTML web pages with high-quality plots and interactive tables.

automating within a single R function (cemitool) the entire module discov- ery process - including gene filtering and functional analyses

, gene set enrichment analysis (**GSEA**) [6] can **asso- ciate the activity of a module with the study phenotypes** (i.e. sample group)

. Over-representation analysis (**ORA**) can be used to reveal if a set of co-expressed genes is enriched for genes belonging to known pathways or functions. In

workflow

* gene expression file containing the genes as rows and the samples as columns. This file is the only required input for CEMiTool’s analyses.
* unsupervised filtering method based on the inverse gamma distribution to select the genes used in the analyses.

* soft-thresholding power β used to determine a similarity criterion between pairs of genes
* genes are then separated into modules using the Dynamic Tree Cut pack- age

* If an optional file containing **gene interactions** (e.g. protein-protein interaction data) is provided, the package will return network graphs composed of **interact- ing genes within the same module.**
* if the user provides a sample **annotation file**, perform gene set enrichment analysis (**GSEA**), allowing users to visualize **which modules are induced or repressed in the different phenotypes**. Finally,
* given an optional file containing **gene sets**, CEMiTool will perform an over rep- resentation analysis (**ORA**) based on the hypergeometric test to determine the **most significant module functions**

over representation analysis (ORA) via the clusterProfiler R package [8].

sample annotation file describing the phenotypes (i.e. disease, healthy, treated, etc) of samples, CEMiTool performs a gene set enrichment analysis using the fgsea (Fast Gene Set Enrichment Analysis) R package; genes from co-expression modules will be treated as gene sets and the z-score normalized expression of the samples within each phenotype will be treated as rankings on the analysis; assess if the activity of a module is altered across dif- ferent phenotypes.

provide a gene interaction file to visualize the interactions between the genes in each co-expression module. This

customize their module graphs according to different interaction databases

Although WGCNA provides a function named pickSoft- Threshold that can automatically select the soft-threshold β value, we have created an alternative algorithm, which is based on the concept of Cauchy sequences [17], that improves the automatic selection of this parameter, allow- ing for more reliable and consistent results. Moreover, our algorithm allows a lower threshold for R2 (R2 >0.8) when compared to WGCNA’s default threshold (R2 >0.85). This, in turn, allows for lower values of β**.Oncea β value is chosen, subsequent steps for creating modules follow standard WGCNA procedure**

### DysPIA: A Novel Dysregulated Pathway Identification Analysis Method - Wang et al. 2021

We adopted the idea of Correlation by Individual Level Product into analysis and performed a fast enrichment analysis.

most of the DC methods relied on the Pearson correlation coefficient (PCC)

The **Correlation by Individual Level Product** (CILP) was proposed in Lea et al. (2019) article to **identify factors associated with interindividual variation in correlation**

can be used to estimate the dysregulated status of each gene pair between case and control samples. In

a unique score for each gene pair in each sample, rather than a summary PCC statistic for a group.

These conventional methods ofpathway analysis focused on gene marginal effects in a pathway and ignored gene interactions that may contribute to a phenotype of interest. For



* first generation methods: Classical pathway enrichment analysis methods, such as DAVID, based on overrepresented statistical tests (such as Fisher’s exact test and hypergeometric test) to assess whether DE/DV genes were overrepresented in a predefined pathway.

* second-generation methods:  ranking all genes according to their DE levels, and then used the weighted Kolmogorov–Smirnov statistic to test whether genes from a prespecified pathway were significantly overrepresented toward the top or bottom of the ranked gene list (GSEA and later ones)

* third-generation: the pathway topology was incorporated into the analysis. However, this gene regulation information (edge) was just used to adjust the gene (node) value in these methods

  

For two genes in a pathway, neither of them may have an effect on a phenotype of interest. However, when they were jointly considered, they may have a significant effect on the studied phenotype due to the gene–gene interaction. In

the measures of these methods are mainly based on individual gene levels, they can be deemed as node (gene)- centric methods

the gene-pair relationships have not been fully considered.

differential co-expression (DC)-based pathway analysis aimed to identify pathways with more gene regulation differences related to the phenotype ofinterest

we proposed a novel method called
Dysregulated Pathway Identification Analysis (**DysPIA**) to overcome these shortcomings. 

**A pathway is represented by the regulated gene pairs, but not just a set of genes,** which is used in the traditional pathway analysis. 

We adopted the idea of **CILP** into analysis and performed a **fast GSEA-like enrichment analysis**.

1.  calculate a **Dysregulated Gene Pair Score** (DysGPS) for each gene pair of interest (**individual-level-based statistic** instead of population-level). 
   * z-score normalization
   * for each gene-pair Xi and Xj from CILP-like statistic (product of these two genes’ standardized expression values)
   * dysregulated gene-pair score (DysGPS) calculated based on a two-sample Welch’s t-test between groups
2. calculate the **Dysregulated Pathway Score** (DysPS) based on a GSEA-like formula (pre-ranked pathway enrichment analysis pipeline)
   * Rank the N gene pairs in the combined-background set in descending order based on DysGPS to form a gene- pair list
   * Calculate DysPS for each pathway
     * gene pairs get different weights based
       on whether they are in the pathway. 
   * Randomly permute the sample labels, recalculate the DysGPS in the background, and recalculate the DysPS for each pathway
3.  permutation-based significant **P-values** estimated and the BH-FDR adjustment performed. 



### Advantages of CEMiTool for gene co-expression analysis of RNA-seq data - Cheng et al. 2020

we evaluate three co-expression analysis packages (WGCNA, CEMiTool, and coseq) using published RNA-seq datasets derived from ischemic cardiomyopathy and chronic obstructive pulmonary disease. Results show that the packages produced consensus co-expression clusters using default parameters. CEMiTool package outperformed the other two packages and required less computational resource and bioinformatics experience.

groups of co-expressed genes are clustered based on several methods, including hierarchical clustering and K-means clus- tering, which are discussed in detail elsewhere [19].

WGCNA and CEMiTool are both based on hierarchical clustering, while coseq uses K-means clustering. WGCNA

WGCNA and CEMiTool packages are similar in principle, but the latter provides an automated pipeline. CEMiTool

*WGCNA*

WGCNA was used to construct signed weighted gene co-expression
modules from the top 5000 variable genes

Signed correlations cluster positively and negatively correlated genes into modules.

Module detection was based on the hi- erarchical clustering of adjacencies given by the topological overlap measure. A soft thresholding power β -value

. A series of β -values (ranging between 1 and 30) was screened to evaluate the average connectivity degrees of different modules. A β -value was selected by plotting the R2 against soft threshold β. 

the adjacency matrix was transformed into a topological overlap matrix (TOM) to measure the connectivity of genes within the network

The connectivity of genes is defined as the sum of its adjacency in relation to all other genes in the network

TOM is measured between a value of 0 and 1

Based on TOM, a higher value (towards 1) indicates the set of genes are highly connected, and therefore, the strong interconnectivity will create a meaningful co-expression association

Contrary, when the TOM value is closer to 0, it signifies no connections between genes. A minimum module size was set to 40, highly correlated modules were merged by setting merging modules threshold to 0.2, and

*CEMiTool*

CEMiTool is similar to WGCNA but runs on an automated pipeline.

CEMiTool CEMiTool is similar to WGCNA but runs on an automated pipeline.

CEMiTool implements auto- mated gene filtration on gene expression profiles based on the inverse gamma distribution.

### WGCNA: An R package for weighted correlation network analysis - Langfelder 2008

To calculate the adjacency matrix, an intermediate quantity called the co-expression similarity sij is first defined. The default method defines the co- expression similarity sij as the absolute value of the corre- lation coefficient between the profiles

signed co-expression measure can be
defined to keep track of the sign of the co-expression information.

Using a thresholding procedure, the co-expression simi- larity is transformed into the adjacency. An

A weighed net- work adjacency can be defined by raising the co-expres- sion similarity to a power

weighted adjacency aij between two genes is proportional to their similarity on a logarithmic scale,

Adjacency functions for both weighted and unweighted networks require the user to choose threshold parameters, for example by applying the approximate scale-free topology criterion

Several measures of network interconnectedness

the topological overlap measure

WGCNA identifies gene modules using unsupervised clus- tering, i.e. without the use of a priori defined gene sets. The user has a choice of several module detection meth- ods. The default method is hierarchical clustering

One drawback of hierarchical clustering is that it can be difficult to determine how many (if any) clusters are present in the data set

A co-expression module may reflect a true biological signal (e.g. a pathway) or it may reflect noise (e.g. a technical artifacts, tissue contamina- tion, or a false positive). To test whether the identified modules are biologically meaningful, gene ontology information (functional enrichment analysis) can be used. Toward this end, we provide an R tutorial that describes how to interface the WGCNA package with rele- vant external software packages and databases.

### BioNERO: an all-in-one R/Bioconductor package for comprehensive and easy biological network reconstruction - Almeida-Silva et al. 2021

BioNERO, an R package that aims to integrate all aspects of network analysis workflows, 6 including expression data preprocessing, gene coexpression and regulatory network 7 inference, functional analyses, and intra and interspecies network comparisons

Users can infer three types of GCNs (signed, signed hybrid or unsigned), and 58 pairwise gene-gene correlations can be calculated with Pearson’s r, Spearman’s ρ, or biweight 59 midcorrelation (median-based, which is less sensible to outliers).

we implemented three widely used GRN inference algorithms: 74 GENIE3 (Huynh-Thu et al., 2010), ARACNE (Margolin et al., 2006), and CLR (Faith et al., 75 2007). However, choosing the most appropriate number of top edges to keep is a persisting 76 bottleneck, and users often pick an arbitrary number. We implemented a method to simulate 77 different networks by splitting the graph in n subgraphs, each containing the top nth quantiles

two network comparison features in BioNERO, 

1.  consensus module 84 identification and module preservation. **Consensus modules are gene modules that co-occur 85 in networks inferred from independent expression sets**, and they can be used to explore core 86 components of the studied phenotype that are not affected by experimental effects or natural 87 biological variation. 
2. module preservation focuses on the differences, and it can be used to 89 explore patterns of transcriptional divergence within and across species. For

### LeapR: An R Package for Multiomic Pathway Analysis - Danna et al. 2021

leapR package, a framework to rapidly assess biological pathway activity using diverse statistical tests and data sources, allowing facile integration of multisource data

one generally needs to evaluate numerous approaches requiring a unified framework. Furthermore, none support data integration using phosphoproteomics data.	

simplify pathway analysis of multiple different types of
data, including post-translational modification, we have developed a framework, the layered enrichment analysis of pathways in R (leapR), to represent multiple omics types, perform pathway analysis on the individual sets or combined sets, and analyze and represent the results in a biologically meaningful manner.

Multiomic Pathway Enrichment

our approach focuses on a role-agnostic approach to data integration, which treats the different types of data the same in the enrichment process. Though

### Comparing Statistical Tests for Differential Network Analysis of Gene Modules - Arbet et al. 2021

Differential network analysis (DiNA)

how this network structure differs between two or more groups/phenotypes (

3 common applications of DiNA involve

1. testing whether the connections to a single gene differ between groups, 
2. testing whether the connection between a pair of genes differs between groups, or 
3. testing whether the connections within a “module” (a subset of 3 or more genes) differs between groups. 

e focuses on the latter, as there is a lack of studies comparing statistical methods for identifying differentially co-expressed modules (DCMs). Through

The DGCA R package (McKenzie et al., 2016) simply uses the mean (or median) of the differences. A potential problem with this approach is that positive and negative differences can cancel out, thus losing power to detect DCMs where some correlations increase while other correlations decrease between conditions

Many test statistics can be formulated as functions of the difference (or product) in V(gm) between the two groups. For example, the “Dispersion Index” (DI), used by GSCA (Choi and Kendziorski, 2009) and DiffCoEx (Tesson

GSNCA can still be used to test whether the network structure of a module differs between the two groups. Briefly, GSNCA assigns a weight vector to each group of length the absolute differences of the weight vector between the two |Pm| (one weight per gene) and the test statistic is the sum of groups. The ith gene is given a weight wi that is proportional to the sum of the correlations between the ith gene with all other genes. Thus, a gene that is highly correlated with many other genes will be given a larger weight, which indicates the gene may have regulatory importance

The motivation of the PND test is, given a “partially differentially co-expressed module” (a module where some of the correlations, but not all, change between groups), then the higher the exponent p, the less weight is given to the null correlations that do not change between groups. Therefore, we expect the PND test with a large value of p (e.g., p ≥ 4) to be more sensitive for detecting DCMs where only a small proportion of the module correlations change between conditions.

In contrast to unconditional correlation, TOM captures shared relationships or “neighbors” between the two genes,

The intuition behind TOM is that if the two genes i and j are connected to a common set of genes, then the similarity between the two genes, S(gm)
ij , should increase (i.e., the greater
the number and strength of the connections that are shared by genes i and j, the larger the TOM value will be for those two genes).

To calculate TOM, one must first calculate the correlation matrix, then convert to an adjacency matrix (aij), and then calculate the TOM matrix

“signed” versions since we want to be able to detect correlations that change from positive to negative between groups, when



### Mechanism-Centric Approaches for Biomarker Detection and Precision Therapeutics in Cancer - Yu et al. 2021

MECHANISM-CENTRIC COMPUTATIONAL APPROACHES FOR BIOMARKER DISCOVERY

* Gene Co-expression Network Analysis

*Network Construction: WGCNA and lmQCM*

* WGCNA
  * genes cannot be assigned to multiple modules, exposing WGCNA’s limitation since many genes participate in multiple biological processes and often perform multiple functions. An

* lmQCM
  * m identifies densely connected subnetworks (i.e., quasi-cliques) using a greedy search algorithm which allows module overlaps (Ou
  * can also identify smaller modules, which can highlight more specific and interpretable biological connections as compared to much larger modules of WGCNA that frequently contain over a thousand genes

*Network Mining: Centered Concordance Index, Eigengenes, and Hubs*

* CCI
* eigengene
* intramodular connectivity
* betweeness centrality

*Regulatory Network Analysis*

Transcriptional Regulatory Networks

* Network construction: ARACNe
  * information-theoretic algorithm utilizes tissue-specific gene expression profiles to estimate pairwise mutual information between expression levels of TFs/co-TFs and expression levels of their potential (activated or repressed) targets
  * ability to measure not only linear (which would be captured for example by the Pearson correlation) or monotonic (which would be captured for example by Spearman correlation) relationships, but also non- linear associations. Another
  * data processing inequality, which eliminates any “indirect” regulatory relationship through the principle that mutual information on the indirect path cannot exceed mutual information on any part of the direct path
* Network mining: MARINa and VIPER The
  *  ARACNe network can be effectively interrogated (i.e., mined) using MARINa (Lefebvre et al., 2010) and VIPER (Alvarez et al., 2016), two algorithms that identify TFs/co- TFs as driver biomarkers associated with specific phenotypes (e.g.,
  * MARINa (Lim et al., 2009; Lefebvre et al., 2010) requires a differentially expressed signature, defined as a ranked list of genes between any two phenotypes of interest. Then, the activated and repressed targets for each TF/co-TF (as inferred by ARACNe) are assessed for their enrichment in the over- and under-expressed parts of this signature (Lefebvre et al., 2010; Figure 3). Such enrichment is referred to as TF/co- TF transcriptional activity, and if it is statistically significant, the TF/co-TF is referred to as a Master Regulator (MR). As
  * VIPER estimates TF/co-TF transcriptional
    activity on an individual sample-based level, as opposed to a two-phenotype signature-based level required by MARIN
  * MARINa requires carefully selected multiple samples of the same phenotype to construct a differential expression signature, VIPER is able to utilize single-sample analysis by scaling the overall patient cohort (to its average expression for each gene).
  * several advantages of VIPER include estimation of TF/co- TF activity through a so-called mode of regulation (taking into account whether targets are activated, repressed, or their direction cannot be determined), inference of regulator-target interaction confidence, and accounting for target overlap between different regulators (Alvarez

Multi-Omic Regulatory Network

* RegNetDriver

Protein–Protein Interaction Network-Based Analysis

* Network construction and mining Chuang et al.
  * hybrid approach to combine a PPI network with tissue-specific gene expression profiles across patient samples.
  * Tissue-specific gene expression data are then overlaid onto all PPI subnetworks. For each subnetwork, its activity in each sample/patient is defined as a combination of z-scores for the subnetwork genes. This defines patient-specific vectors of subnetwork activities, which are then mined for phenotype associations.
  * Activities of subnetworks are evaluated for their association with specific phenotypes (e.g., metastatic and non-metastatic), where associations can be calculated by mutual information, t-score, or Wilcoxon score and is referred to as the subnetwork discriminative potential/score. Next, the method selects subnetworks with a locally maximal discriminative score and performs significance testing to ensure subnetworks are non- random and robust. In
* Pathway-Based Analysis: pathCHEMO and pathER
  * pathCHEMO was specifically developed to compare poor versus good therapeutic response (as categorical outcomes) in cancer
    *  it evaluates differential behavior of biological pathways on both transcriptomic (RNA expression) and epigenomic (DNA methylation) levels between any two phenotypes of interest 
    * an RNA expression treatment response signature is defined as a list of genes ranked by their differential expression between poor and good treatment response
    *  genes in each pathway are evaluated for their enrichment in either over-expressed, under-expressed, or differentially expressed (which includes both over- and under- expressed) part of this signature. 
    * Enrichment in the over- and under-expressed parts separately allows identification of pathways where the majority of genes exhibit a similar behavior (i.e., are either over- or under-expressed), while enrichment in the differentially expressed part of the signature allows identification of pathways where some genes are over-expressed and some are under-expressed (which depicts a complex interplay ofactivation and repression relationships inside a molecular pathway). 
    * this enrichment is referred to as the RNA expression-based activity level of a molecular pathway. 
    * DNA methylation-based activity for each pathway is estimated in the same manner using a DNA methylation treatment response signature. 
    * Pathways that are enriched in the RNA expression treatment response signature and the DNA methylation treatment response signature are then integrated to select those that are significantly affected on both expression and methylation levels (Figure 6). 
  * pathER applies a pathway-
    based approach on a single-patient level, 
    *  allows the association of pathway activity across a patient cohort to a wide range of therapeutic responses 
    * multivariable regression Cox proportional hazards model to associate pathway activity levels with time-to-therapeutic failure, thus capturing poor, good, and medium therapeutic responses. 



### Network-based classification of breast cancer metastasis - Chuang et al. 2007

we apply a protein- network-based approach that identifies markers not as individual genes but as subnetworks extracted from protein interaction databases.

subnetwork markers are more reproducible than individual marker genes selected without network information, and that they achieve higher accuracy in the classification of metastatic versus non-metastatic tumor

a more effective means of marker identifica- tion may be to combine gene expression measurements over groups of genes that fall within common pathways. Several approaches have been proposed to score known pathways by the coherency of expression changes among their member genes (Pavlidis

a remaining hurdle to pathway-based analysis is that the majority of human genes have not yet been assigned to a definitive pathway

a number of approaches have been demonstrated for extracting relevant subnetworks based on coherent expression patterns of their genes (Ideker et al, 2002; Chen and Yuan, 2006) or on conservation of subnetworks across multiple species

Each subnetwork is suggestive of a distinct functional pathway or complex

identifying markers of metastasis within gene expression profile

not encoded as individual genes or proteins, but as subnetworks of interacting proteins within a larger human protein–protein interaction network

although genes with known breast cancer mutations are typically not detected through analysis of differential expression, such as P53, KRAS, HRAS, HER-2/neu, and PIK3CA, they play a central role in the protein network by interconnecting many expression-responsive genes. Third,

To integrate the expression and network data sets, we
overlaid the expression values of each gene on its correspond- ing protein in the network and searched for subnetworks whose activities across the patients were highly discriminative of metastasis.

a candidate subnetwork was first scored to assess its activity in each patient, defined by averaging its normalized gene expression values. This

Second, the discriminative potential of a candi- date subnetwork was computed based on the mutual information between its activity score and the metastatic/ non-metastatic disease status over all patients. 

Significantly discriminative subnetworks were identified by comparing their discriminative potentials to those of random networks

A subnetwork is defined as a gene set that induces a single connected component in the protein–protein interaction network

Given a particular subnetwork M, let a represent its vector of activity scores over the tumor samples, and let c represent the corresponding vector of class labels (metastatic or non-metastatic). To

To derive a, expression values gij are normalized to z-transformed score

The individual zij of each membergene in the subnetwork are averaged into a combined z-score, which is designated the activity aj. Many

score the relationship between a and c. In this study, we define the discriminative score S(M) as MI(a0,c), the mutual information MI between a0, a discretized form of a, and c

To derive a0 from a, activity levels are discretized into log 2ð# of samplesÞþ 1bc¼ 9 equally spaced bins

A rationale for using MI in cancer classification is to capture potential heterogeneity of expression in cancer patients (Tomlins et al, 2005), that is, differences not only in the mean but in the variance of expression

Given the discriminative score function S, a greedy search is performed to identify subnetworks within the protein–protein interaction net- work for which the scores are locally maximal.

Candidate subnetworks are seeded with a single protein and iteratively expanded. At each iteration, the search considers addition of a protein from the neighbors of proteins in the current subnetwork and within a specified network distance d from the seed. The addition that yields the maximal score increase is adopted; the search stops when no addition increases the score over a specified improvement rate r

To assess the significance of the identified subnetworks, three tests
of significance are performed

For the first test, we perform the same search procedure over 100 random trials in which the expression vectors of individual genes are randomly permuted on the network

The second test indexes each real subnetwork score on a ‘local’ null distribution, estimated from the scores of 100 random subnetworks

. Third, we test whether the mutual information with the disease class is stronger than that obtained with random assignments of classes to patients (



### Convolutional neural network for human cancer types prediction by integrating protein interaction networks and omics data - Chuang et al. 2021

. Here, we collected 6136 human samples from 11 cancer types, and integrated their gene expression profiles and protein–protein interaction (PPI) network to generate 2D images with spectral clustering method. To predict normal samples and 11 cancer tumor types, the images of these 6136 human cancer network were separated into training and validation dataset to develop convolutional neural network (CNN).

Recently, Teppei et al. combined two kinds of biological data, the gene expression profile and human PPI network, to generate 2D representation as the input of the spectral-CNN model17

we integrated PPI network and gene expression profile of 11 cancer types, to generate 6136 network images in 2D representation by using spectral clustering (i.e., Laplacian matrix). Where 1228 network images were used for training and testing in CNN model; and the other 4908 images, gene expression clustering and survival data were used for validation

https:// github. com/ bioxg em/ CNN_ model. git

The human PPIs were collected from five public databases (i.e., BioGRID20, DIP21, IntAct22, MINT23 and
MIPS24), including 16,433 human proteins and 181,868 PPIs. To combine RNA-Seq and PPIs data, we assigned proteins with gene expression using gene name and gene ID, and finally acquired 14,230 proteins and 152,519 PPIs for further analysis.

We first iden- tified differentially expressed genes (DEGs) between tumors and corresponding normal tissues for 11 cancer types by computing gene expression fold change and modified t-statistic (limma package v.3.38.3). Finally, 12,024 genes were considered as DEGs with |fold change| ≥ 2 and adjust p value < 0.01 in at least one cancer type. By these DEGs, we selected a maximum-subnetwork with 6261 DEGs with 28,439 PPIs and combined the gene expression profiles of 5528 tumors and 608 normal tissues. These cancer networks will be processed with dimensionality reduction utilizing spectral clustering, for cancer prediction and classification in CNN model.

we used a spectral cluster- ing approach, Laplacian (L) matrix to reduce dimensionality of complex cancer networks and applied on CNN techniques. The

the matrix cells were assigned value “− 1” when
the two proteins had interaction, otherwise the cells were assigned value “0”; whereas the cells in diagonal were assigned value of node degree (the number of edges connected to the node in the network). Next, we obtained the eigenvalue and eigenvector of Laplacian matrix using linear transformation. To retain the network topology and connectivity, we utilized the smallest and second smallest non-negative and non-zero eigenvalues with their corresponding eigenvectors to map the cancer network (6261 DEGs and 28,439 interactions) into 2D spaces with 100 × 100 cells (Fig. 2B)27,28. After the dimensionality reduction for PPI network, the 1849 unique nodes were displayed in 2D representation and assigned with gene expression value of clinical samples (if numerous genes overlapped into a single node, then their gene expression was averaged and assigned to the node). In total, we generated 6136 images of cancer networks for CNN model to predict tumors, normal tissues and cancer types



### Investigating the relevance of major signaling pathways in cancer survival using a biologically meaningful deep learning model - Feng et al  2021

, deep learning models have recently been proposed in survival prediction, which directly integrates multi-omics data of a large number of genes using the fully connected dense deep neural network layers, which are hard to interpret

investigate potential associations between patient survival and individual signaling pathways, which can help domain experts to understand deep learning models making specific predictions



we built a simplified and partially biologically meaningful deep neural network, DeepSigSurvNet, for survival prediction. In the model, the gene expres- sion and copy number data of 1967 genes from 46 major signaling pathways were integrated in the model. We



The cox proportional hazards model (Cox PH) model [1] is the classic model for survival analysis. The Kaplan–Meier estimator curve [2], CoxPH model and logrank test [3] are widely used to display and compare the survival probability over time of patients in different groups or conditions. 

the interpret- able analysis identified the distinct patterns of these signaling pathways, which are helpful in understanding the relevance of signaling pathways in terms of their applica- tion to the prediction of cancer patients’ survival time. These

Compared with the Cox PH model, the deep learning models showed improved pre-
diction accuracy by flexibly integrating a large number of genomics features without strong parametric assumptions,

To identify the potentially associated signaling pathways of hidden nodes, the
Pearson’s correlation values between the expression of individual genes and the output
of the given hidden nodes were calculated to identify the most linearly correlated genes.
Then, gene set enrichment analysis (GSEA) [13] was employed to link the hidden nodes
with the enriched signaling pathways. 

instead of using multi-omics data of a large num- ber of genes, a set of cancer signaling pathways were modeled using a simplified and partially biological meaningful deep neural network architecture

investigate the relevance or influence of these signaling pathways within the context of survival outcome prediction using a biologically meaningful and simplified deep learn- ing model, DeepSigSurvNet.

To interpret deep learning models’ prediction, a set of inter- pretation and explaining approaches have been proposed, e.g., the smmothgrad [18] and Layer-Wise Relevance Propagation (LRP) approach [19], to identify the features that can influence the model prediction results. Interestingly, the interpretable analysis using the smoothgrad approach identified distinct probability density distribution patterns of these signaling pathways, which can be helpful in understanding the relevance of the signaling pathways in terms of their association with cancer patients’ survival. 

There are 303 pathways in the KEGG database, and 45 of them are annotated as “signaling pathways”. 

46 signaling pathways (45 signaling pathways + cell cycle) are selected

In the ‘input layer’, there were two input features, i.e., normalized gene expression across TCGA samples and integer copy number variation, for each gene.

Then, the genes’ state were connected to the 46 signaling pathways only if a gene was included in a signal- ing pathway (not a full connection layer). The gene connection matrix and pathway connection matrix were used to design the connections.

The output of the 46 sign- aling pathways was used as the input for the convolution and inception [21] layers

The inception [21] module used multiple kernel filter sizes in each layer, instead of stacking more layers sequentially. It can capture informative features via the dimension reduction and reduce the vanishing gradient problem. The

To better model and predict the survival time of cancer patients, three clinical factors (age, gender and stage) and the vital status were concatenated with the genomics data. To

To investigate the relevance of individual signaling pathways in survival time prediction, we employed the smooth- grad approach, which is available in the “iNNvestigate” package [22].

we employed the ‘iNNvestigate’ package to calculate the relevance scores of the individual signaling pathways on individual cancer patients in each of the four types of cancer



### Bayesian network model for identification of pathways by integrating protein interaction with genetic interaction data - Fu et al. 2017

probabilistic graphical model to develop a new method that integrates genetic interaction and protein interaction data and infers exquisitely detailed pathway structure. We

the pathway network as Bayesian network and applied this model to infer pathways for the coherent subsets of the global genetic interaction profiles, and the available data set of endoplasmic reticulum genes

a Bayesian model that inte-
grates high-throughput protein and genetic interaction data to reconstruct detailed biological pathway structures. The model can organize related genes into the corre- sponding pathways, arrange the order of genes within each pathway, and decide the orientation of each intercon- nection. Based

Based on protein interaction network, the model predicts detailed pathway structures by using genetic interaction information to delete redundancy edges and reorient the kept edges in the network. 

our model represents a biological pathway network as a Bayesian network [27], in which **each node presents the activity of a gene product**. Different

introducing protein interaction networks as underlying pathway structures (different from activity pathway network APN)

a scoring func- tion is defined by gene pairwise score,

a pathway network as a Bayesian network that is a directed acyclic graph. The

The activity of a gene is assigned to a node in the network [26].

The edge in the network is an interaction in protein interaction network

it presents the conditional dependency be- tween the nodes connected as well.

we only utilize conditional independence assumptions of the Bayesian network theory to construct a network that can represent independence assumptions hidden in the gene interaction data.

we score it in term of gen- etic interaction quantitative measurement

For every pair of genes, there are four topo- logical structures and their local scores shown

Despite the larger score indicating the more possible local structure for each gene pair, we still need every one of four scores to find the optimal global structure. We computed the four possible scores for each pair of genes

we can com- pute **a local score for every pair of genes in a candidate pathway network** N, and sum up all of the scores for all pairs to define the global score function f(N), to which the Bayesian network posterior probability distribution p(N|D) is proportional

Different from study of Ref. [26], [APN] we do not include
every edge score in f(N), because **the edge in our network represents protein interaction that insures its existence**. Then, it avoids the dilemma how to adjust the balance between the two scores.

We utilized annealed importance sampling [26, 28] to learn the pathway structure by the above distribu- tion p

The annealed importance sam- pling approach can assign weights to pathway
networks sampled by simulated annealing schedules, then to evaluate that converge to the real network structure. The

we propose a Bayesian network model to identify pathway structures by integrating protein inter- action with genetic interaction data. Our approach makes use of the complementarity between protein (physical) and genetic (functional) interaction data to refer the biological pathway structures.

#### Identification of active signaling pathways by integrating gene expression and protein interaction data - Kabir et al. 2018

we present a method to simultaneously predict the set of active signaling pathways, together with their pathway structure, by integrating protein-protein interaction network and gene expression data. We

A number of bioinformatics methods have been proposed for the reconstruction of known signaling pathways by using PPI data. For ex- ample, CASCADE_SCAN generates a specific pathway for a list of protein molecules using a steepest descent method. That is, the method takes the input proteins and then finds their interaction partners iteratively based on some evidences (i.e.,

Pathlinker reconstructs the known sig- naling pathways by taking a subnetwork of PPI that con- sists of the Rs and TFs of interest [13].

methods that combine PPI and genetic
interaction data to identify signaling pathway structure. The activity pathway network (APN) approach utilizes high-throughput genetic interaction data and applies the Bayesian learning method to identify detailed structure of known signaling pathways

All of the above methods aim to restructure the top-
ologies of known signaling pathways.

no open-source methods have been reported that simultaneously and comprehensively identify the set of active signaling pathways and the likely pathway structures for a gene expression profile

we propose an ap- proach to systematically identify the set of active recep- tor-mediated signaling pathways within any given cell, by combining PPI and gene expression data. [SPAGI]

we collected the known R [receptors], K [kinases] and TF signaling molecules (2134 genes/proteins in total) from public data sets [

* R: from a curated database of the Fantom5
  project [24].
* K: from the Uniprot curated database
* TF: from a database of sequence-spe- cific DNA-binding TFs identified by gene ontology (GO) based annotation

PPI data from STRING database (version 10) [26] to obtain all currently known PPIs for the 2134 known R/K/TF signaling molecules - while

we have considered here all the physical and other inferred (e.g., co-expression) interactions when defining PPIs to maximize our abil- ity to detect the full network structure.

We selected PPIs defined by STRING as ‘high confidence’ (i.e. confidence_score > = 700)

This thresholding yielded 16,550 and 19,502 PPIs for mouse and human respectively. After obtaining these highly scored PPIs both for the human and mouse organisms we have merged all the PPIs by assuming that the molecules have one-to-one homology mapping between the organisms

From the combined high scored PPIs, we collected only the PPIs for the signaling pathways that have interactions able to make full paths from R to K to TF

the word ‘path’ is defined as a single R/K/TF prediction, whereas the word ‘pathway’ is defined as the collection of paths that all start from the same R (i.e., all paths de- fined by a single R constitutes a pathway

For each potential signaling pathway, we first calculated the proportion of active molecules (defined as highly expressed genes based on the above high expression threshold) for each path. We then summed all the pro- portions of all the paths for the pathway and divided the total proportion value by the total number of paths of the pathway. This final value was termed the Activity score (As)

It should be noted that as currently applied, the
SPAGI method detects receptor-mediated signaling pathways. Modification of the SPAGI approach could be used to identify other cellular control mechanisms involving PPIs independent of TFs. Also, at this stage it is not clear whether the other pathways highly rankedbythe activityscore aretruly active,aspro- tein expression and protein activation state (e.g., via phosphorylation) within a tissue cannot be deter- mined from gene expression data. 



### Pathway-guided deep neural network toward interpretable and predictive modeling of drug sensitivity - Deng et al. 2020

the biological knowledge of pathways, we reshaped the canonical DNN structure by incorporating a layer of pathway nodes and their connections to input gene nodes, which makes the DNN model more interpretable and predictive compared to canonical DNN

we reshape the canon- ical DNN structure by introducing a layer of pathway nodes and their connections to input gene nodes and drug target nodes. The pathway layer is followed by several fully-connected layers and an output layer. We conducted extensive performance evaluations on multiple independent drug sensitivity data sets, and demonstrated that our model significantly out- performed not only canonical DNN model, but also eight other classical regression models.

Our empirical experiments show that our method successfully makes advantage of both the excellent predictive ability of deep neural network and the biological knowledge of pathways, and achieves pharmacological interpretability and predictive ability in modeling drug sensi- tivity in cancer cells

To explore the pharmacological mechanism of action, we take the pathways to construct a DNN model, by introducing a layer of pathway nodes and their con- nections to input gene nodes. The

The pathway layer is followed by two fully-connected hidden layers and an output layer.

each neuron in the input layer represents a gene (cell line feature or drug target), while the corresponding output is the quantitative drug sensitivity.

for a drug tested on a cell line, the input values of cell line nodes are normalized gene expression levels, and the input values of drug target nodes are the STITCH confidence scores if targeted and zero if not targeted by the assayed drug

The first hidden layer, referred to as the pathway layer, has 323 neurons derived from the KEGG pathway repository. The connections between the input layer and the pathway layer are determined by the associations between genes and pathways. A 323*1278 mask matrix, denoted by M, was used to encode the relation between gene nodes and pathway nodes, in which 1 represents the existence of an association between the gene and pathway node, and 0 otherwise. In the back propagation process, the weights of the edges are iteratively updated using the rule as below:
Wmask = Wfeedback ∗M

we compared the outputs of the pathway
nodes with drug treatments to those without drug treatments (corresponding drug target feature nodes set to 0). we observed remarkable decreases in disease-related pathway nodes upon drug treatments. 

the remarkable output decrease of disease-related pathway nodes during forward propagation upon inputs of drug targets can imply the inhibition effect on corre- sponding pathways induced by drug treatment of cancer cells. In

By integration of our current knowledge presented in the form of pathways maps into the
DNN structure, i.e. the pre-specified connections between the input layer and the pathway layer, our model can extract the features in the pathway level, rather than gene level, to pre- dict higher cellular process and final phenotypical changes via subsequent fully-connected layers. 

, we observed re- markable activity changes of disease-related pathway nodes upon drug input signals, which implies significant impact on disease-related pathways induced by drug treatments with high sensitivity on cancer cells. In contrast, most of irrelevant pathway nodes show trivial changes in their activity. 

### Inference of patient-specific pathway activities from multi-dimensional cancer genomics data using PARADIGM - Vaske et al. 2010

s. A gene is modeled by a factor graph as a set of interconnected variables encoding the expression and known activity of a gene and its products, allowing the incorporation of many types of omic data as evidence. 

The method predicts the degree to which a pathway’s activities (e.g. internal gene states, interactions or high- level ‘outputs’) are altered in the patient using probabilistic inference

a PGM framework based on factor graphs (Kschischang et al., 2001) that can integrate any number of genomic and functional genomic datasets to infer the molecular pathways altered in a patient sample.

A gene family was created whenever the cross-reference for a BioPAX protein listed proteins from distinct genes. Gene families represent collections of genes in which any single gene is sufficient to perform a specific function

We also extracted abstract processes, such as ‘apoptosis,’ that refer to general processes that can be found in the NCI collection.

PARADIGM produces a matrix of integrated pathway activities (IPAs) A where Aij represents the inferred activity of entity i in patient sample j. The matrix A can then be used in place of the original constituent datasets to identify associations with clinical outcomes

We first convert each NCI pathway into a distinct probabilistic model

A pathway diagram from NCI was converted into a factor graph that includes both hidden and observed states. The factor graph integrates observations on gene- and biological process-related state information with a structure describing known interactions among the entities

, we use variables
to describe the states of entities in a cell, such as a particular mRNA or complex, and use factors to represent the interactions and information flow between these entities.

These variables represent the differential state of each entity in comparison with a ‘control’ or normal level rather than the direct concentrations of the molecular entities. This representation allows us to model many high-throughput datasets

It also allows for many types of regulatory relationships among genes

**The factor graph encodes the state of a cell using a random variable**
**for each entity X={x1,x2,...,xn} and a set of m non-negative functions, or factors, that constrain the entities to take on biologically meaningful values as functions of one another. **

PARADIGM models various types of interactions across genes including transcription factors to targets (upper-left), subunits aggregating in a complex (upper- right), post-translational modification (lower-left) and sets of genes in a family performing redundant functions (lower-right).

Each entity can take on one of three states corresponding to activated,
nominal or deactivated relative to a control level (e.g. as measured in normal tissue) and encoded as 1, 0 or −1 respectively. The states may be interpreted differently depending on the type of entity (e.g. gene, protein, etc). For example, an activated mRNA entity represents overexpression, while an activated genomic copy entity represents more than two copies that are present in the genome.

In order to simplify the construction of factors, we first convert the
pathway into a directed graph, with each edge in the graph labeled with either positive or negative influence. First, for every protein coding gene G, we add edges with a label ‘positive’ from GDNA to GmRNA, from GmRNA to Gprotein and from Gprotein to Gactive to reflect the expression of the gene from its number of copies to the presence of an activated form of its protein product. Every interaction in the pathway is converted to a single edge in the directed graph. Using

Using this directed graph, we then construct a list of factors to specify
the factor graph. For



### Prediction and interpretation of cancer survival using graph convolution neural networks - Ramirez et al. 2021

a novel graph convolution neural network (GCNN) approach called Surv_GCNN to predict the survival rate for 13 different cancer types using the TCGA dataset. For each cancer type, 6 Surv_GCNN models with graphs generated by correlation analysis, GeneMania database, and correlation + GeneMania were trained with and without clinical data to predict the risk score (RS)

A novel network-based interpretation of Surv_GCNN was also proposed to identify potential gene markers for breast cancer. The signatures learned by the nodes in the hidden layer of Surv_GCNN were identified and were linked to potential gene markers by network modularization.

The Cox-PH method employs both quantitative vari- ables such as molecular expression levels, age, and weight and cate- gorical variables including sex and different treatment methods. Furthermore, the Cox regression model extends survival analysis methods to assess the effect of several risk factors on survival time simultaneously. However, the Cox-PH method tends to suffer in high- dimensional data, and regression cannot learn any complex non-linear functions.

A previous remedy to this downside is the application of support vector machines in survival analysis [13,14]. Other machine learning methods have also been applied for survival analysis to extract significant molecular predictors for early diagnosis and optimal treat- ment outcomes [11,15]

deep learning has been used to predict survival outcomes
with the loss function from the Cox regression model

a 1-layer artificial neural network (ANN) to predict sur- vival outcomes

Hoa et al. also created a multilayer model integrating biological pathways and clinical data to predict the Prognostic Index (PI) of each patient

Convolutional Neural Networks (CNNs) have been used to predict the survival of patients

Recurrent Neural Networks were used for time-series data- sets to predict survival outcomes

The challenge of using the CNN method to analyze gene expressions lies in two facts. 

1. Gene expression data is 1-dimensional while CNN needs 2-dimensional (2D) input. 
   * Different embedding methods have been applied to transform the 1D data into 2D. 

2. Further, the convolutional approach works well on the Euclidean manifold. However, molecular expression data, in essence, represents the outcome of molecular interactions, which are in **non- Euclidian manifold** represented by interaction networks [24].
   *  On the other hand, the one-dimensional gene expression data is easily **mappa- ble onto a graph and graph convolutional neural networks (GCNN) is a promising approach for the non-Euclidean manifold**, suggesting the



survival analysis on 13 TCGA cancer types spanning
5,963 samples with 3 different graphs: 

1. a statistical relationship graph using correlation, 
2. a database-driven graph from GeneMania including gene-gene and protein–protein interactions, and 
3. a merged graph including both correlation and GeneMania to examine which graph can generate the best outcome. 

We also incorporated clinical data into our model to increase prediction accuracy and

the GCNN performs a similar operation to the traditional convolu-
tional neural networks but it learns features from neighboring nodes in a graph

interpreting the model was separated into two parts: finding the significant nodes in the hidden layer contributing to the RS and dissect the relationship between gene expression levels to the hidden layer nodes

























**aracne + VIPER**

* MI based
* needs providing TF list
* 
* VIPER sample-based; can test DE
* do not use PPI data 

wgcna/cemi

* co-expr based
* do not use PPI data (post-hoc only)

causalpath

* starts with the PPI - proteo regulation
* comparison-based: too sparse; correlation-based: too rich

size of the results ??? 

size of the regulons in aracne (how  many genes)



methods that rely on single-patient/sample
mining (e.g., VIPER, the PPI network-based method by Chuang et al., and pathER) rely on dataset scaling to define its single- sample signatures (defined by comparing each gene to the average of its expression in the dataset of interest) making interpretation of any findings from such analyses dataset-specific

what are the limitations of enrichment based analyses

