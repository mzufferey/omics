### Aibar et al. 2017

scenic, a computational method for simultaneous gene regulatory network reconstruction and cell-state identification from single-cell rna-seq data (http://scenic. aertslab.org). on a compendium of single-cell data from tumors and brain, we demonstrate that cis-regulatory analysis can be exploited to guide the identification of transcription factors and cell states. scenic provides critical biological insights into the mechanisms driving cellular heterogeneity.



SCENIC workflow consists of three steps (Fig. 1a, Supplementary Fig. 1 and see Online Methods). In the first step, sets of genes that are coexpressed with TFs are identified using GENIE3 (ref. 8) (Supplementary Fig. 1a). Since the GENIE3 modules are only based on coexpression, they may include many false positives and indirect targets. To identify putative direct-binding targets, each coexpression module is subjected to cis-regulatory motif analysis using RcisTarget (Supplementary Fig. 1b and see Online Methods). Only modules with significant motif enrichment of the correct upstream regulator are retained, and they are pruned to remove indirect targets lacking motif support. We refer to these processed modules as regulons.

we developed the AUCell algorithm to score
the activity of each regulon in each cell (Supplementary Figs. 1c and 2, and see Online Methods). For a given regulon, comparing AUCell scores across cells makes it possible to identify which cells have significantly higher subnetwork activity. The resulting binary activity matrix has reduced dimensionality, which can be useful for downstream analyses.

SCENIC workflow. SCENIC is a workflow based on three new R/bioconductor packages: (i) GENIE3, to identify potential TF targets based on coexpression; (ii) RcisTarget, to perform the TF- motif enrichment analysis and identify the direct targets (regu- lons); and (iii) AUCell, to score the activity of regulons (or other gene sets) on single cells. We

GENIE3. GENIE3 (ref. 8) is a method for inferring gene regula- tory networks from gene expression data. In brief, it trains ran- dom forest models predicting the expression of each gene in the data set and uses as input the expression of the TFs. The different models are then used to derive weights for the TFs, measuring their respective relevance for the prediction of the expression of each target gene. The highest weights can be translated into TF-target regulatory links8. Since GENIE3 uses random-forest regression, it has the added value of allowing complex (e.g., non- linear) coexpression relationships between a TF and its candidate targets

The input to GENIE3 is an expression matrix. The preferred
expression values are gene-summarized counts (which might or might not use unique molecular identifiers, UMIs25). Other meas- urements, such as counts or transcripts per million (TPM) and FPKM/RPKM are also accepted as input. However, note that the first network-inference step is based on coexpression, and some authors recommend avoiding within-sample normalizations (i.e., TPM) for this task because they may induce artificial covaria- tion26

The output of GENIE3 is a table with the genes, the potential regulators, and their ‘importance measure’ (IM), which represents the weight that the TF (input gene) has in the prediction of the target. We

GRNBoost. GRNBoost is based on the same concept as GENIE3: inferring regulators for each target gene purely from the gene expression matrix. However, GRNBoost does so using the gra- dient-boosting machines (GBM)28 implementation from the XGBoost library29. A GBM is an ensemble learning algorithm that uses boosting30 as a strategy to combine multiple weak learners, like shallow trees, into a strong one. This contrasts with random forest, the method used by GENIE3, which uses bag- ging (bootstrap aggregation) for model averaging to improve regression accuracy. GRNBoost uses gradient-boosted stumps (regression trees of depth 1)31 as the base learner. GRNBoost’s main contribution is casting this multiple regression approach into a Map/Reduce32 framework based on Apache Spark2

RcisTarget. RcisTarget is a new R/Bioconductor implementation of the motif enrichment framework of i-cisTarget and iRegulon RcisTarget identifies enriched TF-binding motifs and candidate transcription factors for a gene list. In brief, RcisTarget is based on two steps. First, it selects DNA motifs that are significantly over- represented in the surroundings of the transcription start site (TSS) of the genes in the gene set. This is achieved by applying a recovery-based method on a database that contains genome-wide cross-species rankings for each motif. The motifs that are anno- tated to the corresponding TF and obtain a normalized enrich- ment score (NES) > 3.0 are retained. Next, for each motif and gene set, RcisTarget predicts candidate target genes (i.e., genes in the gene set that are ranked above the leading edge). 

AUCell. AUCell is a new method that allows researchers to iden- tify cells with active gene regulatory networks in single-cell RNA- seq data. The input to AUCell is a gene set, and the output is the gene set ‘activity’ in each cell. In SCENIC, these gene sets are the regulons, which consist of the TFs and their putative targets. AUCell calculates the enrichment of the regulon as an area under the recovery curve (AUC) across the ranking of all genes in a par- ticular cell, whereby genes are ranked by their expression value. This method is therefore independent of the gene expression units and the normalization procedure

Cell clustering based on gene regulatory networks. The cell regulon activity is summarized in a matrix in which the columns represent the cells and the rows the regulons. In the binary regu- lon activity matrix, the coordinates of the matrix that correspond to active regulons in a given cell will contain a “1,” and “0” oth- erwise. The equivalent matrix, which contains the continuous AUC values for each cell regulon, is normally referred to as the AUC activity matrix. Clustering of either of the regulon activ- ity matrices reveals groups of regulons (jointly, a network) that are recurrently active across a subset of cells. The binary activity matrix tends to highlight higher order similarities across cells (and therefore highly reduces batch effects and technical biases); on the other hand, the AUC matrix allows researchers to observe more subtle changes. 



### Tovar et al. 2015

880 microar-
ray expression profiles from several experimental datasets that are available on the Gene Expression Omnibus site (http://www. ncbi.nlm.nih.gov/geo/GEO)



First, we generated a network for every known human TF in the primary breast cancer gene expression dataset by using the Algorithm for the Reconstruction of Accurate Cellu- lar Networks (ARACNE) (Basso et al., 2005; Margolin et al., 2006). ARACNE is a computational algorithm widely used to identify statistical relationships among genes, by calculating the mutual information (MI) between gene pairs from microarray expression data (Basso et al., 2005; Margolin et al., 2004)



MARINais designed to infer transcription factors that control the
transition between two phenotypes A and B, as well as the mainte- by the activation or repression of specific TFs then their targets nance of the latter phenotype. If the A→B transition is supported should be among the most differentially expressed genes between the two cellular phenotypes, with activated and repressed targets at opposite ends of the expression range. MARINa estimates the importance and biological relevance of a TF on a given phenotype by computing the statistical significance of the overlap between its regulon and the gene expression signature using sample per- mutation to estimate the distribution of the enrichment score (ES) (Subramanian et al., 2005) in the null condition (Lefebvre et al. 2010)



An important step for this algorithm is the selection of the Transcription Factors, since they will determine the rest of the cal- culation. A proper annotation of transcription factors is crucial for an accurate description of the process under investigation. Here, we used the HGU133A annotation file, in which we found 1142 TFs (Supplementary Material 1). This list was compared with other three lists. Those lists are available in Shimoni and Alvarez (2013), Vaquerizas et al. (2009) and http://www.bioguo.org/AnimalTFDB/ Animal Transcription Factor DataBase, respectively. We want to stress that all four lists show consistency among them



shadowingand synergy ofTMRsover their target gene sets (Lefebvre et al., 2010). 



Causal network (CN) analysis

(the Ingenuity Knowledge Base
(IKB)). IKB reports a series of experimentally observed cause- effect relationships related to transcription, expression, activation, molecular modification, binding events and transport processes



Network informa- tion was supplied with differential expression analysis between tumoral and healthy samples, which defined the input for the CN study. Differentially



Causal network (CN) analysis was performed with the Inge-
nuity Pathway Analysis method ((IPA®,



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



### Paull et al. 2021 - MOMA

, using a network-based approach, identified 407 master regulator (MR) proteins responsible for canalizing the genetics of individual samples from 20 cohorts in The Cancer Genome Atlas (TCGA) into 112 transcriptionally distinct tumor subtypes

. MR proteins could be further organized into 24 pan-cancer, master regulator block modules (MRBs), each regulating key cancer hallmarks and predictive of patient outcome in multiple cohorts. Of

Of all somatic alterations detected in each individual sample, >50%were predicted to induce aberrant MR activity, yielding insight into mechanisms linking tumor genetics and transcriptional identity and establishing non-oncogene dependencies.

step1: gene expression profiles from 20 TCGA cohorts (Table S1) were first transformed to protein activity profiles by using the Virtual Proteomics by Enriched Regulon Analysis (VIPER) al- gorithm

Candidate MR proteins were then identified by Fisher’s integration of p values from (1) their VIPER-measured activity, (2) functional genetic al- terations in their upstream pathways by Driver-Gene Inference by Genetical-Genomic Information Theory (DIGGIT) analysis (Chen et al., 2014), and (3) additional structure and literature- based evidence supporting direct protein-protein interactions between MRs and proteins harboring genetic alterations, via the Predicting Protein-Protein Interactions (PrePPI) algorithm (Zhang et al., 2012) (steps 2 and 3) (Figure



vector of integrated (Log_10 p)2 (MOMA scores) to weigh each MR’s contribution in a tumor subtype clustering step (step 4) (Figure S1D). Finally, genomic saturation analysis upstream of top candidate MRs identified those most likely to control the subtype transcriptional identity (step 5) (Figure S1D). Finally, this was followed by identification and functional characteriza- tion of MR block sub-modules, termed MR-Blocks (MRBs), recurring across multiple subtypes (step 6) (Figure S1E). See STAR methods for a detailed description of each step

VIPER has been extensively validated as an accurate method-
ology to measure a protein’s activity, on the basis of the enrichment of its tissue-specific activated and repressed tran- scriptional targets (regulon) in over and under-expressed genes (Alvarez et al., 2016)—i.e., akin to a highly multiplexed gene- reporter assay. To generate accurate regulons for 2,506 regulatory proteins annotated as transcription factors (TFs) and co-transcription factors (co-TFs) in Gene Ontology (GO) (Ash- burner et al., 2000; The Gene Ontology Consortium, 2019), we used Algorithm for the Reconstruction of Accurate Cellular Net- works (ARACNe) (Basso et al., 2005); see STAR methods for ARACNe and VIPER accuracy metrics

For each candidate MR, we first identified candidate upstream modulator proteins by using the Conditional Inference of Network Dynamics (CINDy) algorithm (Giorgi et al., 2014) and then assessed whether the presence of genomic alterations in their encoding genes was associated with differential MR activity (activity Quantitative Trait Locus analysis [aQTL]). These two steps comprise the DIGGIT algorithm, which was highly effective in elucidating key driver mutations missed by prior analyses in GBM

. Minimum cohort size reflected the need to generate accurate regulatory network models by using the ARACNe algorithm



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



##### Madhamshettiwar et al. 2012

comparative evaluation of nine state-of-the art gene regulatory network inference methods encompassing the main algorithmic approaches (mutual information, correlation, partial correlation, random forests, support vector machines) using 38 simulated datasets and empirical serous papillary ovarian adenocarcinoma expression-microarray data. We



### Yuan et al. 2019

we propose an approach, which use biweight midcorrelation to measure the correlation between factors and make use of nonconvex penalty based sparse regression for gene regulatory network inference (BMNPGRN). BMNCGRN incorporates multi-omics data (including DNA methylation and copy number variation) and their interactions in gene regulatory network model. The experimental results on synthetic datasets show that BMNPGRN outperforms popular and state-of-the-art methods (including DCGRN, ARACNE and CLR) under false positive control. Furthermore, we applied BMNPGRN on breast cancer (BRCA) data from The Cancer Genome Atlas database and provided gene regulatory network.



we propose an approach, which use
Biweight Midcorrelation to measure the correlation between factors and make use of Nonconvex Penalty based sparse regression for Gene Regulatory Network inference (BMNPGRN). BMNPGRN integrates heterogeneous multi-omics data to infer gene regulatory network with gene expression data, DNA methylation data and copy number variation data. In order to infer gene regulatory network. Firstly, we combine biweight midcorrelation coefficient algorithm, which is an efficient algorithm for computing correlation coefficient, with ‘differential correlation strategy’ to learn associations among DNA methylation sites. Then, nonconvex penalty based sparse regression is used to find gene-related biological factors, and the parameter of method is determined by cross-validation. Meanwhile, nonconvex penalty based sparse regression is used under stability selection which can control false positives effectively [29]. Finally, BMNPGRN identifies gene regulatory network based on the probabilities of biological factors (i.e. gene, DNA methylation sites, and CNV)



Our proposed approach BMNPGRN has advantages
over existing gene regulatory network inference methods. Firstly, BMNPGRN can find more DNA methylation sites which are associated with gene regulatory network. Such method provide deeper insight into gene regulation mechanism. Secondly, BMNPGRN can effectively control false positives using stability selection strategy. Furthermore, BMNPGRN can more accurately find biological factors using nonconvex penalty based sparse regression. Finally, to the best of our knowledge, BMNPGRN is the first method which is applied to breast cancer data obtained by high-throughput sequencing technology



##### Brunner et al. 2021

a novel stochastic dynamical systems model that predicts gene expression levels from methylation data of genes in a given GRN. 

https://github.com/kordk/stoch_epi_lib.

we developed a dynamic interaction network model [25] that depends on
epigenetic changes in a gene regulatory network (GRN). Dynamical systems integrate a set of simple interactions (i.e., transcription factor (TF) binding to a promoter region and subsequent gene expression) across time to produce a temporal simulation of a physical process (i.e., gene regulation in a given GRN). Therefore, the predictions of a dynamical systems model (e.g., TF binding and unbinding events, gene expression levels) emerge from a mechanistic understanding of a process rather than the associations between data (e.g., predicting an outcome from a set of predictor variables). A dynamical systems model can predict gene expression using epigenetic data and a GRN by simulating hypothesized mechanisms of transcriptional regulation. Such models provide predictions based directly on these biological hypotheses, and provide easy to interpret mechanistic explanations for their predictions. The dynamical systems approach offers a number of unique characteristics. First, a stochastic dynamical system provides us with a distribution of gene expression estimates, representing the possibilities that may occur within the cell. Next, the mechanistic nature of the approach means that the model can provide a biological explanation of its predictions in the form of a predicted activity level of various gene-gene regulatory interactions. Finally, a dynamical systems approach allows for the prediction of the effects of a change to the network. To our knowledge, there are no studies that have taken a dynamical systems approach to predicting gene expression from methylation data and a GRN.





### Boltz et al. 2019

determined collective influencers (CI), defined as network nodes that damage the integrity of the underlying networks to the utmost degree

Morone and Makse22 introduced an optimization method to determine an optimized set of nodes termed collective influencers (CI). Such nodes were obtained via optimal percolation theory through the investigation of their propensity to damage the underlying network, strongly emphasizing the role of weakly con- nected nodes. Wondering if collective influencers in protein-protein interaction network carry biological signif- icance, we expected that collective protein influencers were enriched with e.g. disease or essential genes. Within the human interactome, CI proteins were indeed enriched with essential, regulatory, signaling and disease genes as well as drug targets, strongly suggesting that such well-defined protein groups have significance. Furthermore, we found that CI proteins were evolutionarily conserved as CI proteins in evolutionarily conserved networks in different organisms.
Results

, a set of collective influencers (CI) is defined as the minimum set of nodes that, upon deletion, destroy the largest connected component of the underlying network22

. Calculating a score for each protein that reflects its propensity to damage the underlying largest connected component, proteins with the largest score were removed in each step. The procedure stopped when the largest connected component disappeared, providing a list of removed nodes as collective influencers (CI). W



The determination of collective influencers is based on optimal percolation, aiming at the determination of a minimum set of nodes that fragments the underlying network.



The collective influence theory for optimal percolation is based on the message passing equations of the per-
colation process. For

ref22 = Morone, F. & Makse, H. A. Influence maximization in complex networks through optimal percolation. Nature



##### Azad et al. 2021

methodologies with a fully Bayesian approach in discovering novel driver bio-markers in aberrant STPs given high-throughput gene expression (GE) data.

‘PathTurbEr’ (Pathway Perturbation Driver) uses the GE dataset derived from the lapatinib (an EGFR/HER dual inhibitor) sensitive and resistant samples from breast cancer cell lines (SKBR3).

Differential expression analysis revealed 512 differentially expressed genes (DEGs) and their pathway enrichment revealed 13 highly perturbed singalling pathways in lapatinib resistance, including PI3K-AKT, Chemokine, Hippo and TGF-β singalling pathways.

the aberration in TGF-β STP was modelled as a causal Bayesian network (BN) using three MCMC sampling methods, i.e. Neighbourhood sampler (NS) and Hit-and-Run (HAR) sampler that potentially yield robust inference with lower chances of getting stuck at local optima and faster convergence compared to other state-of-art methods. Next,

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



a fully Bayesian statistical modelling approach for analysing the statistical aspect of the perturbed STP structure yielded from the MCMC sampling algorithms of BN structure learning (see previous subsection). We have used p1-model, initially proposed by Holland and Leinhardt [18], is a special class of exponential families of distributions, for this study that offers robust and flexible parametric models, which are used to evaluate the probability that a gene to be hub in the perturbation network inferred in previous subsection

, the exponential family of distribution (i.e. p1-model) is a common choice for Pr(u) as it explicitly allows the distributions parameters to be tied with dif- ferent network statistics (e.g. global density of nodes, in-degree, out-degree, etc.) that control the formation of that network at the first place

we have employed a fully Bayesian approach to infer the posteriori ofparameters of the above p1-model. For that purpose, MCMC sampling methods i.e. Gibbs sampling approach were adopted, which

. Our approach relies on a hierarchical Bayesian model 



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



##### Tripathi et al. 2017

sgnesR (Stochastic Gene Network Expression Simulator in R) is an R package that provides an interface to simulate gene expression data from a given gene network using the stochastic simulation algorithm (SSA)



##### Zhang et al. 2020

we propose a new framework with a new metric to identify driver modules with low-frequency mutation genes, called iCDModule. Inspired by the gravity model, we integrate the coverage and mutual exclusivity in mutation information, define a new metric between gene pairs, called mutation impact distance, to help identifying potential driver genes sets, including those have extremely low mutation rates but play an important role in functional networks. A genetic network is constructed by combining the defined mutation impact distance and then the driver module identification problem is formalized as the maximum clique solution problem, and an improved ant colony optimization algorithm is used to solve it. iCDModule is applied to TCGA breast cancer, glioblastoma, ovarian cancer to test performance.

MDPFinder [9], Multi-dendrix [12], ComMDP and SpeMDP [13] used integer linear programming for settling the problem of maximum cover-exclusive sub-matrix to detect driver pathways.

Hotnet [16], Hotnet2 [17], Hierarchical Hotnet [18], use a thermal diffusion method on PPI network and the diffusion value is used to extract pathways or modules with high con- nectivity. MEMo [19] use interaction networks and functional relationship graphs to derive the largest clique in similar graphs, and combine the mutual exclusivity to process the largest clique. Babur et al. [20] introduce a method based on seed growth on the genetic network. It applies TCGA datasets to detect pan-cancer pathways or modules, and determines the growth strategy through a properly defined by mutual exclusivity score. BeWith [21] put the interaction density and mutual exclusivity as the optimization goal, using an ILP method to solve it. MEMCover [5] combines mutation data and interaction data to identify mutual exclusivity mutation gene in the same or different cancer tissues

we propose a new approach to de novo identify driver modules with low-frequency mutation genes, called iCDModule. The

The mutation score of each gene is calculated based on the mutation damage rate and coverage characteristics to measure the contribution of the mutation to cancer. 

Inspired by the gravity model [22], based on the mutation scoring function and mutation exclusivity, a new metric is defined between gene pairs, namely mutation impact distance. Construct a network of mutation genes with nodes corresponding to the mutation genes. If the mutation impact distance between genes is greater than the average of all non- zero mutation impact distances and they have edge in the PPI network, generate edge between them. The significance of the edges between gene pairs takes into account coverage, exclusivity, and functional correlation.

The problem of driver module identification is transformed into a maximum clique solution problem, and an improved ant colony optimization algorithm is used to settle this problem. Finally,



### Zhang et al. 2011

a linear model to explain similarities between disease phenotypes using gene proximities that are quantified by diffusion kernels of one or more PPI networks. We solve this model via a Bayesian approach, and we derive an analytic form for Bayes factor that naturally measures the strength of association between a query disease and a candidate gene and thus can be used as a score to prioritize candidate genes

This method is intrinsically capable of integrating multiple PPI networks

The Bayesian regression approach can achieve much higher performance than the existing CIPHER approach and the ordinary linear regression method.

Most existing computational methods for inferring causative genes from candidates are formulated as a one-class novelty learning problem that is usually solved with the guilt-by-association principle,  which suggests to compute a score from functional genomics data to quantify the strength of association between a query dis- ease and a candidate gene, and then rank candidate genes according to their scores to facilitate the selection of susceptibility genes [4].

Pearson’s correlation coeffi-
cient of similarities between phenotypes and closeness of genes in a single protein-protein interaction (PPI) network can be used as a concordance score to facilitate the priori- tization of candidate genes

However, PPI networks are far from complete.

we propose a Bayesian
regression approach that can be used with either a single PPI network or multiple networks to prioritize candidate genes. We adopt a linear model to explain disease simi- larity using gene proximity, and we solve this model via a Bayesian approach, which yields an analytic form of Bayes factor for measuring the strength of association between a query disease and a candidate gene. We then use Bayes factors as scores to prioritize candidate genes. We show the validity of assumptions of this approach, and we demonstrate the effectiveness of this approach on five PPI networks via large scale leave-one-out cross- validation experiments and comprehensive statistical analysis. We further show the capability of our approach in integrating multiple PPI networks

known associations between disease phenotypes and genes extracted from the Online Mendelian Inheritance in Man (OMIM) database

he Human Protein Reference Database (HPRD) contains human protein-protein interactions that are manually extracted from the literature by expert biolo- gists [28].

Biological General Repository for Interaction Data- sets (BioGRID) contains protein and genetic interactions of major model organism specie

the Biomolecular Interaction Network Database (BIND) contains both high-throughput and manually curated interactions between biological molecules [30]. 

the IntAct molecular interaction database (IntAct) contains protein-protein interaction derived from literature [

he Molecular INTeraction database (MINT) contains information about physical interactions between pro- teins [32]. From

We adopt a linear regression model to explain disease similarities in the phenotype similarity profile using gene similarities in one or more gene proximity profiles, and we solve this regression model via a Bayesian approach [34].



##### Gosline et al. 2012

SAMNet, for Simultaneous Analysis of Multiple Networks, that is able to interpret diverse assays over multiple perturbations. The algorithm uses a constrained optimization approach to integrate mRNA expression data with upstream genes, selecting edges in the protein–protein interaction network that best explain the changes across all perturbations

The result is a putative set of protein interactions that succinctly summarizes the results from all experiments, highlighting the network elements unique to each perturbation.

the human dataset measured cellular changes in four different lung cancer models of Epithelial-Mesenchymal Transition (EMT), a crucial process in tumor metastasis. SAMNet

SAMNet, for Simultaneous Analysis of Multiple Networks, an algorithm that uses a network flow model to integrate two distinct high-throughput experiments across multiple conditions.

SAMNet uses a constrained optimization formulation based
on the multiple commodity flow problem to model multiple experiments simultaneously as ‘‘commodities’’ that must transit from a common source to a common sink through a shared protein interaction network. Each edge in the interaction network has a particular capacity, and therefore must be ‘shared’ by all commodities. This constraint forces the algorithm to select interactions that are unique to each cellular perturba- tion, thus avoiding the selection of common stress pathways, a common pitfall of other optimization approaches.

SAMNet is a powerful tool for modeling diverse sources of high throughput data across multiple experiments



### Gao et al. 2017

identifying cancer pathway modules through coordination between coverage and exclusivity.

CovEx, to predict the specific patient oriented modules by 1) discovering candidate modules for each considered gene, 2) extracting significant candidates by harmonizing coverage and exclusivity and, 3) further selecting the patient oriented modules based on a set cover model.

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



##### Gabr et al. 2013

, we consider the problem of finding causal orderings of nodes in such protein interaction networks to discover signaling pathways. We adopt color coding technique to address this problem. Color coding method may fail with some probability. 

. Our key contribution in this paper is elimination of the key conservative assumptions made by the traditional color coding methods while computing its success probability. We do this by carefully establishing the relationship between node colors, network topology and success probability. 



### Fu et al. 2017

d probabilistic graphical model to develop a new method that integrates genetic interaction and protein interaction data and infers exquisitely detailed pathway structure. We modeled the pathway network as Bayesian network and applied this model to infer pathways for the coherent subsets of the global genetic interaction profiles, and the available data set of endoplasmic reticulum genes.

The protein interaction data were derived from the BioGRID database. Our

. In order to automatically identify detailed pathway structures using high-throughput genetic interaction data, the activity pathway network (APN) was developed [26]. However, these available approaches cannot fully take advantages of the comple- mentarity between protein and genetic interaction data to infer the biological pathway structures.

we present a Bayesian model that inte-
grates high-throughput protein and genetic interaction data to reconstruct detailed biological pathway structures. The model can organize related genes into the corre- sponding pathways, arrange the order of genes within each pathway, and decide the orientation of each intercon- nection. Based

Based on protein interaction network, the model predicts detailed pathway structures by using genetic interaction information to delete redundancy edges and reorient the kept edges in the network.

Similar to APN [26], our model represents a biological pathway network as a Bayesian network [27], in which each node presents the activity of a gene product. Different from APN that drew network sample from complete network, our method introducing protein interaction networks as underlying pathway structures. 

a scoring func- tion is defined by gene pairwise score, which can avoid the unadjusted balance between gene pairwise score and edge score in the APN. T

our model is able to improve computational efficiency of stochastic simulation algo- rithm and overcome the limitation of APN that some edges in the results are difficult to interpret. In our model, each edge in the network can capture physical docking, and represent functional dependency



##### Erdogdu et al. 2017

Formulate the induction and control of gene regulatory networks (GRNs) from gene expression data using Partially Observable Markov Decision Processes (POMDPs)

partial observability would be a more natural and realistic method for handling the control of GRNs

We propose a method for the construction of POMDP model of GRN from only raw gene expression data which is original and novel. Then, we introduce a novel approach to decompose/factor the POMDP model into sub-POMDP’s in order to solve it efficiently with the help of divide-and-conquer strategy.

we focus on the GRN control problem. The problem
requires maintaining certain expression level for a single gene or certain expression levels for a group of genes

One of the important aspects of GRN control problem is that
it is not possible to obtain complete state information. 

n this paper, we are proposing to model the GRN control prob-
lem in a more natural and realistic way, mainly as a Partially Observable Markov Decision Process (POMDP). In fact, there are some aspects of the problem that are not fully observable, and unfortunately all the above mentioned models assume full observ- ability and hence simplify the problem. We argue that it is only possible to solve the GRN control problem in a realistic setting if partially observability is properly accounted in the model



### Li et al. 2018

We propose a network-based method to detect cancer specific driver modules (CSDM) in a certain cancer type to other cancer types.

We construct the specific network of a cancer by combining specific coverage and mutual exclusivity in all cancer types, to catch the specificity of the cancer at the pathway level.

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
Algorithm



##### Li et al. 2020

To identify driver modules with rarely mutated genes, we propose a functional similarity index to quantify the functional relationship between rarely mutated genes and other ones in the same module. Then, we develop a method to detect Driver Modules with Rarely mutated Genes (DMRG) by incorporating the functional similarity, coverage and mutual exclusivity

If a rarely mutated gene is in- cluded in a driver module, the functional relationships between this rare gene and other ones in the same module can be high. Thus, we propose a functional similarity index to quantify these relationships with network propagation, which can amplify weak similarities between different genes. This functional similarity can successfully capture functional relationship between driver genes with any mutation frequencies in the same driver module, especially for rarely mutated genes. Then, we develop a method to detect Driver Modules with Rarely mutated Genes (aka, DMRG) by incorporating the functional similarity, coverage and mutual exclusivity. 



### Gysi et al. 2020

a method for the systematic comparison of an unlimited number of networks, with unlimited number of transcripts: Co-expression Differential Net- work Analysis (CoDiNA).

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

the ‘‘differential clustering algorithm’’ for revealing conserved and diverged co-expression patterns. 

Existing approaches for comparative gene expression analyses emphasize mostly conserved co-regulation patterns, rather than differences in expression patterns [8,9,11].

To better capture differential expression patterns, we developed a novel approach, termed the differential clustering algorithm (DCA), for systematically characterizing both similarities and differences in the fine structure of co-regulation patterns

The DCA is applied to a set of orthologous genes that are
present in both organisms. As a first step, the pair-wise correlations between these genes are measured in each organism separately, defining two pair-wise correlation matrices (PCMs) of the same dimension (i.e., the number of orthologous genes) (Figure 2A). Next, the PCM of the primary (‘‘reference’’) organism is clustered, assigning genes into subsets that are co-expressed in this organism, but not necessarily in the second (‘‘target’’) organism. Finally, the genes within each co-expressed subgroup are re-ordered, by clustering according to the PCM of the target organism. This procedure is performed twice, reciprocally, such that each PCM is used once for the primary and once for secondary clustering, yielding two distinct orderings of the genes.

The results of the DCA are presented in terms of the rearranged PCMs. Since these matrices are symmetric and refer to the same set of orthologous genes, they can be combined into a single matrix without losing information. Specifically, we join the two PCMs into one composite matrix such that the lower-left triangle depicts the pair-wise correla- tions in the reference organism, while the upper-right triangle depicts the correlations in the target organism (

An automatic scoring method is then applied to classify clusters into one ofthe four conservation categories: full, partial, split, or no conservation of co-expression

##### Kao et al. 2004

For complex transcriptional networks, more sophisticated tools are required to deconvolute the contribution of each regulator. Here, we demonstrate the utility of network component analysis in determining multiple transcription factor activities based on transcriptome profiles and available connectivity information regarding network connectiv- ity. We

NCA takes advantage of the connectivity information to decompose DNA microarray data to determine both TFA and the control strength (CS) of each regulatory pair

everal linear decompositions of the matrix log [Er] have been used in the study of gene expression array, such as singular value decomposition (2) and independent component analysis (3). Although these decomposition tech- niques have strong statistical foundations, their molecular basis is difficult to pinpoint. The solution obtained by NCA is based not on any hypothesis
of relationship between the TFAs, but on the structure of [CS], namely, the connectivity structure of the network linking TFs and genes. Specifically, such constraints involve, for example, setting to zero the elements CSij when gene i is not regulated by TFAj, but can also include constraints on the polarity of the regulation (induction or repression). We demonstrated (1) that if the underlying transcriptional network satisfies the following properties, such decomposition becomes unique up to some normalization factors: (i) The connectivity matrix [CS] must have full-column rank. (ii) When a node in the regulatory layer is removed along with
all of the output nodes connected to it, the resulting network must be characterized by a connectivity matrix that still has full-column rank.
(iii) The log [TFAr] matrix must have full row rank. In other
words, each regulatory signal cannot be expressed as a linear combination of the other regulatory signals. This criterion requires M ⬎ L as a necessary but not sufficient condition



### Li et al. 2020

To explore the dynamics of the mammalian cellular aging
network, we employ non-linear differential equations (Tyson and Novák, 2010) to describe the dynamics of each genes expression in the network. A

A biological system is naturally subject to intrinsic and extrinsic fluctuations. Therefore, we added an additional fluctuation term in the ODEs to characterize the stochastic behaviors of the mammalian cellular aging process. We use the Langevin dynamic approach to simulate the gene circuit dynamics. From the resulting dynamic trajectories of the gene expressions, we collected the statistics and quantified the underlying potential landscape

, we presented a mathematical model to describe the dynamic features of the mammalian cellular aging process. We built the underlying gene regulatory network by integrating the information from previous experimental studies. The genes and wirings in the gene regulatory network were formed, and the dynamics of gene expression was described by nine non- linear ordinary differential equations. 

we have provided a framework to reveal the
underlying mechanism of fast-aging and slow-aging in mammals based on landscape and flux theory. We predict the key genes and interactions in the fast-aging and slow-aging processes.



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



### Guo et al. 2018

identify the personalized-sample driver genes from the cancer omics data due to the lack of samples for each individual. To circumvent this problem, here we present a novel single-sample controller strategy (SCS) to identify personalized driver mutation profiles from net- work controllability perspective

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

to bridge the personalized driver mutation discovery problem and the structural network controllability prob- lem. Therefore,



##### Ma et al. 2016

several methods for active learning of casual networks has been developed recently14–18. The active learning methods utilize both observational and experimental data to discover causal networks. These methods typically first construct a draft of the causal network, generally represented as an unori- ented or partially oriented graph, from observational data. Then, the methods select a variable for experimen- tation/manipulation to further refine the graph. The experimental data obtained from the targeted experiment is used to update the draft of the causal network. The process of variable selection, experimentation, and causal network update is repeated until some termination criterion is satisfied, e.g. all edges in the causal network are oriented. Since randomized controlled experiments are costly, active learning methods employ various heuris- tics when selecting variables for experimentation in order to minimize the required number of experiments. It

Among the active learning methods examined in this study, ODLP variants achieved the best local pathway
reconstruction quality with low cost on the 5 transcription factors examined. In



##### Maere et al.  2008

the distance measures used in traditional clustering algorithms have difficulties in detecting one of the most prominent features of perturbational data, namely partial correlations between expression profiles. Biclustering methods on the other hand are specifically designed to capture such partial correlations. However, most biclustering algorithms do not provide measures for pair-wise expression correlation between genes, but rely on emergent properties of groups of genes and conditions (modules) in order to identify statistically significant subpatterns in the data. This reliance complicates the elucidation of less modular regions in the underlying transcriptional network.

gene expression profiles are discretized into three categories (upregulated, downregulated, unchanged) based on p-values for differential ex- pression. For each pair of profiles, we then assess the probability that the observed overlap of upregulated and downregulated fields is generated by chance. The resulting correlation p-values are corrected for multiple testing and translated to edges in a coexpression net- work, which is then clustered into (overlapping) expression modules using a graph clustering procedure that identifies densely connected components in the network. Relevant condition sets are then determined for all modules and the modules are screened for enrichment of Gene Ontology categories and transcription factor binding sites. Finally, a regulation pro- gram is learned for each module in an attempt to explain the expression behavior of the module’s genes as a function of the expression of a limited set of regulators (transcription factors and signal transducers). 

ENIGMA, that addresses some of these issues. ENIGMA leverages differential expression analysis results to extract expression modules from perturbational gene expression data. 



Our goal was to build a method that: (i) leverages differential expression analy- sis results to extract co-differential expression networks and expression modules from perturbational gene expres- sion data, (ii) is able to detect significant partial coexpression relationships between genes and overlap between modules, (iii) depends on parameters that can be auto- matically optimized or set on reasonably objective grounds. (iv) produces a realistic amount of modules, and (v) visually integrates the expression modules with other data types such as Gene Ontology (GO) information [28], transcription factor (TF) binding data, protein and genetic interactions, in order to facilitate the biological interpreta- tion of the results. 

ENIGMA takes as input a set of perturbational expression data, externally calculated p-values for differ- ential expression (e.g. using the limma package in Biocon- ductor [29]) and other data types if available. ENIGMA uses a novel combinatorial statistic to assess which pairs
of genes are significantly co-differentially expressed (henceforth abbreviated as coexpressed for the purpose of readability). The resulting coexpression p-values are cor- rected for multiple testing and translated to edges in a coexpression network, which is clustered into expression modules (i.e. groups of significantly co-differentially expressed genes) using a graph-based clustering algorithm inspired on the MCODE algorithm [30]. The clustering procedure depends on two parameters that control the density of individual modules and the overlap between modules. The main reason why we chose a two-tier clus- tering approach (data → coexpression network → cluster- ing) is that it allows simulated annealing-based optimization of the clustering parameters to obtain opti- mal coverage of the coexpression network, in terms of module overlap and redundancy. The

In the post- processing phase, ENIGMA determines relevant condition sets for each module, visualizes their substructure and overlap with other modules, screens the modules for enriched GO categories, suggests potential regulators for the modules based on regulator-module coexpression links and enrichment of TF binding sites, and overlays protein and genetic interaction data.



### Mahajan et al. 2021

We combined these epigenomic datasets and integrated them with the reference human protein interactome using a novel network propagation approach.

a network-based approach to integrate multi-omic epigenetic datasets.

At a conceptual level, tumor mutations and features derived from epigenomic datasets drive large-scale cell behavior in similar ways. Individual mutations or chromatin features may exert a weaker, more localized influence, but groups of mutations or epigenetic marks will exert a stronger, systemic influence by influencing multiple components of the underlying molecular network. Thus, both mutations and epigenetic features represent biological priors that are weak or uninformative when considered in isolation but are strong when analyzed together in the context of the underlying topological network structure. A network-based approach can be used to filter out isolated noisy signals and hone in on biologically significant pathways of genes that mutations or epigenomic changes would target at multiple points.

a novel network-based approach to propagate epigenomic signals
derived from our ATAC-seq and CUT&RUN studies and thus uncover network modules that are likely to drive Tfh cell differentiation. Our approach is based on the premise that biologically relevant epigenetic signatures reflect the need for coordinated regulation of genes that encode interacting proteins. A similar conceptual premise has been used to identify subnetworks enriched for mutations in cancer genomics and network-based GWAS (Cowen et al., 2017; Leiserson et al., 2013, 2015; Reyna et al., 2018; Vandin et al., 2011). The

he underlying technique of network propagation represents a powerful way to combine a wide range of signals taking into account the structure of the underlying network. We implemented network propagation using random walk with restart, which is equivalent to an insulated heat diffusion process ((Cowen et al., 2017; Leiserson et al., 2013, 2015; Reyna et al., 2018; Vandin et al., 2011)). We used this approach to unify features derived from the different epigenomic datasets and discover networks driving Tfh cell differentiation. 

we incorporated higher-order relationships between proteins encoded by these genes using the high-quality (i.e., each edge is validated experimentally using multiple independent assays) reference human protein interactome network (Cusick et al., 2009; Das and Yu, 2012).



Each of the epigenetic signals analyzed exhibited distinct trends across the four stages of Tfh cell differentiation. These differences likely reflect complementary regulatory processes underlying the differentiation of this polarized cell type (Fig. 1d). We expect that these processes all involve regulating the expression of an underlying set of components that drive the transition from naive to GC Tfh, and multi-omic integration strategies will be required to uncover these core genes. We focused on the discovery of trajectories using the epigenomic data. The transcriptomic data was used as an orthogonal dataset to validate the identified trajectories.

Next, for each differentiation pattern, peaks from each dataset were aggregated into gene-centric scores using a proximity-weighted count heuristic (Fig. 2b) based on first principles underlying the organization of promoters and enhancers. Gene-centric scores from each dataset were then combined by weighting the datasets equally, as we did not have any a priori information regarding relative information context across the datasets. We refer to each gene’s score as its “peak proximity score” (PPS), which we can calculate for each possible differentiation pattern and use as a surrogate for gene regulation

We ran network propagation through HotNet2 using the combined PPS’s for each protein’s gene as seed values (Figs. 3a, 3b). HotNet2 helps us incorporate the topological structure of the underlying protein interaction network as an informative prior on how proteins are functionally related to one another

The use of network propagation is motivated by the hypothesis that proteins encoded by relevant genes relevant to Tfh cell differentiation are more likely to interact with one another and make up cohesive network submodules (Fig. 3b). Further, our approach can help filter out genes that may not play a role in Tfh cell differentiation but nevertheless show similar epigenetic patterns to relevant genes

network propagation constitutes the final step in a three-step pipeline: identification of the most likely patterns of Tfh cell differentiation, generation of gene-centric scores using a proximity-weighted count heuristic (i.e. calculation of PPS’s), and network propagation using these scores as seeds.

we adapted the HotNet2 algorithm to identify protein network modules enriched for signals derived from ATAC-seq and CUT&RUN datasets obtained from Tfh cells at four stages of differentiation within a healthy human subject’s tonsil. By framing the search for drivers as a search for identifying protein subnetworks with concentrated epigenomic signals, our approach recovered known regulators of Tfh cell differentiation, and we have implicated several novel candidate genes involved in Tfh cell differentiation. 



network propagation represents a flexible and broadly useful approach for studying biological processes that require the integration of disparate epigenomic datasets. This



Calculation of a gene-centric peak proximity scores
For a biological (EBSeq) pattern of interest, peak proximity scores (PPS’s) for a gene (“gene- centric PPS”) were calculated for a biological pattern of interest by calculating a weighted sum of the number of peak clusters from the pattern in the gene’s vicinity using distance-based weights to reflect the varying contributions of promoter (± 10 kb) [3x weight], proximal enhancer (± 25 kb) [2x weight] and distal enhancers (± 250 kb) [1x weight] to gene regulation. Gene-centric PPS scores were first calculated individually for each epigenomic dataset. These were then combined into a final PPS score using equal weights for the different datasets.

Network propagation using random walk with restart
For each pattern, we used network propagation to integrate the gene-centric scores. We used the random walk with restart algorithm (equivalent to an insulated heat diffusion process) as implemented in HotNet2 (Leiserson et al., 2015) to identify high scoring sub-networks. Relevant HotNet2 code is available at https://github.com/raphael-group/hotnet2. We ran network propagation on the union of the high-quality reference human binary and co-complex protein interactomes (Das and Yu, 2012). 



##### Mall et al. 2017

The ability to detect statistical relevant changes in the interaction patterns induced by the progression of the disease can lead to the discovery of novel relevant signatures. Several procedures have been recently proposed to detect sub-network differences in pairwise labeled weighted networks. Methods: In this paper, we propose an improvement over the state-of-the-art based on the Generalized Hamming Distance adopted for evaluating the topological difference between two networks and estimating its statistical significance. The proposed procedure exploits a more effective model selection criteria to generate p-values for statistical significance and is more efficient in terms of computational time and prediction accuracy than literature methods. Moreover, the structure of the proposed algorithm allows for a faster parallelized implementation.



### Matsubara et al. 2019

merges proteome and transcriptome data 

spectral clustering of the network and map resulting eigenvectors from the 2nd and 3d smallest eigenvalues into 2D space which preserves topological distance and cluster structure

map gene expression profiles onto this mapping->generate a set of image-like representation of protein networks processed by deep deconvolutional layers

PPIs from HINT database



##### Wouters et al. 2019

cancer cells from each sample form a distinct cluster per patient, whereas the corresponding normal host cells from various patients cluster together according to their cell type (Tirosh et al. 2016; Puram et al. 2017; Lambrechts et al. 2018). This observation is somewhat counterintuitive because cells with similar gene expression profiles are known to occur in multiple tumors, for instance cells in specific cell cycle stages. 

gene regulatory network inference using SCENIC has been shown to normalize away part of these tumor- specific differences, resulting in one pan-tumor cluster of cycling cells (Aibar et al. 2017). Nonetheless, the unsupervised discovery of common transcriptional states remains a challenge

SCENIC network inference to the single-cell expression matrix

. A transcription factor with its candidate targets is called a regulon. SCENIC yields a regulon-cell matrix with regulon activities across all single cells, and provides therefore an alternative dimensionality reduction. A UMAP visualization based on the regulon-cell matrix reveals three candidate cell states in an unsupervised manner,



##### Xie et al. 2020

methods are urgently needed which can separate the impact of true regulatory elements from stochastic changes and downstream effects. We propose the differential network flow (DNF) method to identify key regulators of progression in development or disease. Given the network representation of consecutive biological states, DNF quantifies the essentiality of each node by differ- ences in the distribution of network flow, which are capable of capturing comprehensive topological dif- ferences from local to global feature domains. 

s a new approach for quantifying the essentiality of genes across networks of different biological states

The existing differential network analy- sis methods mainly fall into two categories. The first category is focused on capturing linear or nonlinear correlation differences in gene expression between two gene regulatory networks (GRNs). For instance, DDN [10] is the first algorithm to detect topological differences by lasso regression in network inference. DISCERN computes a novel perturbation score to capture how likely a given gene has a distinct set of regulators between different condi- tions, which is shown to be robust to errors in network structure estimation. pDNA [12] incorporates prior information into differ- ential network analysis using non-paranormal graphical models, which relaxes the assumption of normality of omics data to find more cancer-related genes



The existing differential network analy- sis methods mainly fall into two categories. The first category is focused on capturing linear or nonlinear correlation differences in gene expression between two gene regulatory networks (GRNs). For instance, DDN [10] is the first algorithm to detect topological differences by lasso regression in network inference. DISCERN computes a novel perturbation score to capture how likely a given gene has a distinct set of regulators between different condi- tions, which is shown to be robust to errors in network structure estimation. pDNA [12] incorporates prior information into differ- ential network analysis using non-paranormal graphical models, which relaxes the assumption of normality of omics data to find more cancer-related genes

The second category is focused on topological differences
between constructed GRNs. For instance, DEC captures the global differential eigenvector connectivity to prioritize nodes in net- works [7]. DiffRank [13] computes the linear combination of differ- ential connectivity and differential betweenness centrality to order genes. DCloc [14] computes the average proportion of changes of each node’s neighborhood as a significance score by iteratively removing edges with different thresholds. DiffNet [15] evaluates topological differences between two networks based on general- ized hamming distance and its statistical significance. TKDS [8] measures the importance of genes by calculating the graphlet vec- tor distance

Two issues remain which the above methods fail to address.
First, all the above methods utilize networks which assume the existence of edges based on co-expression 

The second issue is more subtle. While techniques such as spec-
tral analysis provide a global perspective on connectivity, these approaches fail to encapsulate the flow of information inherent in all biological networks. Network

we propose the differential network flow (DNF)
method to identify key regulators between two networks under different biological conditions. This algorithm is built upon the idea of network flow and information theory. Rewiring of a GRN can be characterized as a dynamic pattern of network flow [19], such a flow-based model captures multiple (from local to global) features of network structure. Information theory is able to quantify the uncertainty in networks, making networks built upon information-theoretic measurements a more acceptable representation of biological systems at the molecular scale [20]. Therefore, DNF is capable of capturing comprehensive topo- logical differences by quantifying the flow in a network.

DNF is built upon the ideas of network flow and information
theory. The novelty of DNF lies in quantifying node-to-node infor- mation entropy according to the network flow in a gene regulatory network, and in characterizing each node as a distribution of net- work flow, which is equal to the distribution of information entropy. The distribution differences of one gene in different net- works represents its essentiality in the biological process responsi- ble for the network’s evolution. Genes are ordered by the magnitude of this difference to establish a ranking

we
employ a three-step process to integrate both transcriptomics and proteomics datasets in GRN-construction. First, a network skeleton is built by differential expression analysis using the tran- scriptomics dataset. Specifically, the skeleton gene sets are selected based on a given criterion, such asjlog2FoldChangej > u, p-value < v, where u represents the fold change of gene expression and v rep- resents the statistical significance of differential expression. Sec- ond, the known corresponding protein–protein interactions in the STRING database (http://string-db.org) are used to establish the gene-gene network for the selected genes. For example, sup- pose p skeleton genes are selected in the first step, then a network skeleton with p nodes is described by an adjacency matrix Ai;j, such that Ai;j >0, i, j =1, ..., p, if protein i and protein j are functionally associated. Finally, the absolute value of the spearman correlation coefficient (scc) of expression is adopted to estimate the strength of connections between adjacent genes, and edges Ai;jfor which scc(i, j) < 0.1 are discarded (see

we construct a pair of GRNs based on a specific transcriptomics dataset (e.g. cancer and control samples) and a generic proteomics dataset



