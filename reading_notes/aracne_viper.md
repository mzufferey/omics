reverse engineering of such networks from genome-wide expression profile

ARACNe (algorithm for the reconstruction of accurate cellular networks)

ARACNe (algorithm for the reconstruction of accurate cellular net- works), a new approach for the reverse engineering of cellular networks from microarray expression profile

ARACNe first identifies statistically significant gene-gene coregulation by mutual information, an information-theoretic measure of relatedness.

It then eliminates indirect relationships, in which two genes are coregulated through one or more intermediaries, by applying a well-known staple of data  transmission theory, the ‘data processing inequality’ (DPI)

relationships included in the final reconstructed network have a high probability of representing either direct regulatory interactions or interactions mediated by post-transcriptional modifiers that are undetectable from gene-expression profiles

availability of a large set of gene-expression profile data representative of perturbations of the cellular systems, leading to the analysis of a broad range of cellular states and associated gene-expression levels. This is necessary because genetic interactions are best inferred when the genes explore a substantial dynamical range.

an equivalent dynamic richness can be efficiently achieved by assembling a considerable number of naturally occurring and experimentally generated phenotypic variations of a given cell type.

ARACNe reconstructed a network suggestive ofa hierarchical, scale-
free organization, with a power-law relationship between the number of genes and their connectivity1

Mutual information. Mutual information for a pair of discrete random variables, x and y, is defined as I(x,y) ¼ S(x)+ S(y) ? S(x,y), where S(t)is the entropy of an arbitrary variable t. Entropy for a discrete variable is defined as the average of the log probability of its states

DPI. First we define two genes, x and y, as indirectly interacting through a third gene, z, if the conditional mutual information I(x,y|z) is equal to zero. The

This inequality is not symmetric, meaning that there may be situations where the triangle inequality is satisfied but x and z may be directly interacting. As a result, by applying the DPI to discard indirect interactions (i.e.,(x, z) relation- ships for which the inequality is satisfied), we may be discarding some direct interactions as well. These are of two kinds: (i) cyclic or acyclic loops with exactly three genes and (ii) sets of three genes whose information exchange is not completely captured by the pairwise marginals

We introduce a percent tolerance for the DPI to account for inaccurate estimates of the difference between two close mutual information values

This has the advantage of avoiding rejection of some borderline edges, resulting in some loops of size three to occur in the predicted topology



Alvarez et al. 2016

Identifying the multiple dysregulated oncoproteins that contribute to tumorigenesis accurate inference of aberrant protein activity in biological samples is still challenging as genetic alterations are only partially predictive and direct measurements of protein activity are generally not feasible.

virtual inference of protein activity by enriched regulon analysis (VIPER), for accurate assessment of protein activity from gene expression data.

We used VIPER to evaluate the functional relevance of genetic alterations in regulatory proteins across all samples in The Cancer Genome Atlas (TCGA).

n addition to accurately infer aberrant protein activity induced by established mutations, we also identified a fraction of tumors with aberrant activity of druggable oncoproteins despite a lack of mutations, and vice versa.

We propose that the expression of the transcriptional targets of a
protein, collectively referred to as its regulon, represent an optimal multiplexed reporter of its activity

no experimentally validated methods to accurately assess the activity of arbitrary proteins in individual samples based on the expression of their regulon genes

regulon analysis, using the master
regulator inference algorithm (MARINa), can help identify aberrantly activated tumor drivers17–21. However, this requires multiple samples representing the same tumor phenotype and cannot be used to assess aberrant protein activity from individual samples.

we introduce a new regulatory-network based approach to infer protein activity from single gene expression profiles (VIPER

VIPER infers protein activity by systematically analyzing expression of the protein’s regulon, which is strongly tumor-context-dependent2

We used the algorithm for the reconstruction of accurate cellular networks (ARACNe23; Online Methods) to systematically infer regulons from tissue-specific gene expression data (

VIPER is based on a probabilistic framework that directly
integrates target ‘mode of regulation’, that is, whether targets are activated or repressed (Fig. 1b and Supplementary Figs. 1 and 2), statistical confidence in regulator-target interactions (Fig. 1b) and target overlap between different regulators (pleiotropy) (Fig. 1d) to compute the enrichment of a protein’s regulon in differentially expressed genes (Online Methods

VIPER uses a fully probabilistic yet efficient enrichment analysis framework, supporting seamless integration of genes with different likelihoods of representing activated, repressed or undetermined targets, and probabilistic weighting of low vs. high-likelihood protein targets. To achieve this, we introduce analytic rank-based enrichment analysis (aREA) a statistical analysis based on the mean of ranks (Fig. 1c and Online Methods). Differential protein activity is quantitatively inferred as the normalized enrichment score computed by aREA.
Systematic

The regulatory networks were reverse engineered by ARACNe4

**The VIPER algorithm tests for regulon enrichment on gene expres- sion signatures.**

1. The <u>gene expression signature</u> is first obtained by comparing two groups of samples representing distinctive phenotypes or treatments

Any method that generates a quantitative measurement of difference between the groups can be used

<u>single-sample-based gene expression signatures</u> can be obtained by comparing the expression levels of each feature in each sample against a set of reference samples by any suitable method, including , including for example Student’s t-test, Z-score transformation or fold change; or relative to the average expres- sion level across all samples when clear reference samples are not available.

2. Then we compute the **enrichment of each regulon on the gene expression signature** using different implementations of aREA
3. Finally, we estimate the significance, including P value and normalized enrichment score

,by comparing each regulon enrichment score to a null model generated by randomly and uniformly permuting the samples 1,000 times. 

when the number of samples is not enough to support permutation with repo- sition (at least five samples per group is required), permutation of the genes in the gene expression signature or its analytic approximation can be used

Analytic rank-based enrichment analysis. aREA tests for a global shift in the positions of each regulon genes when projected on the rank-sorted gene expression signature

we used the mean of the quantile-transformed rank positions as test statistic (enrichment score). The enrichment score is computed twice: first by a one-tail approach, based on the absolute value of the gene expression signature (i.e., genes are rank-sorted from the less invariant between groups to the most differen- tially expressed, regardless of the direction of change); and then by a two-tail approach, where the positions of the genes whose expression is repressed by the regulator (R−) are inverted in the gene expression signature before comput- ing the enrichment score

The one-tail and two-tail enrichment score estimates are integrated while weighting their contribution based on the estimated mode of regulation through a procedure we call ‘three-tail’ approach (see

and then by a two-tail approach, where the positions of the genes whose expression is repressed by the regulator (R−) are inverted in the gene expression signature before comput- ing the enrichment score

The one-tail and two-tail enrichment score estimates are integrated while weighting their contribution based on the estimated mode of regulation through a procedure we call ‘three-tail’ approach (see

The contribution of each target gene from a given regulon to the enrichment score is also weighted based on the regulator-target gene interaction confidence

Finally, the statistical significance for the enrichment score is estimated by comparison to a null model generated by permuting the samples uniformly at random or by an analytic approach equivalent to shuffle the genes in the signatures uniformly at random

Mode of regulation. The mode of regulation (MoR) is determined based on the Spearman’s correlation coefficient (SCC) between the regulator and the target expression, computed from the data set used to reverse engineer the network.



for complex non-monotonic dependencies (for example, for con- text-specific rewiring59–61), assessing the MoR may not be trivial

we first model the SCC probability density for all regulator-target interactions in the network using a three-Gaussian mixture (Supplementary

g (i) clearly repressed targets (MoR−), (ii) clearly activated targets (MoR+), and (iii) non-monotonically regulated targets for which the MoR cannot be reliably estimated (MoRNM).



Regulator-target confidence. We used the mutual information (MI) between regulator and target gene mRNA levels as inference of regulator-target interac- tion confidence. To

To compute a regulator-target interaction confidence score, we first generated a null set of interactions for each regulator by selecting target genes at random from all the profiled genes while excluding those in the actual regulon (i.e., ARACNe inferred). The number of target genes for the null regulon was chosen to match those in the actual regulon. Then we computed a CDF for the MI in the ARACNe regulons (CDF1) and null regulons (CDF2), and estimated the confidence score for a given regulator-target interaction (interaction confidence or IC) as the ratio: IC = CDF1 / (CDF1 + CDF2).

IC used to weight the contribution of each target gene to the enrichment score



Pleiotropy. Pleiotropic regulation of gene expression (genes regulated by several different transcription factors) can lead to false positive results if a non-active regulator shares a significant proportion of its regulon with a bona fide active regulator (

To account for this effect, we extended the shadow analysis procedure originally described in ref. 17 to take full advantage of the probabilistic framework used by VIPER

, we first generated all possible pairs of regulators AB satisfying two conditions: (i) both A and B regulons are significantly enriched in the gene expression signature (P < 0.05), and (ii) they co-regulate (A ∩ B) at least ten genes. Then

Then we evaluate **whether the regulons in each pair are enriched in the gene expression signature mostly due to the co-regulated genes**. This is per- formed by computing the enrichment of the co-regulated genes (A ∩ B) on a subset of the gene expression signature representing only the genes in A (pA) and in B (pB), where pA and pB represent the estimated P values for the enrich- ment computed by aREA.

Then we compute the pleiotropy differential score as PDE = log10(pB) − log10(pA). If pA < pB, we penalize the co-regulated genes for A by PDE PI / NT, where pleiotropy index (PI) is a constant and NT is the number of test pairs involving the regulon A. Conversely, if pA > pB we penalize the co-regulated genes for B by |PDE|PI / NT



######### SUMMARY



ARACNe => reconstruction of cellular networks => infer regulons (= transcriptional targets of a protein)

VIPER => enrichment of a protein's regulon in differentially expressed genes

* target 'mode of regulation' : targets are activated or repressed, statistical confidence in regulator-target interactions 
* target overlap between different regulators (pleiotropy)





usage in Paull et al. 2021



we developed Multi-omics Master- Regulator Analysis (MOMA). MOMA integrates gene expression and genomic alterations profiles to identify MR proteins and MR modules that represent the key effectors of a tumor’s mutational landscape and are thus responsible for implementing its associ- ated cancer cell identity.

gene expression profiles from 20 TCGA cohorts (Table S1) were first transformed to protein activity profiles by using the Virtual Proteomics by Enriched Regulon Analysis (VIPER) al- gorithm (Alvarez et al., 2016) (step 1) (Figure S1B). Candidate MR proteins were then identified by Fisher’s integration of p values from (1) their VIPER-measured activity, (2) functional genetic al- terations in their upstream pathways by Driver-Gene Inference by Genetical-Genomic Information Theory (DIGGIT) analysis (Chen et al., 2014), and (3) additional structure and literature- based evidence supporting direct protein-protein interactions between MRs and proteins harboring genetic alterations, via the Predicting Protein-Protein Interactions (PrePPI) algorithm (Zhang et al., 2012) (steps 2 and 3) (Figure

We used the vector of integrated (Log_10 p)2 (MOMA scores) to weigh each MR’s contribution in a tumor subtype clustering step (step 4

Finally, genomic saturation analysis upstream of top candidate MRs identified those most likely to control the subtype transcriptional identity (step 5)

Finally, this was followed by identification and functional characteriza- tion of MR block sub-modules, termed MR-Blocks (MRBs), recurring across multiple subtypes (step 6

To generate accurate regulons for 2,506 regulatory proteins annotated as transcription factors (TFs) and co-transcription factors (co-TFs) in Gene Ontology (GO) (Ash- burner et al., 2000; The Gene Ontology Consortium, 2019), we used Algorithm for the Reconstruction of Accurate Cellular Net- works (ARACNe)

For each candidate MR, we first identified candidate upstream modulator proteins by using the Conditional Inference of Network Dynamics (CINDy) algorithm (Giorgi et al., 2014) and then assessed whether the presence of genomic alterations in their encoding genes was associated with differential MR activity (activity Quantitative Trait Locus analysis [aQTL]). These two steps comprise the DIGGIT algorithm, which was highly effective in elucidating key driver mutations missed by prior analyses in GBM (Chen

MOMA was used to analyze 9,738 primary samples, from 20 TCGA tumor cohorts (with nR100 samples) (Table S1). Minimum cohort size reflected the need to generate accurate regulatory network models by using the ARACNe algorithm (Basso

Cluster reliability score (CRS) The CRS was introduced in (Alvarez et al., 2018) as a statistically sound way to assess the fit of each sample within a cluster. For each sample, a distance vector V1, representing its distance from all other samples in the same cluster and a vector V2, representing its distance from all other samples in the cohort are computed. The sample distance matrix was computed by taking the weighted VIPER scores for each sample (VIPER activity values multiplied by each MR’s MOMA Score) and calculating the pairwise Pearson correla- tions. The normalized enrichment score of V2 distances, ranked from the largest to the smallest one, in V1 distances, is then assessed using aREA. This produces a p-value that represents the tightness and separation of the cluster being considered in relation to all other samples. A cluster-wide reliability score for each cluster is assessed as the average cluster reliability (NES) of each sample in the cluster, scaled between 0 and 1. Finally, the reliability of the entire clustering solution (global cluster reliability score) is assessed as the average of the cluster-wide reliability score of all clusters in the solution.

MRB analysis The 407 MRs identified by saturation analysis that were also statistically significant in R4 subtypes (recurrence analysis) were clus- tered based on their VIPER-inferred activity, using a Euclidean distance metric and partitioning around medoids (PAM) for k= 2 to 100 clusters (Figure S1E). To compute the Euclidean distance, each MR was associated with a 112-dimensional vector representing its VIPER-inferred activity in each subtype. A Cluster Fitness score was defined as the Average Cluster Reliability Score for all MRs in a cluster. The analysis identified k= 24 as the optimal clustering solution (Figure S5A). Each ‘‘core-set’’ cluster identified by this analysis was then expanded by the m MRs with the best average Euclidean distance to those in the core-set, for m =0, . 100. For each m additional MRs in each MRB, the trace of the covariance matrix of the Tumor Hallmark enrichment across the 24 MRBs was calcu- lated to assess the total variance of the solution. This variance showed optimal increase for m =6(Figure S5B). These optimization steps to ensured uniqueness, specificity, and robustness of the MRB solution.

Gene expression data (counts) from two studies (Rajan et al., 2014; Zhang et al., 2016) were collected. Both studies were analyzed in the same way as follows. Counts downloaded from the GEO portal (GEO: GSE48403, Rajan et al., 2014; and GEO: GSE67070, Zhang et al., 2016) were normalized using the variance stabilizing transformation function available from the DESeq2 package 1.26.0 in R. The metaVIPER approach (Ding et al., 2018), available from the R VIPER package, was then used to generate two interactomes from the TCGAPRAD cohort (this manuscript) and the 2015 SU2C metastatic Castration Resistant Prostate Cancer (mCRPC) cohort (Rob- inson et al., 2015). Regulons were pruned to the top 100 targets with the highest likelihood using the pruneRegulon function of the VIPER package. Gene expression signatures for each individual sample were computed using the method ttest available from the viper function. Enrichment analysis on VIPER-inferred protein activity signatures was computed and resultant NES scores used. Clustering of labeled samples due to similar activation profiles of MRB:14 on patient samples was performed using the hierarchical clustering algorithm available from the ComplexHeatmap package.