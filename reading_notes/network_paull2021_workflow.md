### A modular master regulator landscape controls cancer transcriptional identity - Paull et al. 2021

To study MR modularity and genetic drivers in 9,738 samples
in The Cancer Genome Atlas (TCGA) (Weinstein et al., 2013), on a sample-by-sample basis, we developed Multi-omics Master- Regulator Analysis (MOMA). MOMA integrates gene expression and genomic alterations profiles to identify MR proteins and MR modules that represent the key effectors of a tumor’s mutational landscape and are thus responsible for implementing its associ- ated cancer cell identity.



**The MOMA framework**

The MOMA framework is shown in both a simplified (Figures 1A– 1C) and a detailed (Figure S1A–S1E) conceptual workflow. Briefly, gene expression profiles from 20 TCGA cohorts (Table S1) were first transformed to protein activity profiles by using the Virtual Proteomics by Enriched Regulon Analysis (VIPER) al- gorithm (Alvarez et al., 2016) (step 1) (Figure S1B). Candidate MR proteins were then identified by Fisher’s integration of p values from (1) their VIPER-measured activity, (2) functional genetic al- terations in their upstream pathways by Driver-Gene Inference by Genetical-Genomic Information Theory (DIGGIT) analysis (Chen et al., 2014), and (3) additional structure and literature- based evidence supporting direct protein-protein interactions between MRs and proteins harboring genetic alterations, via the Predicting Protein-Protein Interactions (PrePPI) algorithm (Zhang et al., 2012) (steps 2 and 3) (Figure S1C). We used the
vector vector of integrated (Log_10 p)2 (MOMA scores) to weigh each MR’s contribution in a tumor subtype clustering step (step 4) (Figure S1D). Finally, genomic saturation analysis upstream of top candidate MRs identified those most likely to control the subtype transcriptional identity (step 5) (Figure S1D). Finally, this was followed by identification and functional characteriza- tion of MR block sub-modules, termed MR-Blocks (MRBs), recurring across multiple subtypes (step 6) (Figure S1E). See STAR methods for a detailed description of each step. Somatic genomic alterations considered by the analysis
include single nucleotide variants and small indels (SNVs), as well as somatic copy number alterations (SCNAs) from the Broad TCGA Firehose pipeline, and fusion events (FUSs) reported by Pipeline for RNA-Sequencing Data Analysis (PRADA) (Torres- Garcı´a et al., 2014)(STAR methods). Alternative or complemen- tary algorithms can be easily incorporated into MOMA, for instance to integrate the effect of germline variants, epigenetic alterations, or extracellular signals. VIPER has been extensively validated as an accurate method-
ology to measure a protein’s activity, on the basis of the enrichment of its tissue-specific activated and repressed tran- scriptional targets (regulon) in over and under-expressed genes (Alvarez et al., 2016)—i.e., akin to a highly multiplexed gene- reporter assay. To generate accurate regulons for 2,506 regulatory proteins annotated as transcription factors (TFs) and co-transcription factors (co-TFs) in Gene Ontology (GO) (Ash- burner et al., 2000; The Gene Ontology Consortium, 2019), we used Algorithm for the Reconstruction of Accurate Cellular Net- works (ARACNe) (Basso et al., 2005); see STAR methods for ARACNe and VIPER accuracy metrics. For each candidate MR, we first identified candidate upstream modulator proteins by using the Conditional Inference of Network Dynamics (CINDy) algorithm (Giorgi et al., 2014) and then assessed whether the presence of genomic alterations in their encoding genes was associated with differential MR activity (activity Quantitative Trait Locus analysis [aQTL]). These two steps comprise the DIGGIT algorithm, which was highly effective in elucidating key driver mutations missed by prior analyses in GBM (Chen et al., 2014).
Tumor



**Tumor subtype identificatio**

To identify tumor subtypes representing distinct transcrip- tional tumor identities regulated by the sameMRproteins, weper- formed partitioning around medioids clustering (PAM) (Park and Jun, 2009), based on protein activity profile similarity between samples, after each protein was weighted by its cohort-specific, integrated MOMA Score (STAR

Proteins with more functional mutations in their upstream pathways were deemed more likely determinants of tumor subtype identity and provided greater weight to the clustering solution. Within each cohort, the optimal number of clusterswas determined by using a cluster reli- ability score (CRS) (Figure 2A; STARmethods). Using

MR-based clustering outperformed expression-based clustering in all 20 cohorts (p

MOMA identified 112 subtypes, representing the stratification of cancer into transcriptional identities regulated
by

MOMA identified subtypes and differ-
ential outcome in cohorts that had been previously challenging from a gene- expression analysis perspective.

MOMA identified transcriptional clusters presenting statistically significant outcome differences in 19 out of 20 co- horts (Figures

Despite its unsupervised nature, MR-based clustering recapit-
ulated established molecular subtypes and outcome differences

**Tumor checkpoint MRs**

A tumor checkpoint is defined as a module with the minimumMR repertoire necessary to implement a tumor’s transcriptional identity by canalizing genomic events in its upstream pathways. We thus used saturation analysis to refine the initial ranked-list of subtype-specific proteins produced by MOMA analysis to a

small set of candidate MRs that optimally account for the sub- type’s genetic landscape (STAR methods). By ‘‘accounting for an alteration’’ we mean that it is either harbored by the MR or by proteins upstream of the MR

If driver mutations occurred mostly in or upstream of tumor checkpoint MRs, saturation should be achieved rapidly, with only few MRs. In contrast, if mutations were randomly distrib- uted, saturation should be very gradual. To test this hypothesis, we considered all previously described genomic events (SNV, SCNA, and FUS). To	

To avoid over counting, we consolidated same-amplicon SCNAs upstream of MRs into single regional events, and further refined these by selecting genomic events identified by Genomic Identification of Significant Targets in Cancer (GISTIC) 2.0. We then plotted the fraction of all such events predicted to be in or upstream of the top n candidate MRs, on a sample by sample basis—averaged over all samples in the same subtype (Figure 3A)—and defined the tumor check- point as the MR set needed to achieve a predefined saturation threshold in each subtype (STAR methods). Finally, we identified 407 recurrent MRs (Table S2) occurring in mR4 subtypes, a sta- tistical threshold determined by a null hypothesis model (Fig- ure S3A; STAR methods). Of these, 37 were highly recurrent, occurring in m R 15 subtypes (Figure 3B). The H3/H4 histone chaperone ASF1B emerged as the most pleiotropic MR (m = 31 subtypes), followed by MYBL2 (m = 30), JUP (m = 29), TOP2A (m = 25), and TRIP13 (m = 25)

Consistent with the tumor checkpoint hypothesis, we observed
rapid genomic event saturation in all but 3 subtypes (subtypes S1, S3,and S4 of ovarian cancer). For

At the saturation point, ?50% of all genomic events were ac- counted for,

This supports the role of tumor checkpoints as regulatory bottlenecks responsible for canalizing upstream mutations and suggests that < 50% of all genomic events might be actual passengers

To further assess MOMA’s ability to differentiate between driver and passenger events, we assessed the differential
enrichment enrichment of mutations upstream of MRs in either GISTIC2.0 and CHASM-predicted driver events or all genomic events re- ported by the TCGA Firehose pipeline (

averaged across all MOMA-inferred subtypes of a specific TCGA cancer cohort, differential enrichment of the former was highly statistically significant across

Similar mutational co-segregation differences were detected across virtually all cohort subtypes

Regional (i.e., non-focal) SCNAs have been largely ignored by previous analyses because of their high gene content. However, the DIGGIT analysis module in MOMA is effective at removing regional SCNA genes that are unlikely to modulate MR activity. When DIGGIT-refined regional SCNAs were included, subtypes became highly homogeneous in terms of their mutational reper- toire across patients. Consider,

**Tumor checkpoints are hyperconnected and modular**

Analysis of existing molecular interaction networks confirmed that tumor checkpoints represent hyperconnected modules, compared with equisized protein sets chosen at random from 2,506 regulatory proteins, as a null model. Networks include Hu- manNet 2.0 (Hwang et al., 2019) (p < 5.0 3 10?42, by Kolmo- gorov-Smirnov) (Figure S4A), Multinet (Khurana et al., 2013) (p < 2.0 3 10?37)(Figure S4B), and PrePPI (Zhang et al., 2012) (p = 9.0 3 10?44)(Figure 

We then tested whether subtype-specific tumor checkpoints
might be decomposed into finer-grain MR sub-modules— recurrent across multiple subtypes—representing pancancer core-regulatory structures. Clustering of 407 MRs identified by saturation and recurrence analysis yielded 24 MRBs as an optimal solution (Figure S5A), with each MR assigned only to a single MRB (core-set). Given that individual TFs might perform different functions, depending on interacting co-partners (e.g., MYC/MAX versus MYC/MIZ-1), we used a ‘‘fuzzy’’ clustering al- gorithm to further refine core-sets with additional non-unique MRs (Miyamoto et al., 2008)(Figures S5B and S5C; Table S4; STAR methods)

Each tumor checkpoint is thus deconstructed into a specific combination of activated or inactivated MRBs (Figure 5A), where the MRB activity is computed as the average activity of all of its MRs. Transcriptional targets of individual MRB MRs were statis- tically significantly enriched in Cancer Hallmarks (Drake et al., 2016; Liberzon et al., 2015) and KEGG/Reactome categories

Enrichment of Tumor Hallmarks, KEGG, and Reac- tome categories in genes altered upstream of each MRB was generic and sparser (Table S4), suggesting that functional spec- ificity is manifested after MRB integration, rather than in the up- stream genetics.

**Tumor checkpoint MRs are enriched in essential proteins**

We further assessed whether inferred tumor checkpoint MRs were enriched in essential proteins, based on Achilles Project data (Cowley et al., 2014); see Figure S5E for a conceptual work- flow

cell lines optimally matching MOMA-inferred subtypes were identified by protein activity analysis (STAR methods). Essentiality was then assessed on the basis of Achil- les’ score in matched cell lines.

MRs were highly enriched in essential genes

We then tested MRB-specific essentiality. enriched for cell viability hallmarks ; no essential MRS were found in immune-related MRBs (MRB:10, MRB:19, MRB:22, MRB:23, and MRB:24)—consistent with lack of immune function in cell lines. However,

 the role of many of these MRs in pancancer inflammation was previously reported (Thorsson et al., 2018). This suggests that MOMA can identify MRs that are relevant in a human tumor context but might be missed in viability assays in vitro.

**MRBs improve outcome analysis**

To assess whether MRBs could stratify patient outcome, we used a sparse Lasso COX proportional hazards regression model (Tibshirani, 1997), with MRB activities as predictors. Of the 20 TCGA cohorts, 16 could be effectively stratified, often with highly improved p values compared with those of tumor checkpoint stratification (Figures

To assess whether TCGA-inferred MRBs generalize to other cohorts, we analyzed the METABRIC breast cancer cohort,
Figure

**MRB:2 canalizes driver mutations in prostate cancer** 

To validate the effect of genetic alterations affecting MRB activ- ity, we selected MRB:2, the most recurrently activated across all subtypes (40 out of 112) (Figure 5A). By regularized COX regres- sion, MRB:2 produced some of the largest outcome regression coefficients across TCGA

FOXM1 and CENPF—the 6th and 13th most recurrent MRs (Figure 3B)— were validated as synergistic MRs

). We ranked MOMA-inferred alterations upstream of MRB:2 on the basis of their statistical sig- nificance across all TCGA cohorts and selected those with the strongest MRB:2 association (Figures 6E and 6F; STAR methods), most of which were not identified as drivers by Mut- Sig2.CV (Lawrence et al., 2013) and mutation assessor (Reva

d 6 loss-of-function MRB:2-associated events for experimental validation, including ; ideally suited to detecting activity in- crease in loss-of-function assays. Two short-hairpin RNA (shRNA) hairpins per target were used. Functional; Functional and tumori-genic effects were assessed both in vitro and in vivo (

VIPER analysis following shRNA-mediated silencing of four of
the five candidate genes versus negative controls, revealed sta- tistically significant activity increase of MRB:2 activity, based on its8core-set MRs(

increase in cell migration, as assessed by wound healing assays at the indi- cated time points in relation to control cells infected with scramble shRNAs (Figures 7C, 7D, S7A) This was confirmed by Boyden chamber migration assays (Figures

22rv1 cells engrafted in immune deficient mice

**Pharmacological MRB modulation** 

We then asked whether MRB activity and associated function might be pharmacologically modulated. We

MOMA analysis recapitulated these roles in terms of hallmark enrichments, including androgen and estrogen response, EMT, apical surface and apical junction, and inflammatory response

LNCaP cells treated with the AR antagonist enzalutamide or dimethyl sulfoxide (DMSO) (Handle et al., 2019) confirmed that MRB:14 genes have AR- dependent expression (Figure

MRB:14 activity effectively stratified luminal versus basal samples in BRCA and BLCA TCGA cohorts, by PAM50 classification (Figure S7F), further supporting MRB:14’s role as a positive determinant of hor- mone-signal-mediated luminal state across tissues and loss of luminal identity when inactivated

VIPER analysis of patient-matched biopsies pre- and post- androgen deprivation therapy (ADT) (Rajan et al., 2014) showed pronounced MRB:14 MR activity suppression (Figure S7G). Indeed, metastatic, post-ADT tumors are generally basal-like having undergone EMT, raising the question of whether pro- longed ADT might induce loss of adhesion and metastatic progression (Sun et al., 2012; Tsai et al., 2018). Intermittent testosterone replacement therapy reduced appearance of aggressive tumors (Chuu et al., 2011; Loeb et al., 2017), reflect- ing potential benefit of periodic, AR-mediated cell adhesion reinforcement

To test whether pharmacological activation of MRB:14 MRs
might reduce the migratory, EMT-related potential of aggressive prostate cancer, we used the OncoTreat algorithm (Alvarez et al., 2018) to prioritize 120 FDA-approved and 217 late-stage (phase II and III) experimental drugs, on the basis of their overall ability to activate MRB:14 MRs, by using RNA Sequencing (RNA-seq) profiles of AR-resistant DU145 cells at 24 h after treatment (STAR methods). Four MRB:14-activating drugs were inferred at physiologically realistic concentrations (<10 mM), including fe- dratinib, pevonedistat, ENMD-2076, and lexibulin (Figure 7G), and their effect was assessed in wound healing assays. All four drugs but none of the negative controls significantly inhibited DU145 cell migration at 24 h (Figures 7H and 7I). The latter—triapine, raltiterxed, and dorsomorphin—were randomly selected among drugs with no significant MRB:14 activity effect (Figure





### A modular master regulator landscape controls cancer transcriptional identity - Paull et al. 2021 - *workflow*

**MOMA framework**

STEP 1: VIPER: gene expression -> protein activity profiles

STEP2 and 3

Candidate MR proteins: Fisher’s integration of p values from (log-transf. -> MOMA scores)

(1) their VIPER-measured activity (using ARACNe regulons)

(2) functional genetic al- terations in their upstream pathways by Driver-Gene Inference by Genetical-Genomic Information Theory (DIGGIT) analysis 

 (3) additional structure and literature- based evidence supporting direct protein-protein interactions between MRs and proteins harboring genetic alterations, via the Predicting Protein-Protein Interactions (PrePPI) algorithm 

STEP 4 

MOMA scores to weigh each MR’s contribution in a tumor subtype clustering step 

STEP5

genomic saturation analysis upstream of top candidate MRs identified those most likely to control the subtype transcriptional identity

STEP 6

identification and functional characteriza- tion of MR block sub-modules, termed MR-Blocks (MRBs), recurring across multiple subtypes (step 6)

--

Somatic genomic alterations considered by the analysis: SNVs and SCNAs from the Broad TCGA Firehose pipeline, and fusion events (FUSs) reported by PRADA

ARACNe to generate accurate regulons for 2,506 regulatory proteins annotated as transcription factors (TFs) and co-transcription factors (co-TFs) in Gene Ontology (GO)

For each candidate MR:

1) identifiy candidate upstream modulator proteins by using the Conditional Inference of Network Dynamics (CINDy) algorithm
2) assessed whether the presence of genomic alterations in their encoding genes was associated with differential MR activity (activity Quantitative Trait Locus analysis [aQTL]).

(DIGGIT algorithm)

**Tumor subtype identification**

partitioning around medioids clustering (PAM) based on protein activity profile similarity between samples, after each protein was weighted by its cohort-specific, integrated MOMA Score

Proteins with more functional mutations in their upstream pathways were deemed more likely determinants of tumor subtype identity and provided greater weight to the clustering solution.

a cluster reli- ability score (CRS) to determine the optimal number of clusters

multiple statistically equivalent solutions were identified, the one yielding the best survival stratification was selected

comparison with gene expression based subtypes analysis; comparison to molecular subtypes and outcome stratification

**Tumor checkpoints** (saturation analysis)

saturation analysis to refine the initial ranked-list of subtype-specific proteins produced by MOMA analysis to a small set of candidate MRs that optimally account for the sub- type’s genetic landscape (STAR methods) [ ‘‘accounting for an alteration’’ = either harbored by the MR or by proteins upstream of the MR]

If driver mutations occurred mostly in or upstream of tumor checkpoint MRs, saturation should be achieved rapidly, with only few MRs.

If mutations were randomly distrib- uted, saturation should be very gradual

To test this hypothesis

- considered all previously described genomic events (SNV, SCNA, and FUS)
- avoid over counting: consolidate same-amplicon SCNAs upstream of MRs into single regional events; further refinement by selecting genomic events identified by GISTIC 
-  plott the fraction of all such events predicted to be in or upstream of the top n candidate MRs, on a sample by sample basis
- tumor check- point defined as the MR set needed to achieve a predefined saturation threshold in each subtype 
- sta- tistical threshold determined by a null hypothesis model

assess ability to differentiate driver vs. passenger events: enrichment of mutations upstream of MRs in either GISTIC2.0 and CHASM-predicted driver events or all genomic events re- ported by the TCGA Firehose pipeline (

Regional (i.e., non-focal) SCNAs have been largely ignored by previous analyses because of their high gene content. However, the DIGGIT analysis module in MOMA is effective at removing regional SCNA genes that are unlikely to modulate MR activity. When DIGGIT-refined regional SCNAs were included, subtypes became highly homogeneous in terms of their mutational reper- toire across patients. Consider,

**Tumor checkpoints are hyperconnected and modular** 

Analysis of existing molecular interaction networks: comparison  with equisized protein sets chosen at random from  2,506 regulatory proteins, as a null model

Networks from: Hu- manNet 2.0, Multinet ,PPI

decomposition of subtype-specific tumor checkpoints into finer-grain MR sub-modules recurrent across multiple subtypes—representing pancancer core-regulatory structures: MRs clustering; number of clusters identified by saturation and recurrence analysis; ‘fuzzy’’ clustering al- gorithm to further refine core-sets with additional non-unique MRs (Miyamoto

deconstruct each tumor checkpoint into a specific combination of activated or inactivated MRBs; MRB activity computed as the average activity of all of its MRs

Transcriptional targets of individual MRB MRs were statis- tically significantly enriched in Cancer Hallmarks (Drake et al., 2016; Liberzon et al., 2015) and KEGG/Reactome categories

Enrichment of Tumor Hallmarks, KEGG, and Reac- tome categories in genes altered upstream of each MRB

**Tumor checkpoint MRs are enriched in essential proteins**

enrichment in essential proteins, based on Achilles Project data:

- protein activity analysis for optimal match between cell lines and MOMA-inferred subtypes 
- essentially assessed based on Achilles' score

test MRB-specific essentiality

**MRBs improve outcome analysis**

assess patient outcome stratification with sparse Lasso COX proportional hazards regression model (MRB activities as predictors)

generalize TCGA-inferred MRBs to other dataset (METABRIC)

**MRB:2 canalizes driver mutations in prostate cancer** 

validation of MOMA results in PRAC

experimental validation in vitro and in vivo (shRNA)

VIPER analysis following shRNA-mediated silencing

cell migration assays (wound healing assays and Boyden chamber migration assays)

mice engrafting

**Pharmacological MRB modulation**

assess whether MRB activity and associated function might be pharmacologically modulated

focus on a specific MRB (link to hormonal therapy benefit)

functional annotation of MRB:14 related genes; hormone dependent expression validate in vitro; subtype stratification; VIPER analysis of patient-matched pre- and post-ADT biopsies

test whether pharmacological activation of MRB:14 MRs might reduce the migratory, EMT-related potential of aggressive prostate cancer: OncoTreat algorithm to prioritize drugs, on the basis of their overall ability to activate MRB:14 MRs, by using RNA Sequencing (RNA-seq) profiles of AR-resistant DU145 cells at 24 h after treatment



### A modular master regulator landscape controls cancer transcriptional identity - Paull et al. 2021 - *Methods*

**Activity inference**

Transcriptome-wide expression signatures were computed by two non-parametric transformations. First, each column (tumor sample) was rank transformed and scaled between 0 and 1. Then each row (gene) was rank transformed and scaled between 0 and 1. Finally,

regulatory protein activity was measured by the VIPER algorithm (Alvarez et al., 2016), using tis- sue-matched ARACNE regulons 

**Enrichment analysis**

aREA analysis The analytic Rank-based Enrichment Analysis (aREA) was introduced in (Alvarez et al., 2016) as an analytical methodology to assess gene set enrichment analysis statistics, producing results that are virtually identical to GSEA (Subramanian et al., 2005) without the need for time-consuming sample or gene shufflin

**DIGGIT analysis**

improved version of the original algorithm.

original algorithm combines:

(a) MINDy analysis: identify proteins representing candidate upstream modulators of a MR protein

(b) aQTL analysis: to identify genomic events in candidate upstream modulators associated with statistically significant differential MR activity

(c) conditional association analysis: eliminate genomic events no longer significant given another genomic event

improvement: 

(a) rather than using mutual information, aQTL statistical significance is assessed by aREA-based enrichment analysis of samples, ranked by differential activity of the specific MR, in samples harboring a specific SNV or SCNA events, 

(b) MINDy algorithm replaced by CINDy: more accurate implementation of the con- ditional mutual information foundation of the algorithm,

(c) conditional association analysis eliminated because it produced too many statistical ties when applied to pancancer cohorts

aQTL analysis was performed only for events occurring in >= 4 samples since fewer events are highly unlikely to achieve statistical significance 



**CINDy score** 

3 steps

1) rank proteins by their VIPER statistical significance, integrate across all cohort samples using the Stouffer’s method for p value integration
2) for each statistically significant differentially active protein (= candidate MR), compute the conditional mutual information between the expression of the MRand of its regulon genes, given the expression of any gene harboring a somatic event -> identify mutation-harboring genes encoding for proteins that affect the ability of aMR to regulate its targets
3) for each event type (i.e., SNV, amplified SCNA, or deleted SCNA) all statistically significant CINDy scores for a given MRwere integrated using Stouffer’s method to produce three global CINDy scores; fusion events not analyzed since ARACNe not designed to identify targets of fusion proteins (only aQTL analysis step applied for FUSs)
   

**aQTL score**

3 steps

1.  same as CINDy (could be further improved by integrating across individual subtypes rather than entire cohorts)
2. for each candidate MR, statistical significance of the aQTL event assessed by computing the enrichment of all cohort samples, ranked by the MR’s dif- ferential activity, in samples harboring the event, using aREA
3. for each event type, a global aQTL score (SaQTL) computed representing the integration of all statistically significant MR-event aQTL p-values per MR for that event type, using Stouffer’s method; producing 3 scores (for SNVs, small indels, fusion events). Because < 100 statistically significant CINDy mod- ulators indicates that the dataset is too small for a properly powered CINDy analysis:
   * If >= 100 CINDy-inferred MR modulators were identified in a given cohort, only aQTLs for somatic events harbored by genes with a statistically significant CINDy p-value were integrated. 
   * Otherwise, the p-values of all statistically significant aQTLs were integrated independent of CINDy results. 

**PrePPI score**

PrePPI to identify structure-based protein-protein interactions between proteins encoded by genes harboring a somatic event and each MR protein. 

3 steps

1. same as CINDy and aQTL
2. assign high-confidence interactions in the PrePPI database an empirical
   p-value as follows: 
   1. rank them by likelihood scores
   2. compute p-values as the fraction of interactions with equal or better rank, normalized by the total number of PrePPI interactions in the database
3. for each event type,global PrePPI score (SP) by integrating individual p-values of all statistically significant PrePPI interactions for that event type, using Fisher’s method. This produced three global PrePPI scores.

**Integrated rankings and MOMA scores** 

2 steps

1. For each candidate MR, integrate p-values corresponding to same-type events as assessed by aQTL, PrePPI, and CINDy, by integrating using Stouffer’s method; fusion events: CINDy and PrePPI scores not computed (thus not integrated); aQTL analysis: fusion events considered equivalent to SNVs. This produced 9 integrated p-values for each statistically significant, candidate MR protein
2. rank all proteins in a cohort based on VIPER score, then Stouffer’s method to integrate the 9 p-values for each statistically significant protein (i.e., candidate MR) with its VIPER p-value. This yields a global MOMA p-value which represents the probability that a protein may be a bona fide MR by chance; MOMA score from transformation of MOMA p-value

**Cluster reliability score (CRS)**

a way to assess the fit of each sample within a cluster ((Alvarez et al., 2018)

ample distance matrix computed by taking the weighted VIPER scores for each sample (VIPER activity values multiplied by each MR’s MOMA Score) and calculating the pairwise Pearson correla- tions. 

Activity-based clustering 

clustering of each tissue-specific VIPER activity matrix using k-medoids clustering (distance matrix:  weighted Pearson correlation between VIPER-inferred protein activity vectors; weights: square of the integrated MOMA scores -> increase contribution of high-scoring MRs 

CRS calculated for each sample for each k-value; optimal number of clusters: first local maximum for the Global CRS. 

**Expression-based clustering Similar to Protein Activity-based clustering**

each tissue-specific gene expression matrix was clustered using k-medoids clustering with k set as the same value chosen for the tissue-specific VIPER activity clustering

distance matrix: Pearson correlation between gene expression profiles

**Survival analysis** 

fit a Cox proportional hazards model to each sample grouping defined by the initial clustering

‘‘best’’ survival clusters = the one with the lowest proportion of observed to expected death events

‘‘worst’’ survival = the highest observed/expected ratio

fit a second Cox model exclusively to samples from those two clusters and calculated the significance of survival differences between ‘‘best’’ and ‘‘worst’’ clusters in that model

**Saturation analysis** 

Saturation curves generated by ascertaining the number of functional somatic events upstream of the N most statistically sig-nificant candidate MR proteins, ranked by their global MOMA score

define saturation threshold: 

- assess how many functional somatic events upstream of the first half of all regulatory proteins in that subtype, thus conservatively excluding proteins with a non-statistically significant VIPER activity
- saturation threshold set at 85% of that number

assess how many of the N proteins with the highest VIPER activity were needed to identify that number of somatic events in their upstream pathways

saturation increased so rapidly and significantly, compared to an identical number of randomly selected regulatory proteins (null hypothesis), that increases in event number for N > 100 MRs were not statistically significant. To avoid contaminating functional genomic events with passenger ones, by using non-significant MRs to assess saturation: select a more conservative saturation threshold 

**Driver mutation enrichment**

sample-specific analysis in each cohort to assess the statistical significance of somatic event enrichment, upstream of checkpoint MRs:

- identify activated MRs and their upstream somatic events as for Saturation analysis
- for each sample, compute ratio of all validated CHASM and GISTIC2.0 putative driver events vs  total number of events 
- assess cohort-level significance by comparing the number of samples with a ratio > 1 against a one-tailed binomial null distribution 

**MRB analysis**

*Tumor checkpoints comprise multiple sub- modular structures, termed MR-Blocks (MRBs), which regulate specific tumor hallmarks and are recurrently detected across different subtypes*

MRs identify by saturation analysis that were also statistically significant in >= subtypes (recurrence analysis) clus- tered based on their VIPER-inferred activity, using a Euclidean distance metric and partitioning around medoids (PAM)

to compute the Euclidean distance, each MR was associated with a 112-dimensional vector representing its VIPER-inferred activity in each subtype. 

Cluster Fitness score defined as the average CRS for all MRs in a cluste

cluster identified by this analysis then expanded by  *m* MRs with the best average Euclidean distance to those in the core-set: for each additional MRs in each MRB, the trace of the covariance matrix of the Tumor Hallmark enrichment across the 24 MRBs was calculated to assess the total variance of the solution . variance showed optimal increase for *m* =6



**Achilles essentiality** 

Achilles shRNA DEMETER knockout scores from The Broad Institute for all cell lines in CCLE for all TFs and co-TFs analyzed by MOMA

set appropriate null hypothesis on a gene by gene basis

for each of the MOMA subtypes, match the MR activity vector, weighted by the cohort-specific MOMA score of each MR, to the protein activity profile of each CCLE cell line, using an algorithm included in VIPER (identify the cell lines that best recapitulates subtype-specific MRs as possible dependencies)

assess the essentiality of each MR in cell lines that were significant matches vs those providing clear non-matches (non-parametric rank-based Mann-Whitney-Wilcox test)

significant FDRs after multiple hypothesis correction considered as as essential subtype-specific MRs using null model; stratification of essentiality for each MR across the subtypes where that MR was statistically significantly active. To calculate statistical significance of the enrichment of essential genes, a null model was built by

**MRB:12 validation**

lentiviral-mediated gene silencing

perturbation dataset VIPER analysis

wound healing assays

matrigel invasion assays

xenograft assays

**MRB:14 validation**

analysis of LNCaP cells

drug prioritization: dataset of protein activity profiles of drug response, as inferred from a screening of durgs in prostate cancer cell line

aREA function from VIPER compute a Normalized Enrichment Score (NES) for each drug; 
based on the enrichment of differentially activated proteins, converted to p-values and corrected for multiple hypothesis testing, used as a score to prioritize drugs and statistically significant drugs

Analysis of prostate cells and tumor biopsies: metaVIPER approach from VIPER to generate interactomes; regulons pruned to the top 100 targets with the highest likelihood using the pruneRegulon from VIPER; enrichment analysis on VIPER-inferred protein activity signatures and use resultant NES scores

