---
marp: true
theme: gaia
color: #000
colorSecondary: #333
backgroundColor: #fff
paginate: true
size: 4:3
---
<style>
section {
  font-size: 23px;
}
img[alt~="center"] {
  display: block;
  margin: 0 auto;
}
</style>
<!-- _paginate: false -->
​
# Causal interactions from proteomic profiles: Molecular data meet pathway knowledge


(CausalPath - Babur et al. 2021)



---

### Context

<span style="font-size:20px;">

* classic approach: compile curated pathway models from carefully designed low-throughput controlled experiments
    * well-validated fragments of knowledge, 
    * but extracted from a heterogeneous set of contexts, perturbations, condi- tions, and even organisms: not to making predictions

* data-driven inference: directly infer graphical models, ab initio, from high- throughput measurements 
    * context-specific, predictive models
    * but do not scale in terms of statistical power as the model space is exponentially larger than the observable space


* to alleviate the power issue of the data-driven
approach: get help from prior knowledge when the perturbations in the data are not sufficient to decide between alternative models
* methods that use this strategy, however, use prior knowledge in a reduced form, (e.g. simple interaction networks), omitting mechanistic details and their logical harmony with the new data

</span>

---

### Pathway extraction in perturbation-poor setting

<span style="font-size:20px;">

* the more an experiment lacks extensive perturbations, the more it can benefit from prior knowledge 
* the vast majority of currently available proteomic experiments have either few perturbations (e.g. before/after a stimulation) or only uncontrolled variation (e.g. profiles from disease cohorts), use of prior knowledge in its full potential very important
 $\rightarrow$ **perturbation-poor setting**
* in this setting, model-building activity = selecting parts of the prior knowledge that can best explain the shape of the data = **pathway extraction**

* **CausalPath** is a **pathway extraction method** which uses the rich semantics of curated pathway knowledge, including the type of mechanism, the direction, signs of effect, and post-translational modifications
    * the inferred mechanisms are falsifiable hypotheses that can be experimentally interrogated

</span>

---

### Need for high-confidence priors

<span style="font-size:20px;">

* molecular datasets poorer in perturbations $\rightarrow$ depend more and more on priors for pathway inference, but 
    * priors used in pathway inference methods are generally simple networks (e.g. PPI)
    * they do not inherently suggest any causal relationship between measured molecules
* to extend pathway inference, need for high-confidence priors that contain causality within
* CausalPath for generating causal priors from mechanistic pathway relationships curated from scientific literature, that are then used for the pathway analysis of perturbation-poor proteomic datasets
* CausalPath reduces the pathway inference task into identifying the causal priors that can explain correlations in a given molecular dataset (**pathway extraction**)

</span>

---

### CausalPath 

* experimental data reveal features that change in coordination
* CausalPath automates the search for causal explanations in the literature
    * **maps proteomic profiles to curated human pathways** from multiple resources (integrated in the Pathway Commons database)
    * **detects the potential causal links** in the pathways between measurable molecular features using a **graphical pattern search framework**
    * identifies the subset of the **causal links that can explain correlated changes** in a given set of proteomic and other molecular profiles 

* output: explanations presented as an intuitive network with links to the detailed prior knowledge models and the related literature to create a powerful exploration and analysis platform
* mimics literature search for relationships explaining the data

---

### CausalPath workflow

<u>2 main steps</u> 

1. **detection of causal priors from pathway databases**, performed once and reused in multiple analyses

    * Pathway Commons database
    * BioPAX framework 
        * 12 manually curated *graphical patterns*
        * each pattern captures the control mechanisms over either a phosphorylation of a protein or the expression of a gene

---

Example of a graphical pattern (phosphorylation pattern):


![width:470px center](pictures/example_pattern.png)


---


### CausalPath workflow

<u>2 main steps</u> 

2. **matching causal priors** with supporting correlated changes **in the analyzed data**, performed for every analysis
* *causal conjecture* = pairing of a causal prior with supporting measurements in the molecular dataset that together declare that "one molecular change is the cause of another molecular change"
* 2 forms of molecular changes depending on the experimental setting: 
1. *comparison-based*: up/downregulation for individual features (e.g. "test vs. control")
2. *correlation-based*: positive/negative correlations applying to pairs of features in an uncontrolled study (e.g. uncontrolled cancer cohort)


---

### CausalPath causal conjecture (ternary logical eq.)

<u>Comparison-based analysis</u> 


<span style="font-size:16px;">


$$
\overline{c_{source}\oplus e_{source} \oplus s_{relation} \oplus c_{target}} = true
$$
with $c,e,s \in \{true, false, unknown\}$

* $\oplus$ = ternary operator
  * any operation on $unknown$ will yield $unknown$ result
* $\overline{(...)}$   = logical negation
* $c$ = change of direction of the gene features
  * $true$ in case of upregulation
  * $false$ in case of downregulation
  * $unknown$ if insignificant
* $e$ = the effect of the source feature on its activity
  * $true$ in case of total protein or activating phosphorylation
  * $false$ in case of inactivating phosphorylation
  * $unknown$ in case of a phosphorylation site with unknown effect
* $s$ = the sign of the pathway relation
  * $true$ for phosphorylation and expression upregulation
  * $false$ for dephosphorylation and expression downregulation

</span>

---

### CausalPath causal conjecture (ternary logical eq.)

<u>Correlation-based analysis</u> 

<span style="font-size:16px;">

$$
corr_{source, target}\oplus(e_{source} \oplus s_{relation}) = true
$$

* logical representation of the sign of the correlation replaces the $c$ terms
  * $true$ in case of positive correlation
  * $false$ in case of negative correlation

</span>


<br>
<br>

(NB: correlation-based and comparison-based pseudocodes in Appendix)

---

### CausalPath additional checks


In addition to the logical check by these equations, limit:

1. the phosphoregulations (phosphorylation and dephosphorylation) to the explanation of phosphoprotein changes 
<br>
2. the expressional regulations to the explanation of total protein changes (and optionally mRNA changes)


---

### CausalPath additional statistical measurements

* 2 tests to increase the interpretability of the results
* both relying on data label randomization

1. <u>Network-size test</u>:
    * checks if the correlated changes align with the causal priors in general
    * indicated by **a larger number of interactions in the results than would arise by random chance**
    $\rightarrow$ if significant, additional evidence suggesting that the **data are shaped by the priors**
    <br>
2. <u>Downstream-size test</u>:
    * checks if a protein on the network has **more downstream targets in the results than expected by chance**
    $\rightarrow$ if significant, additional evidence suggesting that **a protein has an influence on the significant number of targets**

---

### CausalPath added value


1. interpretation of complex mechanistic pathway models
<br>
2. site-specific evaluation of phosphoproteomic measurements
<br>
3. a logical test for causality between measurements
<br>
$\rightarrow$ a useful tool for "mechanistic model building"


---

### Limitations

* dependency on observable variance
    * cannot identify a signaling relation that does not significantly vary across the compared samples

* large number of protein phosphorylation sites with unknown functions: downstream cannot be included in the causality network
    * currently, sites with unknown effect and significant change at their signaling downstream reported by CausalPath; users can perform manual curation to increase the coverage of the analysis

* future research
    * work in progress to mine unannotated phosphorylation sites from the literature using NLP tools
    * n-ary systems instead of a binary evaluation of causality (between two proteins) 
    * map a site more accurately to a distinct subset of activities of the protein instead of a binary classification of protein modifications as activating and inhibiting


---

### Conclusion

* context-specific correlations are derived from the data and causality is derived from the literature
* compared with the methods that infer causality from data through mathematical modeling (pathway inference): **wider application area**
    * pathway inference methods have a potential to offer more complete results, but require numerous perturbations and/or time points in the experiments
* CausalPath as a great resource for high-confidence priors
* causality-focused pathway extraction as a mean of transferring knowledge between contexts
    * by detecting variation patterns of proteomic abundances and detecting their consistency with prior knowledge
* importance of using **matching prior-knowledge networks and data**: 
    * gene regulatory networks with transcriptomics 
    * signaling networks with (phospho)proteomics


---

### Conclusion

* CausalPath can be applied to the results of any proteomic and phosphoproteomic experiments to identify differential signaling that is supported by literature knowledge
* requirements:
    * measurement values need to be comparable (normalized) 
    * measurement values need to be associated with related gene symbols
    * phosphopeptide measurements need to specify the phosphorylation sites with respect to their canonical UniProt sequence

$\rightarrow$ CausalPath **generates extremely useful and falsifiable hypotheses** that are unavailable otherwise

$\rightarrow$ since this process systematically considers hundreds of thousands curated mechanisms, it also is **more comprehensive, unbiased, and more consistent** in terms of the generated hypotheses


---
### CausalPath availability

<span style="font-size:20px;">

Data and code availability CausalPath is freely available at http://causalpath.org. Users can upload the proteomic data in a tab-delimited format, along with the analysis parameters, such as how to detect a change in the values. Options include averaging a group of values, getting difference/fold-change of two groups of columns, comparing two groups with a t test, or using correlations in a single group. 

The results are visualized as an interactive network using Cytoscape.js, and the mechanistic details of each interaction can be viewed in SBGN-PD language using a layout algorithm specifically designed for compound graph structures.3

Alternatively, CausalPath can be run locally as a Java application using the sources at https://github.com/PathwayAndDataAnalysis/causalpath. This repository additionally includes examples running CausalPath from R and Python

The generated result networks can be visualized with ChiBE, as well as by uploading analysis output folders to the CausalPath web server at causalpath.org.

</span>

---

---

### Questions/thoughts

* can be used without phosphopeptide measurements ? (available from Karolinska dataset ?)
* outputs of CausalPath difficult to read -> offer a way for easier interpretation / result summary or presentation ? automatize result processing ? (combine with ML methods downstream ? ...) (poor leveraging of the tool in the articles ?!, cf. next slides...)
* use the structure of the network provided by CausalPath as prior in an ANN/for ANN architecture
* add knowledge from scientific articles in addition from databases (e.g. NLP tools to detect association between 2 proteins -> might be more up-to-date than compiled databases) ?

---

### Appendix - CausalPath in Wang et al. 2021

<span style="font-size:16px;">


*Application of CausalPath (Babur et al., 2018) to the protein and phosphoprotein expression data (Figure S7A, Table S5)showed upregulation of the hypoxia pathway in mesenchymal tumors, evi- denced by significant activation of multiple HIF-1 downstream targets (networkpermutationp= 0.0012). Increasedangiogenesis was also evident in mesenchymal tumors, as demonstrated by upregulation of FLT1, MMP14, ENG, and SERPINE1. We observed complex regulation of macrophage activation and po- larization through the upregulation of STAT3, ICAM1, SPI1, and CEBPB. In addition, the M1 polarization marker ARG1 showed increased expression (Arlauckas et al., 2018), along with SERPINE1 and HCK proteins, which promote M2 polarization (Kubala et al., 2018). The elevated inflammatory response in mesenchymal tumors may result in downstream activation of either hypoxia or macrophage polarization through multiple medi- ators, including LANE, IL18, and CD40 (Figure S7A)*

<u>Causative pathway interaction discovery using CausalPath</u> To discover the causative pathway interactions in our proteomic and phosphoproteomic data, we took the normalized expression of protein with < 10% missing values and phosphoprotein with < 25% missing values across all tumor and normal samples as the input to CausalPath (commit 7c5b934). We ran CausalPath in the mode that tests the mean values between test and control groups (value-transformation = significant-change-of-mean), where the test group being the tumors of one subtype and control group being the rest of the tumors. The pathway interaction discovery data source was Pathway Commons v9 (built-in- network-resource-selection = PC). Additionally, we enabled the causal reasoning if all the downstream targets of a gene were active or inactive (calculate-network-significance = true, use-network-significance-for-causal-reasoning = true, permutations- for-significance = 10000). The causative interactions with FDR < 0.05 were extracted and visualized (fdr-threshold-for-data-signif- icance = 0.05 phosphoprotein, fdr-threshold-for-data-significance = 0.05 protein, fdr-threshold-for-network-significance = 0.05). Full result tables were available in Table S5.

</span>

---

### Appendix - CausalPath in Wang et al. 2021

![width:650px center](pictures/wang_example.png)


<span style="font-size:16px;">

*Figure S7. Summary of Pathway Alterations of Mesenchymal Tumors Compared to the Other Tumors. A) Causal explanations for differentially expressed proteomic and phosphoproteomic profiles in mesenchymal tumors versus the rest of the tumors using CausalPath. We highlighted the immune and hypoxia-related interactions enriched in the graph. Genes controlling the immune response might be a potential therapeutic target.*
</span>

---

### Appendix - CausalPath in Satpathy et al. 2021

<span style="font-size:16px;">


*Furthermore, knockdown of these tumor biomarker proteins reduced fitness across 16 LSCC cell lines (https://depmap.org/), suggesting crit- ical roles in key cellular transformation and proliferation pro- cesses (Figure 7C; Table S7).*

<u>CausalPath analysis</u> CausalPath (Babur et al., 2018) searches for known biological mechanisms that can explain correlated proteomic changes in terms of causal hypotheses. We set CausalPath parameters to compare tumors and NATs with a paired t test, used 0.1 as FDR threshold for proteomic change significance and network significance, and detected 5917 potential causal relations between proteins. We repeated the same analysis for each NMF subtype separately and identified 4378 (basal-inclusive), 5334 (classical), 3048 (EMT-en- riched), 3744 (inflamed-secretory), and 4332 (proliferative-primitive) relations. We used these CausalPath network results in the preparation of Figure 7C, identifying potential upstream regulators of oncogenic phosphoproteomic changes. Here an oncogenic phosphoproteomic change can be any of the following 4 events: increase of activating phosphorylation of an oncoprotein, decrease of inactivating phosphorylation of an oncoprotein, decrease of activating phosphorylation of a tumor suppressor protein, and in- crease of inactivating phosphorylation of a tumor suppressor protein. We used the OncoKB database for oncoprotein and tumor suppressor classification (excluded proteins that have both annotations), and used PhosphoSitePlus and Signor databases for the activating/inhibiting classification of phosphorylation sites. In the phosphorylation regulation networks, we included only the targetable regulators (activated proteins) and excluded the untargetable regulators (inactivated proteins).

</span>

---

### Appendix - Supp. Methods


<u>Significance for proteomic data change and correlation</u>

* *comparison-based* analyses:
    * two-tailed t-test for calculating the significance of the difference of the means of the two groups in the comparison, requiring the presence of at least 3 non-missing values from all compared groups
* *correlation-based* analyses:
    * Pearson correlation coefficient and its associated significance, requiring at least 5 samples in the calculation
* both tests assume a null model where molecular readouts change independently. 

The correlation-based causal network provides hypotheses
for the signaling network parts that are differentially active across samples, but it does not indicate which parts are activated together or whether they align with previously defined molecular subtypes. 

---

### Appendix - Supp. Methods

<u>Causality</u>

* pathway inference and pathway extraction are closely related to formal notions of **causality inference**

* a probabilistic causal relationship between 2 events A and B indicates the probability that B depends on the status of A (**Suppes' formulation**)

* this notion can generate predictive models, but does not tell if A may cause B  (e.g. an event X causing both A and B, will satisfy it)

* **Pearl's reformulation to make the model predictive** under an intervention scenario: perturbing the status of A will change the probability of B

* Pearl's notion followed here: detect mechanism-based evidence for **activity change of one protein may affect the abundance of a specific peptide** from another protein in pathway databases

---

### Appendix - Supp. Methods

<u>Derivation of prior relations from detailed mechanistic pathways</u>

* BioPAX-pattern framework
* manual definition of 12 BioPAX patterns to capture potentially causal binary relations that involve phosphorylation and expression of proteins

* causal priors extraction: https://github.com/PathwayAndDataAnalysis/causal-priors-extractor, applied on Pathway Commons 



---

### Appendix - Supp. Methods

Algorithm for selection of explanatory subset of prior relations</u>

Using the extracted causal priors, CausalPath determines if there is sufficient proteomic data that indicates differential activity of that prior



* algorithms that use $.sign$ and $.effect$ properties of the variables
    * take values -1 (false), 0 (unknown) and 1 (true)
    * multiplication of these integer values and checking the result value is an alternative formulation to the original logical equations with the ternary XOR operator


---

### Appendix - Supp. Methods

<u>Pseudocode comparison-based setting</u>


![width:700px center](pictures/comparison_based_pseudocode.png)

when RNAseq data used as the targets of expressional control relations, total protein measurements ($mtt$ in the pseudocode) replaced with the RNAseq measurements of the target ($mrt$)

---

### Appendix - Supp. Methods

<u>Pseudocode correlation-based setting</u>

![width:700px center](pictures/correlation_based_pseudocode.png)

when RNAseq data used as the targets of expressional control relations, total protein measurements ($mtt$ in the pseudocode) replaced with the RNAseq measurements of the target ($mrt$)

---


### Appendix - CausalPath user-tunable parameters

1. Site matching proximity threshold
2. Site effect proximity threshold 
3. Gene focus
4. Generation of a data-centric causal network
5. Protein activity
6. Adjusting phosphopeptide measurements with total protein
7. Data type for expressional targets
8. Using custom resources


---

### Appendix - CausalPath user-tunable parameters

1. <u>Site matching proximity threshold</u>
* default: protein phosphorylation sites in the literature have to exactly match the detected site in the phosphoproteomic dataset to use in causal reasoning 
* might be too strict since there can be slight shifts in the literature, or some nearby sites of proteins are likely to be phosphorylated by the same kinase
* this parameter allows a determined inaccuracy in site mapping 
* increasing it will increase the result network size by allowing proximate site matching
* but new relations in the results are likely to have more false positives 

---

### Appendix - CausalPath user-tunable parameters

2. <u>Site effect proximity threshold</u>
* effect of the phosphorylation sites on the protein activity (as in activating or inactivating) curated by pathway databases, mostly by PhosphoSitePlus
* default: exact matching of these curated site effects with the sites in the data 
* can be too strict: nearby sites generally tend to have similar effects
* this parameter lets the analysis use a determined inaccuracy while looking up site effects
* increasing this parameter will increase the result network size by reducing the portion of phosphorylation sites with unknown effect
* as for site matching, increase in coverage at the cost of increase in false positives

---

### Appendix - CausalPath user-tunable parameters

3. <u>Gene focus</u>

* this parameter lets the analysis use a subset of the literature relations, focusing on the neighborhood of certain proteins indicated by their gene symbols, hence reducing the number of tested hypotheses
* may be useful in two ways
   * removes irrelevant parts of the prior relations, providing a means of complexity management
   * may increase statistical power for differential abundance detection by decreasing the total number of tested peptides


---

### Appendix - CausalPath user-tunable parameters

4. <u>Generation of a data-centric causal network</u>

* CausalPath result networks are gene-centric: genes are represented with nodes, and all other measurements related to a gene are mapped on the gene’s node
* when a data row can map to multiple genes, however, this creates a redundancy
* e.g. the problem exists for mass spectroscopy when a phosphopeptide can be resolved to multiple homologous proteins, if their sequences are identical around the phosphorylation site
* alternative view: CausalPath can generate data-centric views where nodes represent data rows unresolved to particular proteins
    * this view does not support mapping other available omics data onto the network
    * the relations are duplicated when more than one data of the same gene can be explained by the same relation


--

### Appendix - CausalPath user-tunable parameters

5. <u>Protein activity</u>

* the users can insert their own hypotheses as to whether a protein is activated or inhibited in the case of a comparison-based analysis
* the input has to be a gene symbol associated with a Boolean parameter indicating the hypothesized direction of activity change

6. <u>Adjusting phosphopeptide measurements with total protein</u>

* When both phosphoprotein and total protein measurements are available in a study, an optional adjustment can be done to the phosphopeptide values to reflect their relative abundance to the total protein, before applying CausalPath

--

### Appendix - CausalPath user-tunable parameters

7. <u>Data type for expressional targets</u>

* default: CausalPath restricts its logical reasoning within proteomic data
* users can opt to use mRNA data for the targets of expressional relations
* also possible to use mRNA and protein data together by using this parameter multiple times

8. <u>Using custom resources</u>

* default: CausalPath resources embedded in its code base; subject to change with new versions of the software
* for reproducibility and customizability, the following resources can be overriden: 
    * list of priors relations using the "custom-prior-relations-file" parameter
    * list of known site effects using the "custom-site-effects-file" parameter, 
    * list of recognized HGNC symbols using the “hgnc-file” parameter. 
    



---

### Appendix - Context (Barsi et al. 2021)

<span style="font-size:20px;">


* **knowledge-driven** methods use curated lists of gene sets/pathways and use statistical methods for overrepresentation/ enrichment analyses
  * more appropriate for hypothesis generation
  * in most cases these gene sets too general to identify real causal information from data
<br>
* **data-driven** methodologies use ML to predict biological phenotypes
  * good predictive performance but limited ability for generalization and mechanistic insight
<br>
* **causal reasoning** tools connect prior-knowledge networks with gene expression or proteomics data
    * can identify contextualized, sample-specific signaling network alterations and thus causal effects explaining the observed data
    * can be used for hypothesis generation
    * future benchmarking needed

</span>


---

### Appendix - Context (Barsi et al. 2021)



![width:800px center](pictures/barsi_fig1.png)


---

### Appendix - CausalPath overview (Barsi et al. 2021)

<span style="font-size:20px;">

* CausalPath uses kinase/phosphatase— substrate and transcription factor—regulated gene relationships from the Pathway Commons database to create **graphical patterns** 
<br>
  * graphical patterns = **causal associations** (e.g."KinaseA is active when phosphorylated on site P1. Active KinaseA phosphorylates ProteinB on site P2.") 
<br>
* graphical patterns are matched with **measurements** (e.g."KinaseA is phosphorylated on site P1, and ProteinB is phosphorylated on site P2")
<br>
* which yields **causal conjectures** like "KinaseA phosphorylates ProteinB in the given dataset", identifying the potential causal way of signaling.
<br>
* **statistical significance** assessed with a data label permutation-based approach


</span>

---

### Appendix - CausalPath result highlights (Barsi et al. 2021)


<span style="font-size:20px;">


* importance of using the **correct type of prior knowledge with the corresponding omics modality**: 
    * using gene regulatory networks with proteomics data: inferred causal networks not statistically significant; significant causal associations when using same prior-knowledge network with gene-expression data

* higher abundance of transcriptomics datasets but used prior-knowledge networks defined on the level of protein activities (pathways) in most cases 

* association between gene expression and protein abundance/activity can be modest, using gene-expression data with pathway networks can lead to incorrect interpretation of the results

* importance of using **matching prior-knowledge networks and data** (gene regulatory networks with transcriptomics; signaling networks with (phospho)proteomics)

* **correct integration of different types of prior-knowledge networks** and data types also promises to identify causal associations in multi-omics datasets

* benchmarking bottleneck: lack of high-quality data where causal associations are already known (perturbation data; but off-target effects can complicate the evaluation)

</span>

