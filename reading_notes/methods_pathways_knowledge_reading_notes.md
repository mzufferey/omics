### Barsi et al. 2021

Knowledge-driven methods use in most cases extensive, curated lists of gene sets form connected biological processes or pathways and use statistical methods (with or without explicit pathway information) to find overrepresentation/ enrichment of these gene sets in biolog- ical datasets.

they are more appropriate for hy- pothesis generation. However, in most cases the used gene sets are too general to identify real causal information from data. 

, data-driven meth- odologies, including machine learning models, focus on predictive performance.

 predictive models can be important in different fields of biology from drug discovery to patient stratification. Also,

‘‘causal reasoning tools’’ connect prior-knowledge networks (like signaling pathways or gene regulato- ry networks) with genome scale gene expression or proteomics measurements and use statistical tools to identify contex- tualized, sample-specific signaling network alterations and thus causal ef- fects explaining the observed data

CausalPath uses kinase/phosphatase— substrate and transcription factor—regu- lated gene relationships from the Pathway Commons database to create graphical patterns.

These graphical patterns are causal associations like ‘‘KinaseA is active when phosphorylated on site P1. Active KinaseA phosphorylates ProteinB on site P2.’’ These kinds of graphical patterns are matched with measurements like ‘‘KinaseA is phosphorylated on site P1, and ProteinB is phosphorylated on site
P2’’, leading to causal conjectures like ‘‘KinaseA phosphorylates ProteinB in the given dataset’’, identifying the potential causal way of signaling. 

tests the statistical significance of the derived results using a data label permu- tation-based approach.

highlight the importance of using the cor- rect type of prior knowledge with the cor- responding omics modality.

When they used gene regulatory networks with prote- omics data, the inferred causal networks were not statistically significant, while us- ing the same prior-knowledge network with gene-expression data resulted in sig- nificant causal associations

given the higher abundance of transcriptomics datasets (compared to phosphoproteomics, for example), gene expression data are more frequently used in modeling studies. How- ever, the used prior-knowledge networks are defined on the level of protein activities (pathways) in most cases. 

association between gene expression and protein abundance/activity can be modest, using gene-expression data with pathway networks can lead to incorrect interpretation of the results

.suggest the crucial impor- tance of using matching prior-knowledge networks and data, like gene regulatory networks with transcriptomics and signaling networks with (phospho)prote- omics. Correct integration of different types of prior-knowledge networks and data types also promises to identify causal associations in multi-omics datasets.

A bottleneck for this benchmarking is high- quality data where causal associations are already known. For this purpose, perturbation data (where the general cause of changes is given by the used perturbation, i.e., drug, genetic manipula- tion etc.) looks most suitable,9 but off- target effects of perturbations (drugs, small interfering RNA [siRNA]) can compli- cate method evaluation. Nevertheless,





##### Huckstep et al. 2021

Literature-derived, pathway-oriented databases such as the Kyoto Encyclopedia of
Genes and Genomics (KEGG) (Kanehisa & Goto, 2000), Reactome (Joshi-Tope et al., 2005), WikiPathways (Slenter et al., 2018), and the SIGnaling Network Open Resource (SIGNOR) (Perfetto et al., 2016) have

study comparing PPI databases found 375 resources (Bajpai et al., 2020). 

Resources such as the Human Protein Reference Database (HPRD) (Peri et al., 2004; Keshava Prasad et al., 2009), the Biological General Repository for Interaction Datasets (BioGRID) (Stark et al., 2006), and the Search Tool for Retrieval of Interacting Genes/Proteins (STRING) (Szklarczyk et al., 2019) are incredibly useful resources that enable the analysis of interaction data in various contexts. PPI databases tend to be quite large as a result of their derivation from high-throughput experimental data. 

Popular phosphorylation focused databases such as PhosphoSitePlus (Hornbeck
et al., 2015), PHOSIDA (Gnad et al., 2007), Phospho.ELM (Diella et al., 2004) and qPhos (Yu et al., 2019) host many phosphosites and act as repositories for both low and high-throughput data.

qPhos is amongst the largest such databases and contains quantitative information for almost 200,000 non-redundant phosphorylation sites as well as the cell-type and temporal information

An extensive review of phosphoproteomics resources can be found here (Savage & Zhang, 2020). 

Moving away from databases acting as phosphosite repositories, many databases holding signalling information in the form of a K-S networks exist such as RegPhos (Huang et al., 2014), PhosphoNet (Safaei et al., 2011) and Phosphonetworks (Hu et al., 2014). Though they are limited to K-S interactions, they have been proven useful in providing phosphoproteomics mechanistic insight either on their own or while integrated into other databases (Rohrs et al., 2018; McGuire et al., 2017; Tong et al., 2019).

how well equipped each database is to uncover mechanisms of phosphorylation and connect phosphorylation events to a signalling network or pathway. The pathway databases we compare here are Reactome, KEGG, WikiPathways and SIGNOR. The PPI databases we compare are BioGRID and HPRD. PhosphoSitePlus was also included as it contains the largest subset ofphosphoproteomic signalling-related proteins. We first analysed the proteome coverage of each database followed by the phosphorylation coverage of a subset of the above databases. Next, we explored the consistency between the database’s phosphorylation annotations and the amino acid residue found in UniProt’s(Bateman et al., 2017) protein sequence. Finally, we assessed the capability of each database in mapping experimental phosphoproteomics datasets of varying backgrounds.



##### Rubel et al. 2020

existing methods suffer from low recall in recovering protein interactions in ground truth pathways, limiting our confidence in any new predictions for experimental validation. We present the Pathway Reconstruction AUGmenter (PRAUG), a higher-order function for producing high-quality pathway re- construction algorithms. PRAUG modifies any existing pathway reconstruction method,



While these resources are useful for exploration, they are not
intended to make predictions about new proteins and interactions that may be associated with a pathway of interest. Another



the Pathway Reconstruction AUGmenter (PRAUG), a
higher-order function which maps any pathway reconstruction method to an augmented method that improves pathway recon- struction performance. PRAUG is designed based on the observa- tion that pathway reconstruction methods typically perform well when predicting the proteins in a pathway. PRAUG takes as input a pathway reconstruction method and provides a method which uses a traversal on the input method’s predicted nodes to explore pro- tein interactions. Despite

pathway reconstruction methods as inputs for PRAUG

Overview of PRAUG. Pathway reconstruction methods take as input an interactome ?? and a pathway comprised of sources (blue diamonds) and targets (orange rectangles). Given a pathway reconstruction method M, PRAUG defines a new method fM that calls M and performs a traversal on the resultant node set.

The Pathway Reconstruction Problem is solved by recovering both proteins and interactions that connect receptors to transcriptional regulators, and previous work has shown that recovering interac- tions is a challenging task [22]. However, current pathway recon- struction methods are relatively successful at recovering pathway involved proteins.

reconstruction methods can recover the pro- teins in a pathway; however recovering the interactions is more challenging.



PRAUG is a straightforward traversal-based algorithm, and it may not be immediately clear why this approach can be so successful. In this section we justify PRAUG’s depth-first traversal, provide empirical upper bounds on the Pathway Reconstruction Problem given methods which assume a traversal paradigm, and illustrate the effect of parameter selection on PRAUG methods compared to their original counterparts



### Wong et al. 2021



Biofactoid (biofactoid.org) software suite, which crowdsources structured knowledge in articles from authors. Biofactoid is a web-based system that lets scientists draw a network of interactions between genes, their products, and chemical compounds and employs smart-automation to translate user input into a structured language using the expressive power of a formal ontology. 

In contrast, there are few efforts and little technology to support direct submission by authors of biological pathway information and related knowledge in computable form

Biofactoid (biofactoid.org), a web-based software system that
empowers authors to capture and share structured human- and machine-readable summaries of molecular-level interactions described in their publications

Structured knowledge newly acquired in this way becomes part of the global pool of pathway knowledge and can be shared in resources such as Pathway Commons (Rodchenkov et al., 2020), Network Data Exchange (NDEx) (Pratt et al., 2015), and STRING (Szklarczyk et al., 2017), to enhance information discovery and analysis



Biofactoid enables molecular-level detail of biological processes reported in articles
to be shared in a structured format accessible to humans and computers. Interactions (including binding, post-translational modification, and transcription/translation) involving molecules of various types (proteins, nucleic acids, genes, or chemicals, e.g., metabolites and drug compounds) can be represented. 

Biofactoid data represented in the formal BioPAX data
exchange language is integrated with external pathway and interaction knowledge using technologies for BioPAX processing and analysis, including Paxtools (Demir et al., 2013) and cPath2 (Cerami et al., 2006). This enables author-contributed information to be more easily searched, visualized and analyzed across different data sources.



### Babur et al. 2021 - SUPP

<u>Causality</u>

Both **pathway inference** and **pathway extraction** are closely related to **formal notions of causality inference**,
specifically to the **Suppes’ and Pearl’s probabilistic formulations**. 

A **probabilistic causal relationship** between two events, say from event A to event B, indicates t**he probability that B depends on the status of A**, as described by Patrick Suppes [2]. 

While using this notion of causality can generate predictive models, i**t does not tell if A may cause B**. For instance, there can be an event X that is causing both A and B, and this will still satisfy Suppes’ formulation. 

To **make the model predictive under an intervention scenario**, Judea Pearl provided a reformulation: **perturbing the status of A will change the probability of B** [3]. 

We follow Pearl’s notion and detect **mechanism-based evidence for activity change of one protein may affect the abundance of a specific peptide from another protein in pathway databases**, as described in the next section and in Supplementary Information.

<u>Derivation of prior relations from detailed mechanistic pathways</u>

* BioPAX-pattern framework
* manual definition of **12 BioPAX patterns to capture potentially causal binary relations that involve phosphorylation and expression of proteins**



<u>Algorithm for selection of explanatory subset of prior relations</u>

Using the extracted causal priors, CausalPath determines if there is sufficient proteomic data that indicates differential activity of that prior

Logical equations that check conditions of causality



Step 2 of CausalPath workflow: **causal conjecture generation** 

* for comparison-based analysis: integrates 4 types of information such that the prior information forms a causal bridge between a pair of observed proteomic changes
* each type of information has a fixed number of possible values
* valid combinations of these values are detected using Eq. 1



it is important to not to assume simple phosphorylation/ dephosphorylation when reading CausalPath graphs as they can be more complex indirect mechanisms

CDK1 has 2 measured features: total protein (t) and phosphorptein (p). The graph does not show which relation is related to which feature, but in this case, it is inferrable.

<u>Previously published methods for pathway analysis of proteomic datasets</u>
There is no method comparable to CausalPath for its ability to identify causal relations from pathway databases that can explain given proteomic datasets. 

There are, however, methods developed for other forms of pathway analysis for proteomics. These methods generally use prior information in the form of networks stripped from mechanistic details, and **aim to build a network structure that most fits to the profiling data** at hand using an optimization method. 

During the development of CausalPath, instead of a network optimization, **we intentionally focused on better usage of prior information by processing mechanistic details in pathway databases and using them in logical reasoning for causality**. CausalPath, in that sense, is not competing with these alternative methods, but it is complementary to them. CausalPath can easily be paired with any other network optimization method for further complexity management and for using strong priors in network optimization. 

**Temporal Pathway Synthesizer** (TPS) TPS uses PPI and kinase-substrate networks in the process of inferring signaling relations from temporal
post-perturbation proteomic data [1].

**PARADIGM** is one of the earliest pathway analysis methods developed originally for RNA expression and copy number variations and later extended to other data types including proteomics [5]. It uses the pathway models from NCI Pathway Interaction Database (PID), converts them into a factor graph, and predicts each entity’s activity level using an expectation-maximization algorithm

**pCHIPS**
This is a network propagation method for proteomic and other data, based on the TieDIE [6] algorithm, where the purpose is to link differentially active kinases (indicated by proteomic data) to the differentially active transcription factors (indicated by RNAseq measurements of targets) [7]. Proteomic changes on the kinases are propagated downstream, differential transcription factor activities are propagated upstream on the signaling network, and the overlap is identified as a possible linking path or combination of paths

**SigNetTrainer** [8] is a set of algorithms that score the fitness of molecular readouts to a given set of per- turbations and a given directed and signed interaction graph, and solves several interesting problems using integer linear programming, such as finding an optimal subgraph that is most consistent with the measure- ments, or finding minimal set of new relations that will make the network and the measurements consistent.

**PHONEMeS** [(9)] builds models in the form of Boolean networks that best fit to a given set of phospho- proteomic perturbation data.

**Chasman** et al. demonstrate their network inference method on identification of yeast adaptive pathways to NaCl stress [10]. They compile a background network with directed and undirected relations, identify a set of genes/proteins by differential expression, phosphoproteome changes and stress fitness contribution, and find optimal paths from signaling proteins to gene regulation proteins employing integer programming (IP)

**PhosphoPath and PTMapper** Both methods are implemented as a Cytoscape plugin to visualize kinase-substrate relations on the protein-
protein interaction (PPI) network 

**PCST**  maps proteomic and transcriptomic data on the proteins on a weighted PPI and protein-DNA interaction network, then identifies a minimal subnetwork that connects the mapped molecules, prioritizing the most reliable connections [13].

**PHOTON**
This method maps proteomic data from a perturbation study onto the proteins on a weighted PPI network, then calculates a score for each protein based on the weighted average of the observed proteomic changes on its neighbors on the network [14].



These patterns can be extended to either **cover new data sources or to capture new relationship types**. 

* For new data sources, it is very likely that the listed patterns will be applicable as is, however not guaranteed. When needed, the graph structures in the new data source should be investigated and new graph structures should be captured by adding new patterns to the framework. 

* The procedure is similar for capturing new types of relationships. Developers need to study the structure
  of the existing BioPAX data to understand in what forms the information is encoded. Then



<u>Phosphorylation</u>

* Pattern 1:  the regulator proteins activate a Conversion that adds phosphorylation to a protein, or inhibit a Conversion that removes phosphorylation from the downstream protein

* Pattern 2: the regulator proteins are not modeled as regulators of a Conversion, but modeled as inputs and outputs of the same Conversion

* Pattern 3: some proteins in a complex are phosphorylated, and the input complex is also designated to be the controller of the Conversion	 

* Pattern 4: the controller transmits its effect through a controlling small molecule

  

<u>Dephosphorylation</u>

* Pattern 1: the regulator proteins activate a Conversion that removes phosphorylation from a protein, or inhibit a Conversion that adds phosphorylation to the downstream protein
* Pattern 2:  the regulator proteins are not modeled as regulators of a Conversion, but modeled as inputs and outputs of the same Conversion
* Pattern 3: some proteins in a complex are phosphorylated, and the input complex is also designated to be the controller of the Conversion	
* Pattern 4:  the controller transmits its effect through a controlling small molecule

<u>Upregulation</u>

* Pattern 1: the regulator activates a TemplateReaction
* Pattern 2: a Conversion is used instead of a TemplateReaction

<u>Downregulation</u>

* Pattern 1: the regulator inhibits a TemplateReaction

* Pattern 2: a Conversion is used instead of a TemplateReaction.

### Babur et al. 2021

a computational method to infer causal mechanisms in cell biology by analyzing changes in high- throughput proteomic profiles on the background of prior knowledge captured in biochemical reaction knowl- edge bases

automate scientific reasoning processes and illustrates the power ofmapping from experimental data to prior knowledge via logic programming

hundreds of pathway and interaction databases (pathguide. org). 

data-driven inference approach leverages the recent developments in proteomics and other molecular tech- nologies to directly infer graphical models, ab initio, from high- throughput measurements of controlled perturbations and natu- ral variation.1–3 

more established approach is to compile extensive, interconnected pathway models through the curation of reac- tions based on carefully designed low-throughput controlled ex- periments. 

Both the classic and the data-driven inference approaches
have inherent limitations. The classic curation approach uses well-validated fragments of knowledge, but these are extracted from a heterogeneous set of contexts, perturbations, condi- tions, and even organisms. The resulting models, even when carefully restricted to a particular context, are not well suited to making predictions. The purely data-driven inference ap- proaches, on the other hand, create context-specific, predic- tive models, but they do not scale in terms of statistical power as the model space is exponentially larger than the observ- able space. 

get help from prior knowledge when the pertur- bations in the data are not sufficient to decide between alterna- tive models.4–10 The methods that use this strategy, however, use prior knowledge in a reduced form, such as simple interac- tion networks, omitting mechanistic details and their logical harmony with the new data

CausalPath, which uses the rich semantics of curated pathway knowledge, including the type of mechanism, the direction, signs of effect, and post-translational modifications. The inferred mechanisms are falsifiable hypotheses that can be experimentally interrogated.

CausalPath maps proteomic profiles to curated human path-
ways from multiple resources that are integrated into the Pathway Commons database,11 detects the potential causal links in the pathways between measurable molecular features using a graphical pattern search framework, and identifies the subset of the causal links that can explain correlated changes in a given set of proteomic and other molecular pro- files. 

These explanations are presented as an intuitive network with links to the detailed prior knowledge models and the related literature to create a powerful exploration and analysis platform

Themethod takes into account hun- dreds of thousands of curated mechanisms, which would be infeasible to do manually.

<u>The CausalPath workflow has two main steps:</u> 

(1) detection of causal priors from pathway databases, performed once and reused in multiple analyses, and 

(2) matching causal priors with supporting correlated changes in the analyzed data, performed for every analysis. 



Existing kinase-substrate databases and transcription fac-
tor-target databases are valuable sources for causal priors, but they capture only a small part of the known biology; hence, they are limited for comprehensive causal reasoning



other databases that take a more detailed modeling approach for biochemical processes, such as Re- actome. The Pathway Commons database provides integra- tion of such detailed models collected from publicly available resources in the format of the BioPAX modeling language



Detailed process models provide a great opportunity to identify causal relations between the molecular measure- ments, but they require sophisticated algorithms to reason over them. 

To detect the causal prior relations, i.e., structures that
imply causal relationships between proteins in the Pathway Commons database, we used the BioPAX-pattern software12, and manually curated 12 graphical patterns

. Each graphical pattern captures the control mechanisms over either a phosphorylation of a protein or the expression of a gene. 

We assessed the overlap of these prior relations with the ‘‘canonical pathways’’ gene sets in MSigDB to understand its coverage. This collection has 2,815 gene sets curated from the databases BioCarta, KEGG, NCI-PID, Reactome, and WikiPath- ways. 



a ‘‘**causal conjecture**’’ as a pairing of a causal prior
with supporting measurements in the molecular dataset that together declare that ‘‘one molecular change is the cause of another molecular change.’’



Alternative to the above controlled perturbation setting where we compare control and test samples, a causality network
can be generated based on the correlations of measurements in an uncontrolled cohort—as is common in cancer biology. In

e.g. “Measured peptide levels of NCK2-pT110 is positively correlated with the peptide levels of MAPK14-pS180-T182”, for a **correlation-based causality hypothesis**

The **change** here can be detected in two different forms based on the experimental setting: it can be up/downregulation for individual features in a ‘‘test versus control’’ comparison setting, or it can be positive/negative cor- relations applying to pairs of features in an uncontrolled study, as is common in cancer biology. We

We call an analysis in the former setting ‘‘**comparison-based**’’ and the latter ‘‘**correla- tion-based**.’’ 



To formalize and generalize the example of causal conjecture
detection in comparison-based analysis, we can formulate it with a ternary logical equation

**comparison-based**
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



**correlation based**
$$
corr_{source, target}\oplus(e_{source} \oplus s_{relation}) = true
$$

* logical representation of the sign of the correlation replaces the $c$ terms
  * $true$ in case of positive correlation
  * $false$ in case of negative correlation



<u>equations in the code</u>

* .sign and .effect properties of variables that takes values -1, 0 and 1, corresponding to false, unknown and true, respectively
* multiplication of these integer values and checking the result value is an alternative formulation to the original logical equations where ternary XOR (⊕) operator is used
* when the RNAseq data is used at the targets of expressional control relations, the total protein measurements (mtt in the pseudocode) are replaced with the RNAseq measurements of the target (mrt)





In addition to the logical check by these equations, we limit:

1. the phospho regulations (phosphorylation and dephos- phorylation) to the explanation of phosphoprotein changes and
2. the expressional regulations to the explanation of total protein changes (and optionally mRNA changes)



On top of the logic-based detection of causal interactions,
we provide two types of statistical measurements to increase the interpretability of the results

1) ‘‘**Network-size test’**’ checks if the correlated changes align with the causal priors in general, which is indicated by **a larger number of interactions in the re- sults than would arise by random chance**, which we test by data label randomization
2) ‘‘**Downstream-size test’**’ checks if a protein on the network has **more downstream targets in the re- sults than expected by chance** using the same randomization procedure. 

Significant values from these two tests provide additional evidence suggesting the data are shaped by the priors or that a protein has an influence on the significant num- ber of targets, respectively, which consequently increases our confidence in the results



Experimental data reveal protein features that change in coordination, and CausalPath automates the search for causal explanations in the literature.



context-specific correlations are derived from the data and causality is derived from the literature. 

Compared with the methods that infer causal- ity from data through mathematical modeling (pathway infer- ence), this method has a much wider application area. Pathway inference methods have a potential to offer more complete re- sults, but they require numerous perturbations and/or time points in the experiments, whereas the pathway extraction strat- egy is applicable to any simple comparison, or a set of profiles from a cohort with some variance to explain.

CausalPath is a great resource for high-confidence priors, which can inform the hybrid pathway inference methods that benefit from prior data



 Causality-focused pathway extraction approaches provide a means of transferring knowledge between contexts. CausalPath does this by detecting variation patterns of proteomic abundances and detecting their consistency with prior knowledge. A limitation of this approach is its dependency on observable variance; therefore, it cannot identify a signaling relation that does not significantly vary across the compared samples.



The added value that our method brings to the field of pathway extraction is three-fold:

(1) interpretation of complex mechanistic pathway models, 

(2) site-specific evaluation of phosphoproteo- mic measurements, and 

(3) a logical test for causality between measurements.



When all these are combined, pathway extrac- tion becomes a useful tool for ‘‘mechanistic model building.’’ 

We expect future research will take these ideas further, poten- tially addressing these two challenges: 

(1) instead of a binary evaluation of causality (between two proteins), n-ary systems can be developed, and 

(2) instead of a binary classification of protein modifications as activating and inhibiting, a site can be more accurately mapped to a distinct subset of activities of the protein.



t evidence of known phospho regulations is more consistently observed in the proteomic data compared with the known expression regulations. 



One major limiting factor in this analysis is a large number of pro- tein phosphorylation sites whose functions are not known; hence, their downstream cannot be included in the causality network. We are actively working to mine these data from the literature using natural language processing tools.28 In the mean- time, CausalPath reports those sites with an unknown effect that also have significant change at their signaling downstream. Users have the option to review this list of modification sites and manually curate them to increase the coverage of the analysis.



CausalPath can be applied to the results of any proteomic and phosphoproteomic experiments to identify differential signaling that is supported by literature knowledge. To use CausalPath, the measurement values need to be comparable (normalized) and need to be associated with related gene symbols, and phos- phopeptide measurements need to specify the phosphorylation sites with respect to their canonical UniProt sequence



, CausalPath can be run locally as a Java application using the sources at
https://github.com/PathwayAndDataAnalysis/
causalpath. This repository additionally includes examples running Causal- Path from R and Python.34 The generated result networks can be visualized with ChiBE,29,30 as well as by uploading analysis output folders to the Causal- Path web server at causalpath.org.



As molecular datasets get poorer in perturbations, they have to depend more and more on priors for pathway inference. But this scheme has one problem: priors used in pathway inference methods are generally simple networks such as protein-protein interactions—they do not inherently suggest any causal relationship between measured molecules. To extend pathway inference to more common proteomic datasets using priors, we need high-confidence priors that contain causality within. In this study, we propose a method for generating causal priors from mechanistic pathway relationships curated from scientific literature, then we use them for the pathway analysis of publicly available perturbation-poor proteomic datasets. Our approach reduces the pathway inference task into identifying the causal priors that can explain correlations in a given molecular dataset, which we distinctly name as pathway extraction.

—CausalPath—uses curated mechanistic human pathways from multiple resources that are integrated into the
Pathway Commons database6, detects the causal links in the pathways between measurable molecular features using a graphical pattern search framework, and identifies the subset of the causal links that can explain correlated changes in a given set of proteomic and other molecular profiles. We

This approach, in essence, mimics the literature search of a biologist for relationships that explain their data.	



-> it generates extremely useful and falsifiable hypotheses that are unavailable otherwise.

-> Since this process systematically considers hundreds of thousands curated mechanisms, it also is more comprehensive, unbiased, and more consistent in terms of the generated hypotheses

Alternative to the above controlled perturbation setting where we compare control and test samples, a causality network
can be generated based on the correlations of measurements in an uncontrolled cohort—as is common in cancer biology.



The correlation-based causal network provides hypotheses
for the signaling network parts that are differentially active across samples, but it does not indicate which parts are activated together or whether they align with previously defined molecular subtypes. The original TCGA study on HGSOC samples iden- tifies four molecular subtypes based on RNA expression, termed as immunoreactive, differentiated, proliferative, and mesen- chymal.19 To understand if we can gain mechanistic insight into the previously defined subtypes, we compared each sub- type to all other samples using a t test with Benjamini-Hochberg FDR control on measurements, but we were unable to generate substantial results within a 0.1 FDR threshold, probably due to the large proportion of missing values in the phosphoproteomic dataset combined with the loss of statistical power due to smaller cohort size for each subtype. Then we tried to constrain the search space with the neighborhoods of some of the genes with differential measurements, and relax the FDR threshold at the same time for further exploration. Six SRC family kinases (SFKs) have proteomic evidence for activation in the immunore- active subtype; hence, we limited the search to the neighbor- hood of SFKs (SRC, FYN, LYN, LCK, HCK, and FGR), set the FDR threshold to 0.2 for phosphoproteomic data, and identified 27 relations



we detected that the
breast cancer phosphorylation network is significant in size (p < 0.0001), while the expression network is not (p = 0.5521), suggesting that the known phosphorylation relations have a much higher impact on the proteomic correlations than known expressional relations



we compared the PAM50 expression subtypes of breast
cancer to see if we could get causal explanations of their prote- omic differences. We were again challenged by decreased sam- ple sizes and missing values, but we detected that luminal A and luminal B subtypes have significant differences from the basal- like subtype.

ESR1 is significantly more active in luminal breast cancers, suggested by both its protein levels and the changes in its downstream (Figure





##### Triantafillou et al. 2017

we apply state-of-the art causal discovery methods on a large collection of public mass cytometry data sets, measuring intra- cellular signaling proteins of the human immune system and their response to several perturbations. We show how different experimental conditions can be used to facilitate causal discovery, and apply two fundamental methods that produce context-specific causal predictions.

r, causal knowledge confers more information as compared to correlation



Computational causality has developed a language to describe, quantify and reason with causal claims. The
most common framework of computational causality is causal Bayesian networks (CBNs), that use a simple assumption to connect causal relationships to associative patterns1. CBNs use directed acyclic causal graphs to describe the causal relationships and connect them to associations expected to hold or vanish in the joint proba- bility distribution. Causal effects can also be computed using CBNs using do-calculus, a formal system for causal reasoning that includes an operation for interventions1. Algorithms for automatically identifying CBNs from limited or without experiments have also been proposed

In this work, we attempt to discover novel causal relationships from a large collection of public mass cytome-
try data of immune cells perturbed with a variety of compounds.

We find that (a) results are highly consistent on data sets that include different donors, experimental
cell-stimulation time-points, or cell types (b) different causal methods often disagree with each other and with known pathways, (c) validity in experimental data is inconclusive. Our

These results indicate that (a) de novo discovery of causal pathway relations is still a challenging task for current causal discovery methods, despite the previous positive results4 and (b) current causal discovery methods do identify reproducible findings across similar data sets. 

Causal discovery makes assump- tions on the nature of causality that connect the observable data properties (i.e., the joint probability distribution of the observed variables) to the underlying causal structure. The most popular causal assumption, famous for inspiring Bayesian networks, is the Causal Markov (CM) assumption. The CM assumption states that every vari- able is independent of its non-effects given its direct causes

In computational causality, the causal structure of a set of variables is often modeled directed graphs. Causal
Bayesian Networks are the most popular framework for causality. In CBNs, causal relationships are modeled using directed acyclic graphs: A directed edge from X to Y denotes that X causes Y directly in the context of measured variables; no measured variable included in the model mediates the relationships. The causal structure is assumed to be acyclic. In addition, no pair of variables in the graph can have an unmeasured common cause. Several exten- sions try to relax these assumptions7–9. For simplicity, we present the theory of causal discovery using conditional independence tests using CBNs. 

Given the causal graph, one can easily identify direct causes and non-effects of a variable. The CM assumption
connects a given causal graph with a set of conditional independencies (CIs)

Conditional independencies can be tested in the data using appropriate tests of independence. Computational
causal discovery tries to reverse engineer the causal graph using tests of independence: given a pattern of condi- tional independencies, a causal discovery algorithm typically tries to identify the causal graph that is associated to this pattern through the CM assumption. To do so efficiently, algorithms assume that all observed independencies in the data are a product of the causal structure, rather than being “accidental’’ or “fine-tuned’’ properties of the model parameters. This assumption is known as the Faithfulness assumption2, and it is employed by most causal discovery algorithms. E

Most of the times, CM and faithfulness do not uniquely associate a CI pattern with a single causal model.



The conditional independence pattern can be identified using appropriate statistical tests. This method, called Local Causal Discovery (LCD), has been proposed for automatically mining causal relationships in large data sets14. This simple structure has also been successfully exploited in genome-wide association studies (GWAS), where Mendelian randomization is used as a known uncaused entity



LCD tries to identify this structure in data, testing for all pairwise dependencies and the conditional inde-
pendence of A and T given S. 

Compared to trying to reconstruct the entire network of measured variables, considering only a small subset
of variables at each run has many advantages for de novo causal discovery. Causal discovery algorithms are noto- riously prone to error propagation17, and a single error in a conditional independence test can affect seemingly remote parts of the learnt network. Employing a local approach allows testing all possible conditional independ- encies, and minimizes the probability that an independence test of unrelated variables will affect the output. Thus, error propagation does not affect the precision of CLCD

We have so far presented CLCD assuming the “true” causal network is a CBN, i.e. there are no confounders
and no feedback loops. However, violation of these assumptions does not have an effect on the validity of causal relationships identified with CLCD. Extensions of CBNs use bi-directed edges to model the presence of a con- founder between two variables. Graphs that also include bi-directed edges are called mixed graphs. Conditional independencies implied by the CMC can be identified in mixed graphs using an extension of d-separation. the equivalent of d-separation for mixed graphs. T

In mass cytometry experiments, the presence of an activator can be modeled as an external, binary variable,
that is set by the experimenter and is not influenced by protein phosphorylation within the cell. All causal models where the activator is caused by a phosphorylated protein can then be excluded, and the conditional independ- ence can only be explained if the mediation of S is causal. For

CLCD only uses a small portion of the available data, as it disregards all but the lowest inhibitor dosage. Data from different dosages of the same inhibitor correspond to different distributions and cannot be pooled together for CLCD. Some methods exist for integratively analyzing data from multiple experiments18–20. The experiments are required to be surgical interventions18,19 (meaning that the targets of the experiments must be known and their values completely set by the experimental conditions). In



in the method BACKSHIFT21, different experimental conditions (or environments) with unknown targets are used to uncover causal structure. Different experimental conditions are modeled as so-called
“shift interventions”, meaning that the values of unknown targets of each intervention are shifted by realizations from a random variable modeling this intervention. Thus, the experiments are not required to be surgical inter- ventions. In contrast, the underlying causal structure is assumed to persist in the experimental condition, affect- ing the system in addition to the induced intervention effect. The method then exploits the presence of different environments to identify the causal structure of the measured variables, as well as the location and strength of the interventions in each experiment. To this end, the method assumes a linear causal model that possibly includes cycles and confounders. An

In this work, we examined the performance of two causal discovery methods on a collection of mass cytom-
etry data. The methods are based on fundamental causal principles, and use multiple data sets and/or different experimental conditions to increase robustness. However, we found that even basic methods often disagree with each other and with background knowledge (such as the KEGG pathways).





##### CausalPath in Satpathy et al. 2021

*Furthermore, knockdown of these tumor biomarker proteins reduced fitness across 16 LSCC cell lines (https://depmap.org/), suggesting crit- ical roles in key cellular transformation and proliferation pro- cesses (Figure 7C; Table S7).*

<u>CausalPath analysis</u> CausalPath (Babur et al., 2018) searches for known biological mechanisms that can explain correlated proteomic changes in terms of causal hypotheses. We set CausalPath parameters to compare tumors and NATs with a paired t test, used 0.1 as FDR threshold for proteomic change significance and network significance, and detected 5917 potential causal relations between proteins. We repeated the same analysis for each NMF subtype separately and identified 4378 (basal-inclusive), 5334 (classical), 3048 (EMT-en- riched), 3744 (inflamed-secretory), and 4332 (proliferative-primitive) relations. We used these CausalPath network results in the preparation of Figure 7C, identifying potential upstream regulators of oncogenic phosphoproteomic changes. Here an oncogenic phosphoproteomic change can be any of the following 4 events: increase of activating phosphorylation of an oncoprotein, decrease of inactivating phosphorylation of an oncoprotein, decrease of activating phosphorylation of a tumor suppressor protein, and in- crease of inactivating phosphorylation of a tumor suppressor protein. We used the OncoKB database for oncoprotein and tumor suppressor classification (excluded proteins that have both annotations), and used PhosphoSitePlus and Signor databases for the activating/inhibiting classification of phosphorylation sites. In the phosphorylation regulation networks, we included only the targetable regulators (activated proteins) and excluded the untargetable regulators (inactivated proteins).

**CausalPath in Wang et al. 2021**

*Application of CausalPath (Babur et al., 2018) to the protein and phosphoprotein expression data (Figure S7A, Table S5)showed upregulation of the hypoxia pathway in mesenchymal tumors, evi- denced by significant activation of multiple HIF-1 downstream targets (networkpermutationp= 0.0012). Increasedangiogenesis was also evident in mesenchymal tumors, as demonstrated by upregulation of FLT1, MMP14, ENG, and SERPINE1. We observed complex regulation of macrophage activation and po- larization through the upregulation of STAT3, ICAM1, SPI1, and CEBPB. In addition, the M1 polarization marker ARG1 showed increased expression (Arlauckas et al., 2018), along with SERPINE1 and HCK proteins, which promote M2 polarization (Kubala et al., 2018). The elevated inflammatory response in mesenchymal tumors may result in downstream activation of either hypoxia or macrophage polarization through multiple medi- ators, including LANE, IL18, and CD40 (Figure S7A)*

<u>Causative pathway interaction discovery using CausalPath</u> To discover the causative pathway interactions in our proteomic and phosphoproteomic data, we took the normalized expression of protein with < 10% missing values and phosphoprotein with < 25% missing values across all tumor and normal samples as the input to CausalPath (commit 7c5b934). We ran CausalPath in the mode that tests the mean values between test and control groups (value-transformation = significant-change-of-mean), where the test group being the tumors of one subtype and control group being the rest of the tumors. The pathway interaction discovery data source was Pathway Commons v9 (built-in- network-resource-selection = PC). Additionally, we enabled the causal reasoning if all the downstream targets of a gene were active or inactive (calculate-network-significance = true, use-network-significance-for-causal-reasoning = true, permutations- for-significance = 10000). The causative interactions with FDR < 0.05 were extracted and visualized (fdr-threshold-for-data-signif- icance = 0.05 phosphoprotein, fdr-threshold-for-data-significance = 0.05 protein, fdr-threshold-for-network-significance = 0.05). Full result tables were available in Table S5.

**CausalPath in Babur et al. 2020**

We next analyzed differential platelet protein modifications measured above with CausalPath.24,36 CausalPath computationally identifies pairs of protein phosphorylation changes with likely cause-effect relations using pathway information accumulated in Pathway Commons25 and other databases. To do this, CausalPath evaluates graph patterns represented with BioPAX,37,38 and identifies potential causal relations (“causal priors”) between site-specific phosphorylations (i.e., “change in phosphorylation of protein X at site A causes change in phosphorylation of protein Y at site B”). CausalPath then identifies parts of proteomic datasets explainable with causal priors (Figure 4A). For example, as the phosphoproteomics experiments above measure increases in p38 (MAPK14) Y182 (activating site) and STAT1 S727 (p38 substrate site) phosphorylation, CausalPath infers that actived p38 phosphorylated STAT1 S727 (Figure 4A). More details on CausalPath are described in Figure S7. From 1,887 and 2,077 significantly increasing and decreasing phosphorylation events measured over Conditions #1 and #2, CausalPath mapped 290 significant phosphorylation changes through 319 inferred relations among 148 proteins (Figure 4B; Figures S8-9) 

As summarized in Figure 4B, CausalPath identified and mapped many causally related phosphorylation events on established platelet GPVI (GP6) effectors, including SFK-mediated phosphorylation of FcR? (FCER1G), and, subsequently Syk, BTK, PLC?2 (PLCG2), PKC? (PRKCD) and p38 (MAPK14). CausalPath also inferred the activation of PKC? (PRKCA) and PKA (PRKACA) from substrate phosphorylation patterns. Other targets not yet characterized for roles in platelet function had notable changes in activating phosphorylation sites and fit into CausalPath models – including many components of Ras/MAPK signaling (i.e., KSR1, SOS1). CausalPath also noted conflicting relations in datasets, where the activation status of kinases did not match associated changes in substrate phosphorylation, including dephosphorylation of several PKA and MAPK substrates (Figures S10-11). While many significantly differentially phosphorylated proteins fit into causal contexts, the majority had no associated causal priors in Pathway Commons, including heavily modified targets such as MYCT1, DENND2C and BIN2 (Figure S12-13).
Validation

-> To examine pathway relations inferred by CausalPath, platelets were preincubated with a panel of kinase inhibitors (detailed in Table S3) prior to CRP-XL stimulation and Western blot analysis with antisera against phosphorylation sites of interest.

-> We next screened effects of inhibiting key nodes in CausalPath models (i.e., Syk, BTK, p38) as
well as more specialized targets (i.e., KSR1, SOS1, and PFKFB3, highlighted in orange, Figure 4B) over essential platelet function responses. As seen in Figure S14 and summarized in Figure 5C, pretreatment of platelets with inhibitors targeting critical GPVI effectors altered the ability of platelets to adhere to a surface of immobilized CRP-XL. Inhibitors of KSR1 (APS-2-79),39 SOS1 (NSC-658497)40 and PFKFB3 (AZ PFKFB3 67)41 also significantly limited platelet adhesion to CRP-XL-coated cover glass.

-> As summarized in Figure 5D and S15, experiments above tested 22 directed pairs of inferred signaling relations, validating 16 direct relations (green arrows), and providing evidence for an additional 19 relations in a multi-step manner (Figure S15), where inhibition of kinases and other nodes had significant effects on platelet adhesion, secretion and integrin activation. Together, these modeling and pharmacological experiments validate signaling relations and begin to place established and putative GPVI effectors together into models of platelet function

Overall, we found that several relations inferred by CausalPath were readily testable to validate
and clarify how more complex phosphorylation patterns fit into GPVI signaling and platelet function

-> We suspected that, like established GPVI effectors (i.e., Syk, PLC?2), other less well described
nodes that integrated into CausalPath models may also regulate platelet function. While

-> Secondarily, CausalPath organizes data that does not fit models established in literature,
highlighting areas for exploration (Figure S12-13). Consequently, we noted that >40 proteins associated with Rab GTPase regulation are modified in platelets in an unspecified manner following GPVI activation. S



##### Carroll et al. 2019

We hypothesized that proteomic data would complement mutation status to identify vemurafenib-sensitive tumors and effective co-treatments for BRAF- V600E tumors with inherent resistance

Linear and nonlinear regression models using RPPA protein or RNAseq were evaluated and compared based on their ability to predict BRAF-V600E cell line sensitivity (area under the dose response curve)

CausalPath software was used to identify protein- protein interaction networks that could explain differential protein expression in resistant cells.

Orthogonal partial least squares (O-PLS) predicted vemurafenib sensitivity with greater accuracy in both melanoma and non-melanoma BRAF-V600E cell lines than other leading machine learning methods, specifically Random Forests, Support Vector Regression (linear and quadratic kernels) and LASSO-penalized regression.

use of transcriptomic in place of proteomic data weakened model performance. 

CausalPath analysis of resistant cell lines CausalPath software [15] was used to identify net- worksofproteinsfromthe RPPA data set thatwere significantly enriched in the resistant cell lines (IC50 AUC < 0.2) compared to the sensitive cell lines. For analysis of predictive protein interactions, proteins with a VIP > 1 were examined (87 of the original 232 proteins met this criteria), and significant change in the mean expression of each protein/phosphorylated protein between the two groups was determined with 10,000 permutations and a FDR of 0.2 for total and phosphorylated proteins. This relaxed discovery rate is consistent with prior use of this algorithm with a constrained subset of proteins [15].

To further analyze these proteins, we next examined their involvement in cellular signaling pathways. CausalPath is a computational method that uses biological prior knowledge to identify causal relationships that explain differential protein ex- pression and phosphorylation [15]. Cell

Cell lines were sorted into sensitive and resistant groups based on IC50 AUC, and CausalPath was used to identify protein-protein in- teractions (PPIs) that explained significant changes in mean expression of the predictive total and phosphopro- teins (VIP score > 1) in the resistant cohort of cell lines

This computational method identified that the resistant subset had increased expression of EGFR and HER3- Y1289, which could be explained by the biological prior knowledge that EGFR transphosphorylates HER3 in EGFR-HER3 heterodimers 

While CausalPath identified expression patterns from PPIs, it is limited by the input proteins represented in the dataset, (i.e., it can- not find the relationship A➔ B➔ C if only A and C are measured). Because the important proteins in the O-PLS model (VIP score > 1, Fig. 3c) do not include the complete cell proteome, CausalPath could not identify a full pathway, but did identify several protein interactions in the PI3K pathway, suggesting that this pathway may also be of interest (Fig. 5a). Manual

Manual curation of 29 pro- teins in the PI3K pathway present in the RPPA dataset are shown in a heatmap in Fig. 5b, with their projections along the principal component space of the O-PLS model in Supplemental Fig. S2. The pathway curation includes receptors, adaptor proteins, and downstream signaling cascade proteins, many of which have a VIP score greater than 1 (Additional file 9: Fig. S2A bolded). Examination of the projections of phosphorylated pro- teins present from this dataset shows that the majority of them project along the negative predictive component space, indicating that elevated levels correlated with more resistant cell lines (Additional file 9: Fig. S2B or- ange). Therefore,

Therefore, through CausalPath analysis and man- ual pathway curation, we have identified that ErbB family signaling and downstream PI3K pathway activa- tion are upregulated in cell lines that are resistant to vemurafenib.





##### Niederdorfer et al. 2020

In this study, we therefore set out to investigate the use of a single signaling knowledge network to predict synergistic drug combinations for four cancer cell lines derived from gastric, colorectal or prostate cancer, by calibrating the general model to cell-specific models using their baseline activity data.

Our results show that model calibration with a protein activity profile combining information from literature-curated and omics-inferred data increases predictive sensitivity. 

Network refinements based on literature knowledge accounting for subtle biologically founded mechanisms improve the model’s predictive performance.

Model refinement was guided by biological insights that potentially underly false negative predictions. For this, the literature was searched for biological mechanisms regarding investigated drug combinations and for molecular mechanisms of single drugs that could explain single drug effects. We



##### Tsirvouli et al. 2020

a relatively large manually curated logical model can be efficiently enhanced further by including components highlighted by a multi-omics data analysis of data from Consensus Molecular Subtypes covering colorectal cancer. 

The model expansion was performed in a pathway-centric manner, following a partitioning of the model into functional subsystems, named modules. The resulting approach constitutes a middle-out modeling strategy enabling a data-driven expansion of a model from a generic and intermediate level of molecular detail to a model better covering relevant processes that are affected in specific cancer subtypes,

a tumor-data driven middle-out approach toward refining a logical model of a biological system can further customize a computer mode



***$\rightarrow$ Niederdorfer and Tsirvouli -> as a proof of principle that knowledge driven analyses yield improvement ??***

