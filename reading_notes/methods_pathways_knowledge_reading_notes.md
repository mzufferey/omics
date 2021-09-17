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

