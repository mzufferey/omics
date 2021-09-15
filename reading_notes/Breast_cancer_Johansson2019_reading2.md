an unbiased analysis of breast tumor proteomes, inclusive of 9995 proteins quantified across all tumors, for the first time reca- pitulates BC subtypes. A



poor-prognosis basal-like and luminal B tumors are further subdivided by immune component infiltration, suggesting the current classification is incomplete. Proteome-based



Genes included within prognostic mRNA panels have significantly higher than average mRNA-protein correlations, and gene copy number alterations are dampened at the protein-level; underscoring the value of proteome quantification for prognostication and phenotypic classification. Furthermore,



in-depth quantitative profile of the proteomes of 45 breast tumors, 9 represented from each of the 5 PAM50-based molecular classifications. We

the first to recapitulate the current mRNA-based molecular classifi- cations with an unsupervised analysis of whole-proteome data.





multiple layers of systems measurements collected on the same tumors, including those of mRNA expression, genome copy-number alterations, single-nucleotide polymorphisms, phosphoprotein levels, and metabolite abundances.





novel immunohistochemical biomarker candidates to more reliably stratify difficult-to-classify patients for treatment options, provide



(www. breastcancerlandscape.org).
Results



, 13,997 protein products of 12,645 genes were identified
at a 1% protein false-discovery rate (FDR) based on 248,949 identified unique peptides (Fig. 1, Supplementary Fig. 1A, B, Supplementary Data S1). The subset of 9995 proteins quantified (with a median of 12 unique peptides/protein and 24 PSMs/ protein for quantification) in each of the 45 tumors, based on gene symbol centric quantification (denoted proteins henceforth), is used for all quantitative proteome analyses (i.e.,







Unsupervised hierarchical clustering of proteome profile



the luminal B and HER2 subtypes are intermixed, indicating similarities in the molecular phenotype. The

(supported by tumor-transcript profiles of)



Gene ontology enrichment analysis reveals that proteins considered luminal markers, basal markers, or members of the HER2 amplicon, localized to the mitochondria or Golgi apparatus, related to proliferation, transcription, adipose tissue, erythrocytes, immune response, or the extracellular matrix are closely correlated and coregulated with members of their respective groups. O



PAM50 subtype agreement with correlation-based hierarchical clustering of tumor protein expression profiles considering only the 37 PAM50 gene members in the quantified proteome demonstrates the patient-stratifying information contained within the entire proteome is derived from a smaller subset





core tumor clusters (CoTC)

Core sets of tumors whose proteomes are representative of a proteome-based grouping 

defined using unsupervised clustering based on high-variance protein (n = 1334) abundance profiles

for normal-like and luminal A: assignments overlap with PAM50

basal-like tumors PAM50: divided into 2 groups

combine HER2 and luminal B

a separate group of luminal B PAM50



 (Supplementary Methods), producing six consensus core tumor clusters (CoTC) (Fig. 2c, Supplementary Fig. 4C–K). CoTC assignments overlap with PAM50 normal-like and luminal A classifications, but divide PAM50  into two groups, and combine HER2 and luminal B while maintaining a separate group of luminal B PAM50 subtype tumors (Fig. 2c, d, Supplementary Fig. 3I). Unsupervised clustering of CPTAC breast tumor proteomes5, using the overlapping high-variance proteins (632 of 1334), identifies three tumor clusters that resemble CoTC1 (basal-like), CoTC3 (luminal A), and CoTC6 (a mix of luminal B and HER2) (Supplementary Fig. 5A, B). Ofnote, the CPTAC patient cohort does not have a defined normal-like tumor subtype.



Defining proteome based tumor subtypes - CoTCs
Cluster Detection - overview In order to identify the strongest subtypes in the protein dataset, we first identified and filtered out unmodulated genes (genes with low standard deviation). From this list we then removed genes that could be associated to adipose, plasma, erythrocytes and immune pathways from the lists defined above, to obtain as pure list of tumor proteins as possible. Next we identified and ignored all outlier samples that could weaken the cluster generation. The remaining samples / genes were then subjected to consensus clustering in order to identify consistent subtypes in the dataset. All code can be found in “Generate_Clustering_and_Network.R” script



Identification of modulated/unmodulated genes To define proteins with variation between tumors we calculated a modified “quantile” standard deviation (sdQ) for each protein by ignoring the lowest and highest values of each protein. The distribution of sdQ was modeled as a mixture of gaussian distributions and we used an expectation maximization method (EM) to estimate the different mixture components using the package mixtools (https://CRAN.R-project.org/package=mixtools). The EM process converged in a 2 distribution solution which we assumed to represent the “Modulated” and “unModulated” proteins. Using this model we estimated the number of “Modulated” and “unModulated” proteins and selected a sdQ threshold which optimized the number of Modulated minus unModulated proteins. As the EM process inevitably produces slightly different thresholds every time it is executed, we performed 10 iterations and rounded the mean of the iterations to the nearest 0.5 in order to have a fully reproducible solution.
Tumor





---

three main classical subtypes are defined by expression of the oestrogen receptor (ER), progesterone receptor (PR; ERPR positive breast cancer) and the epidermal growth factor receptor ErbB2/Her2 (Her2 positive). The triple negative (TN) form (where none of the three markers is expressed) has an especially poor prognosis.

These subtypes matured into four accepted subtypes: Luminal A, Luminal B, Her2-enriched and basal-like breast cancer. While they do not perfectly reflect the clinical subtypes, most luminal tumours are ER/PR-positive, most Her2-enriched ones harbour the gene amplification, and most basal tumours are triple negative

Another attempt to bring the molecular classification of
breast cancer into the clinical practice has been to identify surrogate markers that would allow the identification of the molecular subtypes using the more familiar immunohistochemical approach. Accordingly, the combined evaluation of ER, PgR, HER2 and Ki67 immunoreactivity would approximate the molecular classification of luminal A, luminal B, HER2-enriched and basal-like breast cancers [25]. It should be noted, however, that the 

immunohistochemically and molecularly defined classes do not overlap completely

* Some basal-like breast cancers (according to the molecular classification) will not show the expected triple- negative (ER, PgR and HER2 negative) immunophenotype
* not all the immunohistochemically triple- negative breast cancer will be classified as basal-like by gene expression profiling [26].



the assignment of basal-like and ERBB2+ (HER2+) breast cancers to the poor-prognosis groups by first-generation gene signatures is determined mainly by high expres- sions of proliferation-related genes [



expression of genes involved in cell proliferation is the most heavily weighted compo- nent in calculating the recurrence score [4] and is the basis of genomic grading [5,6]



The breast is comprised of two main types of tissue: (a) glandular tissues that house
lobules (i.e., milk-producing glands) and ducts, and (b) stromal tissues, which provide the “supporting framework” of the breast (i.e., fatty and connective tissues)



The invasive form of breast cancer is seen to breach the duct and lobular wall, invading the breast's supporting tissues. Despite that, this does not imply its metastatic capability, as some breast cancer can be invasive but do not spread further to other organs or lymph nodes.

The other form of breast cancer is carcinoma in situ or non-invasive breast cancer, which can be subcategorized into two types: either ductal carcinoma in situ (DCIS) or lobular carcinoma in situ (LCIS). Generally, DCIS is considered as pre-invasive or non-invasive breast cancer and constitutes one in five new breast cancer cases reported. However, untreated DCIS may spread over time by invading into adjacent breast tissue and progress into invasive cancer. 



Perou and team studied the gene expression of different breast cancer types using cDNA
microarray profiling and sorted them into few intrinsic gene subtypes: basal-like breast cancer (BLBC) (estrogen receptor (ER) negative, progesterone-receptor (PR) negative, and human epidermal growth factor receptor-2 (HER2) negative), HER2-enriched, normal breast-like, luminal subtype A (ER+ and low grade) and luminal subtype B (ER+ and often high grade)[25– 27] and the poorly described triple-negative subtype or “claudin-low” subtype [25–27] (Figure 1). Each of these molecular subgroups has a distinctive prognosis and chemotherapy sensitivity.



the growth of breast cancer cells can be induced by hormones, estrogen,
and progesterone. Thus, selecting specific therapy for breast cancer varies depending on these hormone receptors' expression[28,29]. Generally, 70% of breast cancer patients express ER and receive hormonal therapy treatment[28,30,31]. Studies also indicated that those with luminal-like cancers seem to have better long-term survival compared with other cancers, but those with basal-like and HER-2-positive tumors show higher susceptibility and sensitivity towards chemotherapy[26,32]. In contrast, patients diagnosed with triple-negative cancer have exceptionally shorter survival in first metastatic occurrence than other breast cancer types[33



<u>ER-positive luminal breast cancer</u>

* the most common type of invasive breast cancer

* further classified into subtypes with different outcomes (i.e., luminal A and luminal B) using gene expression profiling with microarrays
  * Luminal A:
    *  high expression of ER-activated genes,
    * low- grade histological appearance, and slower growth, 
    *  resulting in a better prognosis than lumB
    *  associated with the expression of luminal epithelial cytokeratins (CK) 8 and 18, other genes such as hepatocyte nuclear factor 3 alpha (FOXA1), B cell lymphoma 2 (BCL2), erbB3, and erbB4[36]. 
  *  luminal B tumors 
    *  high proliferate rate and tumor aggressiveness than lumA

* The main difference between these two luminal subtypes is the greater level of proliferation- related genes expressions such as avian myeloblastosis viral oncogene homolog (v-MYB), gamma glutamyl hydrolase (GGH), lysosome-associated transmembrane protein 4-beta (LAPTMB4), nuclease sensitive element binding protein 1 (NSEP1) and cyclin E1 (CCNE1) in luminal-B tumors
* luminal tumors form a continuum, thus the determination of these tumors into two subtypes based on proliferation may be arbitrary



<u>The HER2 is overexpressed in 15–30% of invasive breast cancer,</u> 

* mainly due to the overexpression of HER2/HER2 signaling-related genes and genes located in the HER2 amplicon on chromosome 17q12

* can proliferate faster than luminal cancers, they are often highly responsive to targeted therapies aimed at the HER2 protein, resulting in remarkably improved outcomes[43]. 



<u>The normal cell-like subgroup,</u>

* 5–10% of breast cancer patients[44]. 
*  featured by similar gene expression to normal breast epithelium. 
* proliferation rate is often low and responds to adjuvant chemotherapy[44,45].



 <u>BLBC</u> 

* an aggressive subtype that mainly occurs in young women 

* displays exceptionally high metastasis rates to the brain and lung

* characterized by high histological grade, poor tubule formation, central necrotic zones, pushing borders, high mitotic indices, and proliferation rates[47]. 

* does not benefit from anti- estrogen hormonal treatment or trastuzumab as estrogen receptor (ER), progesterone receptor (PR), and HER2 were not expressed in this subtype. 





<u>Luminal Cancers</u>

* immunophenotypic pattern similar to the epithelial, luminal component of the milk ducts of the normal mammary gland, mainly expressing low molecular weight luminal cytokeratins (CK7, CK8, CK18, etc.), E receptors and associated genes (LIV1 and cyclin D1)19

* low association with proliferative genes. Three groups are distinguished from the IHC point of view: Luminal A, Luminal B and HER2 luminal. 
  * <u>Luminal A subtype</u>
    * (E+, Pr+, HER2- and Ki-67<14%)
    * the most common and least aggressive subtype, 
    *  very good prognosis, 
    * very low expression of proliferative genes. 
    * E-positive, Pr-positive (at least 20%) and HER2-negative receptor tumors, with a low Ki-67, less than 14%
    * That is how tubular carcinoma and infiltrating cancers (both ductal and lobular) of histological grade I and II (Figure 1) are. 
    * By expressing E receptors, these carcinomas are
      susceptible to being treated with hormone therapy (tamoxifen or aromatase inhibitors), in addition to surgical or radio/CT treatment that may be required. 
  * <u>Luminal B subtype (E+, Pr+/-, HER2- and Ki-67: 14-30%)</u> 
    * E-receptor positive, although these are usually expressed in less quantity, they can be Pr-positive or not, HER2-negative with an intermediate proliferation index, greater than 14%, but less than 25-30% and they are generally of intermediate/high histological grade. 
    * most BRCA2 cancers belong to this group. 
    * can benefit from hormone therapy along with chemotherapy (CT)20
    * The elevation of Ki-67 makes them grow faster than Luminal A and so have a worse prognosis.



<u>HER2 positive cancers</u> 

* HER2 proto-oncogene (or cerb-B-2) is found on chromosome 17 and is over-expressed in many epithelial tumors. It encodes a protein in the membrane of malignant cells with tyrosine kinase activi- ty. 

* 15-20% of breast carcinomas, tumor cells have additional copies of the HER2 gene and are often associated with alterations in other genes such as TOP2A, GATA4, angiogenesis genes, and proteolysis19. 

* The Ki-67 is always high. 

* usually of high histological grade and have a high proportion of mutations (40 to 80%) in p53 (a gene capable of detecting and repairing damaged DNA and causing cell death, its mutation increases the probability of suffering from cancer). 

* it is a more aggressive, fast growing subtype. 

* could be treated with specific drugs, targeting the HER2/neu protein: anti-HER2 monoclonal antibody (Trastuzumab or Herceptin, Pertuzumab) in addition to surgery and chemotherapy (CT) treatment if necessary.

* 2 subgroups
  * 1- HER2 Luminal (E+, Pr+, HER2 and Ki-67: 15-30%). 	
    * intermediate proliferation index, a Ki-67 of between 15-30%, and being HER2-positive, 
    *  also expresses hormonal receptors whether being of lower degree, 
    * therefore, apart from Herceptin, the hormone therapy is an additional therapeutic alternative. 
    * We could see it as a Luminal B HER2-positive3. It is also called a triple-positive tumor. 
  * 2- enriched HER2 (HER2+, E-, Pr-, Ki-67> 30%) 
    * very aggressive tumors with high Ki-6721 that do not have E or Pr receptors, 
    * therefore, they do not respond to hormonal therapy. 
    * Although about half respond to targeted therapy (Herceptin) and show a better response to CT, the prognosis is poor.



<u>Triple-negative receptor cancers (E-, Pr-, and HER2-)</u>

* The word triple-negative is a nomenclature based on IHC and refers to the phenotype of a group of tumors that are E, Pr, and HER2-negative, 

* 15% of breast cancers. 

* Using the gene expression profile, it was found that it is a very heterogeneous group.

* has been classified into several additional subgroups that include Basal-like subtypes (BL1 and BL2) (50-70%), Claudin- low, mesenchymal (MES) (20-30%), luminal androgen receptor (LAR), immunomodulator (IM) among others

* why basal-like ?  unlike Luminales, have gene expression patterns similar to the deepest, myoepithelial, basal component of the normal mammary gland, expressing high molecular weight cytokeratins (CK 5/6, CK 17) and are characterized by the absence of E receptor and HER2 gene expression. With additional IHC tests (which are not used in daily practice), such as CK 5/6 and EGFR, they can also be divided into two major subgroups: if the tumor expresses both, it is Triple-negative receptor phenoty- pe, but it is called Basal-like type (70-80%). If it does not express these two markers, it will be Non-Basal- like, actually a quintuple negative receptor. Although it is true that most triple-negative receptor tumors fall within the Basal-like spectrum (BL1 and BL2), these two terms are not synonymous and there is at least 30% disagreement between both classifications (IHC vs molecular)



<u>TN cancers</u> 

* common in young women 

* would have a worse prognosis in African American women

* usually have mutations in the p53 oncosuppressor gene

* poorly differentiated cancers, histological grade III with a particularly high mitotic index, usually with lymphocytic infiltration, areas of tumor necrosis, central fibrosis and circumscribed contours, without a stromal reaction

* aggressive behavior, with a high rate of local recurrence and early metastases, despite their high sensitivity to chemotherapy (CT) 

* 85% of BRCA1 tumors belong to this group

* By immunophenotype (E-, Pr- and HER2-), several low-grade carcinomas are also considered here, such as atypical medullary carcinoma, cystic adenoid, apo- crine, metaplastic variant squamous, adenosquamous, fibromatosis-like, etc. that have a better prognosis.

* In the rest, total survival is low, because endocrine therapies and Herceptin are ineffective in this group of tumors. At the moment only chemotherapy (CT)25 is offered. In patients with the BRCA1 mutation, certain PARP122 inhibitors may be effective.





<u>1- Luminal A,</u> 

* 30-70% of all breast cancers. 

* originates in the inner (luminal) cells lining the mammary ducts. 

* hormone- dependent and tends to be estrogen receptor (ER)- and/or progesterone (PR)-positive, HER2-negative 

* displays low proliferative activity based on Ki67 expression. 

<u>2- Luminal B,</u> 

* 10-20% of all breast cancers 

* characterized by ER/ PR-positivity, high Ki67, and/or HER2-positivity. 

1- and 2- patients have the best prognosis, and if they do progress to metastasis is usually to the bone



<u>3- HER2–enriched,</u>

*  ER/PR negativity and positive HER2,

* 5-15% of all breast cancers. 

<u>4-  triple negative/ basal-like (ER-/PR-/HER2-)</u>

* 15-20% of all breast cancers. 

* highest incidence of metastasis and show a high rate of metastasis to the lungs.

* often more aggressive and have a poorer prognosis (at least within the first five years after diagnosis) compared to the ER- positive subtypes (Luminal A and Luminal B tumors)

3- and 4-  are basal- like tumours due to the tumor cells expressing similar features to those of the outer (basal) cells surrounding the mammary ducts. 

 



Here we present an integrated genomic strategy based on the use of gene expression signatures of oncogenic pathway activity (n = 52) as a framework to analyze DNA copy number alterations in combination with data from a genome-wide RNA-mediated interference screen. We identify specific DNA amplifications and essential genes within these amplicons representing key genetic drivers, including known and new regulators of oncogenesis. The genes identified include eight that are essential for cell proliferation (FGD5, METTL6, CPT1A, DTX3, MRPS23, EIF2S2, EIF6 and SLC2A10) and are uniquely amplified in patients with highly proliferative luminal breast tumors, a



we used a series of experimentally derived gene expression signatures that are capable of measuring oncogene or tumor suppressor pathway activ- ity, aspects of the tumor microenvironment and other tumor char- acteristics, including proliferation rate, as a framework by which to integrate multiple forms of genomic data

By further analyzing functional data from a genome-wide RNA-mediated interference (RNAi) screen9, we identified genes that are essential for cell viability in a pathway-dependent and, in some cases, subtype-dependent manner.





. Several studies illustrated that the presence of tumor infiltrating lympho- cytes was often associated with better prognosis several cancer types, including breast cancer12–21. Thus, inclusion of immune signatures in the molecular subtyping may provide additional information beyond routine prognosis in breast cancer. However, until now, no attempt has been made to use these immune signatures to stratify breast cancer



by using the gene expression profiles of immune-related genes with favorable prognosis, the
k-means clustering was applied on the breast cancer samples to establish a robust molecular classification. Then, the associations between the molecular clusters and prognosis, clinicopathological factors, immune-related fea- tures and tumor infiltrating levels were assessed. The





<u>Luminal cancers</u> 

* express Estrogen (ER) and/or Progesterone (PR) Receptors but lack Human Epidermal Growth Factor Receptor-2 (HER2); 

* by far the most common (>70% of cases) 

* have the best prognosis among the three major subtypes [2]. 

* treated with endocrine therapies like fulves- trant or tamoxifen that target ER or aromatase inhibitors that suppress estrogen production. 

* approximately one-third of Luminal cancers acquire hormone resistance, enhan- cing their aggressiveness and leading to local recurrence or distant metastases

<u>HER2- overexpressing cancers</u>; 

* basal-like cancers that lack ER, PR and HER2. 





<u>Triple-negative breast cancer</u> 

* heterogeneous disease characterized by a lack of hormonal receptors and HER2 overexpression. 

* the only breast cancer subgroup that does not benefit from targeted therapies, and its prognosis is poor. 

* defined by a lack of ER and PR expression and a lack of HER2 overexpression. 

* Most triple-negative tumors are included in the so-called basal-like molecular subgroup3, although both categories have up to 30% discordance4



In the present study, we applied probabilistic graphical models to a previously published TNBC cohort5. This
technique allows exploring the molecular information from a functional perspective.



Probabilistic graphical model analysis. A probabilistic graphical model compatible with a high-dimensionality approach to associate gene expression profiles, including the most variable 2000 genes, was performed as previously described18. Briefly, the resulting network, in which each node represents an individual gene, was split into several branches to identify functional structures within the network. Then, we used gene ontology analyses to investigate which function or functions were overrepresented in each branch, using the functional annotation chart tool provided by DAVID 6.8 beta19. We used “homo sapiens” as a background list and selected only GOTERM-DIRECT gene ontology categories and Biocarta and KEGG pathways. Functional nodes were composed of nodes presenting a gene ontology enriched category. To measure the functional activity of each functional node, the mean expression of all the genes included in one branch related to a concrete function was calculated. Differences in functional node activity were assessed by class comparison analyses. Finally, metanodes were defined as groups of related functional nodes using nonsupervised hierarchical clustering analyses.
Sparse k-means classification. Sparse k-means was used to establish the optimal number of tumor groups. This method uses the genes included in each node and metanode, as previously described20. Briefly, classification consistency was tested using random forest. An analysis using the consensus clustering algorithm21 as applied to the data containing the variables that were selected by the sparse K-means method22 has provided an optimum classification into two subtypes in previous studies20. In order to transfer the newly defined classification from the main dataset to other datasets, we constructed centroids for each defined subgroup, using genes included in various metanodes.



Hormone receptor-positive (HR-positive) breast cancer: which is
defined as ER-positive, PR-positive and HER2-negative (HR+, HER2-). The majority of breast cancer cases (up to 70%) fall in this category.19 These tumors have the best prognosis and, in general, respond well to hormone therapy; however, the response rates vary according to a gradient of ER expression.11 HER2-positive breast cancer: which is defined as ER-negative,
PR-negative, and HER2-positive (HR-, HER2+). Approximately 20%– 25% of breast cancers worldwide exhibit amplification of the ERBB2 gene and overexpression of HER2.20 An increase in membrane pres- ence of HER2 is related to high tumor grade, high mitotic count and positive lymph nodes which are poor prognostic factors.21 These tumors respond poorly to chemotherapy11 and endocrine therapy,22 but they are the candidates for anti-HER2 treatment.20 Triple-Negative Breast Cancer (TNBC): which is defined as
ER-negative, PR-negative and HER2-negative (HR-, HER2-). Typically have the worst prognosis and are treated with systemic chemother- apy to which they respond better than other subtypes.2,4 Due to the lack of recognized molecular targets for therapy, TNBC is a subtype of interest for clinical trials with novel treatment approaches. ER-positive, PR-positive and HER2-positive (HR+, HER2+) breast
cancer: which is also referred to as Triple Positive. HR+, HER2+ tumors are typically associated with an intermediate prognosis and can be treated with a combination of hormone therapy, chemotherapy and/or anti-HER2 treatment





Luminal A – Tumors identified as Luminal A subtype are low grade
and have a high expression of luminal epithelial genes, ER, PR, and low expression of HER2.2,4 These tumors are associated with somatic muta- tions in PIK3CA, GATA3, and MAP3K1 genes and often exhibit cyclin D1 overexpression.4 Luminal A tumors generally have a good prognosis and respond well to hormone/endocrine treatment.5 For these reasons, Luminal A patients are considered for chemotherapy de-escalation. Luminal B – Tumors identified as Luminal B are high grade and
have low ER and PR expression compared to Luminal A tumors.4 Luminal B tumors show frequent mutations in the PIK3CA and TP53 genes as well as deviations in the MAPK and retinoblastoma path- ways.4 A subset of Luminal B tumors show overexpression of HER2 and these HER2-positive Luminal B tumors are associated with a worse overall outcome than Luminal B tumors that are HER2-negative.4,28 Compared to Luminal A, Luminal B tumors have
worse





worse prognosis, but in comparison to all the subtypes, its prognosis is intermediate.2,5 As this subtype is more proliferative, treatments that involve both chemotherapy and hormone treatment will be bene- ficial. If the tumor expresses HER2, an anti-HER2 treatment will be an effective measure as well.2,4 HER2-Enriched – Tumors defined as HER2-Enriched are not
determined solely by HR status or HER2 gene expression. HER2-Enriched tumors are driven by the EGFR/HER2 signaling pathway,29 comprised of four major proteins: EGFR, HER2, HER3, and HER4.30 Increased dimerization of HER2 and EGFR results in prolifer- ative activity.30 Expression of HER2-regulated genes serves as useful predictive markers to determine whether a treatment targeting the HER2 pathway will be beneficial.31 HER2-Enriched tumors generally exhibit high expression of HER2 and are negative for luminal epithelial genes (non-luminal). These tumors are highly proliferative and are associated with poor prognosis. Possible treatment plans involve HER2 blockade therapies which include the use of anti-HER2 anti- bodies or small molecule inhibitors.2,12



Basal-Like – Tumors identified as Basal-Like have high expression
of genes typically found in basal cells (i.e., KRT5, KRT17, LAMC1). Basal tumors display specific epidemiological, phenotypic, and molec- ular features, including high grade, proliferation, areas of necrosis, and distinct genetic alterations. There is also increased activation of the WNT pathway and higher expression of EGFR.4 These tumors do not express luminal genes, ER/PR, or HER2 and frequently show muta- tions in the TP53 and BRCA1 genes. Basal-Like tumors are found to be aggressive and associated with a poor prognosis and high risk of relapse.2 Chemotherapy is the recommended treatment for tumors classified as Basal-Like.



The intrinsic subtypes were initially regarded as compatible with clas- sical IHC-based subtype definitions. The 2013 St Gallen meeting advocated for useful, surrogate definitions, created based on immuno- histochemical markers, to describe intrinsic subtypes without a requirement for molecular diagnostics.33 The HER2-positive disease was assumed representative of HER2-Enriched subtype, and TNBC interchangeable with Basal-like. The HR-positive disease was consid- ered equivalent of Luminal, and the Ki-67 marker was introduced to separate Luminal A from Luminal B tumors, as it became evident that this distinction had therapeutic implications. Ki-67 is a nuclear prolif- eration marker, and high levels of this marker correspond with poor outcome.5 The Ki-67 index was established to serve as a prognostic marker of overall survival and a predictive marker of chemotherapy response and recurrence.5,34 The Luminal A subtype was considered ER+ and/or PR+, HER2-, and Ki-67 low, while Luminal B was consid- ered ER+ and/or PR+, HER2- (or HER2+) and Ki-67 high. The surro- gate definitions and clinical characteristics of the four intrinsic subtypes are summarized in Table 1. In the last decade, GEP became increasingly employed by clini-
cians as supplementary prognostic indicators in clinical trials. Numer- ous studies that compared the results from intrinsic and clinical subtyping performed on the same set of tumors manifested consider- able discordance in the two classification approaches' results. The dis- crepancy has been acknowledged by experts who, however, do not comment on the superiority of any approach and continue rec- ommending IHC to define breast cancer subcategories based on ER, PR, HER2, and Ki67 staining. To date, IHC has near-exclusive use in contemporary practice.



Human ERBB2 is a proto-oncogene located on chromosome-17, which encodes the tyrosine kinase receptor HER2, involved in activa- tion of proliferation pathways and affects differentiation, invasive- ness, and survival.5,14 The amplification of ERBB2 was distinguished for its remarkable prognostic value.1







In this article, we combined quantitative proteomics with
miRNA expression analyses in a series of breast cancer samples with appropriate clinical information. We demonstrate that it is possible to perform differential protein expression analyses by LC/MS-MS on tens of FFPE tumors. Protein patterns confirm that estrogen receptor–positive (ERþ) and triple-negative breast cancer (TNBC) subtypes are distinct biologic entities (9). Also, we were able to profile miRNA expression in the same series of samples, which led us to conduct a probabilistic graphical models analysis to integrate miRNA and protein expression data. This



Gatza et al. 2014

we present an integrated genomic strategy based on the use of gene expression signatures of oncogenic pathway activity (n = 52) as a framework to analyze DNA copy number alterations in combination with data from a genome-wide RNA-mediated interference screen. 

We identify specific DNA amplifications and essential genes within these amplicons representing key genetic drivers, including known and new regulators of oncogenesis. The genes identified include eight that are essential for cell proliferation (FGD5, METTL6, CPT1A, DTX3, MRPS23, EIF2S2, EIF6 and SLC2A10) and are uniquely amplified in patients with highly proliferative luminal breast tumors, a

we used a series of experimentally derived gene expression signatures that are capable of measuring oncogene or tumor suppressor pathway activ- ity, aspects of the tumor microenvironment and other tumor char- acteristics, including proliferation rate, as a framework by which to integrate multiple forms of genomic data. O

By further analyzing functional data from a genome-wide RNA-mediated interference (RNAi) screen9, we identified genes that are essential for cell viability in a pathway-dependent and, in some cases, subtype-dependent manner. Our results identify a small number of DNA amplifications as potential drivers of proliferation in poor-outcome luminal breast cancers,

The patterns of pathway activity recapitulated known character-
istics of each subtype, including dysregulation of pathways that can be linked to female hormone receptors, oncogenes and/or tumor suppressor mutation status (Fig. 1). For example, basal-like tumors

which represent ~80% of triple-negative breast cancers, are character- ized by low hormone receptor signaling, mutant p53 signaling and high expression of proliferation pathway activity (Fig. 1). Likewise, HER2-enriched (HER2E) tumors show high expression of the HER2 (ref. 11) and HER2 amplicon (HER2-AMP)12 signatures, whereas luminal A (LumA) tumors show high hormone receptor signaling and wild-type p53 signaling. Highly proliferative LumB tumors, which also show some hormone receptor signaling, are distinguished from less proliferative LumA samples by increased proliferation-associ- ated pathways. Thus, these data robustly recapitulate many previously published pathway and subtype associations.

we sought to identify amplified genes and/or CNAs associated with our previously published 11-gene PAM50 proliferation signature with the hope that these might represent targetable drivers of oncogenesis

we found that basal-like, LumB and HER2E tumors had the highest proliferation levels (Fig.

To identify genes that are specifically amplified in highly proliferative luminal breast cancer, we repeated these analyses using the luminal tumor sub- set (Figs. 3g,h and Supplementary Table 58). Analyzing both popu- lations of tumors identified three classes of proliferation-associated regions (q < 0.05): (i) CNAs associated irrespective of subtype, (ii) CNAs altered in basal-like tumors, and (iii) CNAs altered in highly proliferative luminal tumors. These results allowed us to focus our anal- yses on those genes within regions that are uniquely altered in highly proliferative luminal tumors by censoring proliferation-associated genes that are altered in basal-like breast cancer (e.g., TP53 or INPP4B loss) or that are altered irrespective of molecular subtype (e.g., RB1 loss or MYC amplification). These

an RNAi proliferation screen in which a genome-wide shRNA library (~16,000 genes) had been used to identify essential genes (Fig.

We applied the 52 gene expression signatures to a panel (GSE12777)41 of breast cancer cell lines (Supplementary Fig. 5 and Supplementary Table 59), 27 of which had mRNA expres- sion data and were also part of an RNAi proliferation screen in which a genome-wide shRNA library (~16,000 genes) had been used to identify essential genes (Fig. 4a)9. For each signature, we used a negative Spearman rank correlation to identify pathway-specific essential genes (Fig. 4b and Supplementary Table 60) by comparing the path- way score against the normalized shRNA score across the panel of 27 cell lines. These analyses identified inverse relationships between the abundance of shRNAs targeting key regulatory genes and pathway scores. For

These results confirm that this approach is able to identify essential genes that are known to be functionally associated with pathway activity, thereby suggesting that these data can serve as a biological filter to distinguish pathway-specific essential from nonessential genes



We next sought to distinguish between essential and nonessential genes within regions amplified specifically in highly proliferative luminal tumors. For each subset of tumors, we identified genes in amplified regions that were positively correlated with proliferation and showed an increased amplification frequency 



In this study we utilized gene expression signatures of signaling
pathways to identify patterns that can distinguish the known subtypes of breast cancer. These signatures were developed largely from con- trolled manipulations of the relevant pathways in vitro and are thus based on experimental evidence for pathway activation as opposed to extrapolations of pathway activity achieved from analyses of anno- tated gene lists. Therefore, the use of an experimentally derived path- way signature, as opposed to an analysis of a single genomic alteration, provides a measure of pathway activity irrespective of how the path- way may have been activated. 



a given pathway can be active in a subset of tumors as a result of either an activating altera- tion (i.e., E2F1 or E2F3 amplification) or an independent event that inactivates a negative regulator of the pathway (i.e., RB1 loss and/or mutation), which nevertheless achieves the same end result (i.e., DNA replication and cell proliferation);



Proliferation is one of the most powerful prognostic features in
breast cancers, especially for ER+ cancers38,39. Because proliferation is so important, we used a gene expression signature of proliferation as a means to integrate the DNA copy number data, along with data from a genome-wide RNAi screen of luminal breast cancer cell lines, to identify luminal-specific genetic drivers of proliferation. We iden- tified 12 genes that were amplified uniquely in highly proliferative luminal tumors in the TCGA data set, have a correlation between mRNA expression and DNA copy number and have been shown to be essential for luminal breast cancer cell line viability; we validated 8 of these genes using the independent METABRIC data set. Whereas FGD5, METTL6, DTX3 and MRPS23 amplification was prognostic in luminal tumors, these and many of the other identified genes have been reported previously to regulate tumorigenic characteristics, albeit not necessarily in human breast cancer. For example, FGD5 has been shown to regulate the proangiogenic function of VEGF43, potentially leading to increased proliferation. DTX3 purportedly promotes Notch signaling44,45, whereas EIF6 is a Notch-dependent regulator of cell invasion and migration46, and its inhibition restricts lymphomagenesis and tumor progression47. MRPS23 expression is associated with proliferation, oxidative phosphorylation, invasive- ness and tumor size in uterine cervical cancer48. METTL6 has been reported to contribute to cytotoxic chemotherapy sensitivity in lung cancers49

<u>Gene expression signatures.</u>A panel of 52 previously published gene expression signatures was used to examine patterns of pathway activity and/or microenvironmental states ([Supplementary Table 1](https://www.nature.com/articles/ng.3073#MOESM11)). To implement each signature, the methods detailed in the original  studies were followed as closely as possible. Of these 52 signatures, 22 signatures[10](https://www.nature.com/articles/ng.3073#ref-CR10),[11](https://www.nature.com/articles/ng.3073#ref-CR11),[32](https://www.nature.com/articles/ng.3073#ref-CR32) were originally developed using a Bayesian binary regression strategy  and are comprised of Affymetrix probe sets with positive and negative  regression weights. These signatures were translated to a form that  could be applied to non-Affymetrix expression data. For each signature,  we excluded those probe sets with a negative correlation coefficient.  The remaining probe sets with a positive coefficient were then  translated to the gene level, and replicate genes were merged. To apply a given signature to a new data set, the expression data were filtered to contain only those genes that met the previous criteria, and the mean  expression value was calculated using all genes within a given signature that were present in more than 80% of samples. The list of genes in  each modified signature is shown in [Supplementary Table 2](https://www.nature.com/articles/ng.3073#MOESM11), and the scores for the TCGA data set ([Supplementary Table 3](https://www.nature.com/articles/ng.3073#MOESM11)) and cell line data set are also provided ([Supplementary Table 59](https://www.nature.com/articles/ng.3073#MOESM11)).



Because highly proliferative luminal tumors have a poor prognosis and poor responses to existing therapies38,39, we sought to identify amplified genes and/or CNAs associated with our previously published 11-gene PAM50 proliferation signature with the hope that these might represent targetable drivers of oncogenesis.



To identify those genes that are altered specifically in highly prolif-
erative luminal tumors while excluding those that are associated with proliferation irrespective of subtype, we performed analyses on two subsets of samples: all tumors and all non–basal like tumors (hence- forth called luminal tumors). 

Examining the TCGA breast cancer data set using the PAM50 pro-
liferation signature31, we found that basal-like, LumB and HER2E tumors had the highest proliferation levels (Fig.



Using the PAM50 proliferation signature, we examined the frequency of CNA gains and losses in highly proliferative (top quartile) tumors relative to less prolifera- tive samples irrespective of subtype using the statistical strategies discussed previously 

To identify genes that are specifically amplified in highly proliferative luminal breast cancer, we repeated these analyses using the luminal tumor sub- set (Figs. 3g,h and Supplementary Table 58). Analyzing both popu- lations of tumors identified three classes of proliferation-associated regions (q < 0.05): (i) CNAs associated irrespective of subtype, (ii) CNAs altered in basal-like tumors, and (iii) CNAs altered in highly proliferative luminal tumors. These results allowed us to focus our anal- yses on those genes within regions that are uniquely altered in highly proliferative luminal tumors by censoring proliferation-associated genes that are altered in basal-like breast cancer (e.g., TP53 or INPP4B loss) or that are altered irrespective of molecular subtype (e.g., RB1 loss or MYC amplification). 

Statistical analyses of signature scores. To quantify differences in patterns of signature scores across subtypes, ANOVA followed by Tukey’s post-test for pairwise comparisons was used (as shown in Fig. 1b). To investigate the level of concordance between each of the 52 signatures, the pathway scores calculated for each sample in the TCGA data set (Supplementary Table 3) were analyzed. The R values calculated by Pearson correlation are reported in Supplementary Figure 2 and Supplementary Table 5



Analysis of genome-wide RNAi proliferation data. To identify genes that are required for cell viability in a signature-dependent manner, data from a previously published genome-wide RNAi screen carried out on a panel of breast cancer cell lines were analyzed9. The Gene Active Ranking Profile (GARP)-normalized data were obtained from the COLT database and fil- tered to include only those 27 cell lines for which gene expression data (GSE12777) were also available (acquired February 2013). To identify genes essential for pathway-dependent cell proliferation, a negative Spearman cor- relation was performed comparing predicted pathway activity and GARP score for each sample. A threshold of P < 0.05 was considered significant for all analyses



Breast cancer is the most frequent malignancy in women worldwide and is curable in ~70–80% of patients with early- stage, non- metastatic disease.

Breast cancer is the most frequent malignancy in women and is a heterogeneous disease on the molecular level

Early breast cancer — that is, cancer that is contained in the breast or that has only spread to the axillary lymph nodes — is considered curable. Improvements in multimodal ther- apy have led to increasing chances for cure in ~70–80% of patients. By contrast, advanced (metastatic) disease is not considered curable using currently available thera- peutic options. However,

the histo- logical and molecular characteristics of breast cancer largely influence treatment decisions. The molecular
alterations that drive breast carcinogenesis are many, and several classifications have been developed to group tumours accordingly

The intrinsic classification of Perou and Sorlie1
, reported in 2000, distinguished
four subtypes of breast cancer: luminal A and luminal B (expressing the oestrogen receptor (ER)), basal- like and human epidermal growth factor receptor 2 (HER2)- enriched (without ER expression). This classification shifted clinical management of breast cancer from being based on tumour burden to biology- centred approaches



clinical practice typically uses a surrogate classification of five subtypes on the basis of histological and molecular characteristics (Fig. 1). Tumours express- ing ER and/or progesterone receptor (PR) are consid- ered hormone receptor- positive breast cancers, whereas tumours that do not express ER, PR or HER2 are triple- negative breast cancer (TNBC).



Approximately 10% of breast cancers are inherited and associated with a family history2



. The histological subtypes described here (top right) are the most frequent subtypes of breast cancer ; ductal carcinoma (now referred to as ‘no special type’ (NST)) and lobular carcinoma are the invasive lesions; their preinvasive counterparts are ductal carcinoma in situ and lobular carcinoma in situ (or lobular neoplasia), respectively



The intrinsic subtypes of Perou and Sorlie1
are based on a 50-gene expression signature (PAM50)321 .

The surrogate
intrinsic subtypes are typically used clinically and are based on histology and immunohistochemistry expression of key proteins: oestrogen receptor (ER), progesterone receptor (PR), human epidermal growth factor receptor 2 (HER2) and the proliferation marker Ki67. Tumours expressing ER and/or PR are termed ‘hormone receptor- positive’; tumours not expressing ER , PR and HER2 are called ‘triple- negative’. The



At the molecular level, there is evidence showing that breast cancer evolves along two divergent molecular pathways of progression, mainly related to ER expression, and tumour grade and proliferation (described in the intrinsic classification)

The first pathway —the low- grade-like pathway
— is characterized by gain of 1q, loss 16q, infrequent amplification of 17q12 and a gene expression signature
(GES)(GES) with a majority of genes associated with the ER phenotype, diploid or near diploid karyotypes and low tumour grade. The luminal A group and to some extent the luminal B group fall into this pathway.



The second pathway — the high- grade-like pathway — is charac- terized by loss of 13q, gain of chromosomal region 11q13, amplification of 17q12 (containing ERBB2, encoding HER2) and an expression signature of genes involved in the cell cycle and cellular proliferation61
.
Tumours composed of intermediate to high grade, including HER2-positive tumours and TNBC, fall into this pathway6



The majority of the mutations affecting 100 putative breast cancer drivers are extremely rare64
, therefore, most breast cancers are caused by
multiple, low- penetrant mutations that act cumula- tively. 



Luminal A tumours have a high prevalence of PIK3CA mutations (49%), whereas a high prevalence of TP53 mutations is a hallmark of basal- like tumours (84%). For



Epigenetic alterations are involved in breast carcino- genesis and progression. In breast cancer, genes can be either globally hypomethylated (leading to gene acti- vation, upregulation of oncogenes and chromosomal instability) or, less frequently, focally (locus- specific) hypermethylated (leading to gene repression and genetic instability due to the silencing of DNA repair genes). 



Hormone receptors. The major risk factors for spo- radic breast cancer are linked to hormone exposure. Oestrogen is clearly a promoter of breast cancer, through its binding of the ER located in the nucleus (encoded by ESR1), which is a ligand- activated transcription factor.



drugs blocking
the effects of oestrogen on the mammary gland, such as tamoxifen, or drugs that block the production of oestrogen, such as aromatase inhibitors, have major roles in the treatment of hormone- sensitive breast cancer. As oestrogen interacts with bone, aromatase inhibi- tors can also cause osteoporosis (as menopause does). By contrast, tamoxifen has oestrogen- like effects on the bone, thereby preventing osteoporosis7





HER2. ERBB2 is amplified in 13–15% of breast cancers, causing an activation of the HER2 pathway.



According to the latest edition of the WHO classification, breast carcinomas are divided into 19 different major subtypes, including invasive carcinomas of no special type (70–75%; also known as not otherwise specified (formerly ductal car- cinoma)), lobular carcinomas (10–14%) and the other carcinomas of special type (including 17 different rare histotypes and their subclassifiers)107 (Fig. 5). Breast can- cer of ‘no special type’ is a carcinoma that does not fit into a specific histotype. Some of the special types (such as tubular, cribriform and mucinous) — if at least 90% pure (that is, no mixed histology or <10% of another subtype) — have a very good prognosis. On the other hand, some other special types (such as pleiomorphic lobular carcinoma, high- grade metaplastic carcinoma and micropapillary carcinoma) are associated with the poorest clinical outcome. Another special case is inflam- matory breast cancer, a rare and aggressive form charac- terized by malignant cells blocking the lymph vessels in the skin of the breast (Bo



biomarkers validated for therapy decision making

* ER (IHC)
  * characterization of IHC luminal group
  * poor prognosis if negative 
* PR (IHC)
  * if negative, classified as IHC luminal B
  * strong poor prognosis if negative
* HER2 (IHC)
  * essential to characterize HER2-enriched (ER-negative) and luminal B, HER2+
  * prognostic marker
* Ki67 (IHC)
  * prognostic value in ER+, HER2- 
  * Part of the IHC definition of luminal tumours whereby when Ki67 is low , luminal A tumour likely and when Ki67 high, luminal B tumour likel

* intrinsic subtypes (gene expression profiles)
* first-generation signatures (gene expression profiles)
* PIK3CA mutations
  * predictive marker for PI3KCA inhibitors in lumA and lumB metastatic cancers
* germline BRCA mutations
* PD-L1 (IHC)
  * predictive for immunotherapy





we identified 73 breast cancer pa- tients from the TCGA and CPTAC projects matched whole slide images, RNA-seq, and proteomic data.
