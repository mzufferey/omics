#### Common Model Inputs Used in CISNET Collaborative Breast Cancer Modeling - Mandelblatt et al. 2018



the Cancer Intervention and Surveillance Network (CISNET) breast can- cer models have collaborated to use a nationally representative core of common input parameters to represent key components of breast cancer control in each model. 



The common core of para- meters includes 

- population rates of births and deaths; 
- age- and cohort-specific temporal rates of breast cancer inci- dence in the absence of screening and treatment; 
- effects of risk factors on incidence trends; 
- dissemination of plain film and digital mammography; screening test performance characteristics; 
- stage or size distribution of screen-, inter- val-, and clinically- detected tumors by age;
- the joint distribution of ER/HER2 by age and stage; survival in the absence of screening and treatment by stage and molecular subtype; 
- age-, stage-, and molecular subtype-specific ther- apy; 
- dissemination and effectiveness of therapies over time; 
- and competing non-breast cancer mortality



A key feature of the Cancer Intervention and Surveillance Network (CISNET) collaborative modeling approach is the shared use of a common set of input values.



common input values presently used in the CISNET breast cancer models to estimate trends in US breast cancer incidence and mortality



There are presently six models: 

1. Model D (Dana–Farber); 
2. Model E (Erasmus), 
3. Model GE (Georgetown–Einstein), 
4. Model M (MD Anderson), 
5. Model S (Stanford), 
6. Model W (Wisconsin– Harvard).



Based on the goals of any given analysis, there are also common inputs available for age- and gender-specific utilities and costs for model health states. 

The models either used the common input parameters directly, or as a calibration target depending on individual model struc- tures (Table

common inputs are used with model-specific parameters related to unobservable aspects of breast can- cer history (e.g., tumor growth, proportions, and types of tumors that are non-progressive, sojourn time, lead- time, and how systemic therapy affects survival);



when considering data sources for common parameters, CISNET uses the hierarchy of evidence pro- moted by the US Preventive Services Task Force to select available data of the highest quality for a given para- meter and research question





#### Estimating Breast Cancer Survival by Molecular Subtype in the Absence of Screening and Adjuvant Treatment - Munoz and Plevritis 2018



a modeling approach to estimate histori- cal survival outcomes by estrogen receptor (ER) and human epidermal growth factor receptor 2 (HER2) status

Our approach leverages a simulation model of breast cancer outcomes and integrates data from two sources: the Surveillance Epidemiology and End Results (SEER) databases and the Breast Cancer Surveillance Consortium (BCSC).

- produce ER- and HER2-specific estimates of breast cancer survival in the absence of screening and adjuvant treatment 
- also estimate mean tumor volume doubling time (TVDT) and mean mammographic detection threshold by ER/HER2-status. 



a consortium of independent investiga- tors from the Cancer Intervention and Surveillance Modeling Network (CISNET) reported on the use of several simulation-based models to assess the relative contributions of screening mammography and adjuvant treatment on the reduction in breast cancer mortality for the overall US population. 



In that analysis, all CISNET models began by **recreating incidence and mortality trends in the absence of screening and adjuvant treat- ment interventions.** Then, these interventions were super- imposed based on their dissemination and efficacies across calendar years to assess the effect of the presence of interventions relative to their absence on outcomes.



**possessing molecular-subtype data in the absence of screening and treatment are necessary** to estimate the impact of these interventions by molecular subtype using similar simulation-based approaches

**Assessing molecular-subtype data in the absence of screening and treatment, however, poses a significant challenge**. Given the relatively novel nature of clinically relevant molecular markers, such as estrogen receptor (ER) and human epidermal growth factor receptor 2 (HER2), **historical surveillance data reporting them are rare.** Therefore, using existing methods to infer the impact of screening and treatment on breast cancer trends by ER and HER2 status is not straightforward.

We present a modeling approach to **estimate population-level breast cancer survival by ER and HER2 status, in the absence of screening and treatment**.

Our approach makes use of a previously developed **natu- ral history model of breast cancer**18-20 to **integrate data from 2 distinct sources**: 

1. SEER 
2. Breast Cancer Surveillance Consortium (BCSC). 

Data on women detected with breast cancer between 1996 and 2010 pro- vided by the BCSC includes ER/HER2, mode of detec- tion, and screening histories (note the BSCS data source in Figure 1). Our method not only produces **ER/HER2- specific breast cancer survival cases in the absence of screening and adjuvant**, we simultaneously **produce sev- eral other ER/HER2-specific estimates, including: the distribution of ER/HER2-subtypes by age in the absence of screening, ER/HER2-specific tumor volume doubling times, and tumor size-specific mammography threshold by ER/HER2**. 

This study makes use of a previously developed model (**Model S**, also referred to as BCOS for ‘‘Breast Cancer Outcomes Simulator’’) to **simulate the natural history of breast cancer in the average-risk US population, incor- porating the effect of screening and treatment.**

We **modify Model S to estimate underlying breast cancer pro- gression and survival by ER/HER2 status**.

- we **stratify the natural history of breast cancer by tumor grade**—low (grade I+II) v. high (grade III)—following an approach used in our previous work.26,27.
  - The inclu- sion of grade into our natural history model serves to **leverage this feature’s relationship with ER and HER2 status**. 
  - In addition, it expands the model’s ability to **capture a broader spectrum of tumor aggressiveness**

Model S is well suited for this work because **its underlying natural history model is easily adaptable.**

our methodology for <u>estimating ER/HER2- specific breast cancer survival in the absence of screening and treatment</u> consists of several steps.

1. First, we con- structed **ER-specific and ER/HER2-specific classifiers that can infer these molecular markers based on a patient’s mode of detection, screening history, and tumor features at detection**. 
2. Then, we simulated an **enhanced virtual SEER breast cancer registry that includes patient- level information not found in SEER** such as the mode of detection, screening history, tumor features, and sur- vival in the presence and absence of screening. 
3. We used the molecular classifier (from the first step) to **assign the molecular profile for each individual patient.** This proce- dure may be conceptualized as an **imputation of molecular-specific markers run across a virtual patient registry**.
   - this virtual, enhanced database allows us to evaluate population-level outcomes by ER and HER2 subtypes as if these were measured directly in the general population

List of Steps Involved in the Estimation Procedure, 

1. Train ER/HER2 Classifiers
2. Estimate BCOS’s Natural History Model Parameters SEER
3. Run BCOS to Construct an Enhanced, Virtual SEER Registry
4. Predict ER/HER2 Status for Women in the Virtual SEER Registry
5. Compute Desired Estimates by ER/HER2 Status



*BCSC Data Used to Determine Grade, Dependent on Mode of Detection*

In addition to the inputs required for our natural history model that have been described in previous work,18-20, 27 our estimation also relies on new inputs obtained from data provided by the BCSC and the SEER registry.

*BCSC Data Used to Build ER/HER2 Classifier* 

We also used BCSC data to build 2 classifiers: 

1. to infer ER status; 
   -  a 2-class prediction: (1) ER-positive v. (2) ER-negative; 
   - Alternating Decision Trees (**ADTree**): combines decision trees, voted decision trees, and voted decision stumps based on the concept of boosting, which produces accurate predictions by combining a series of ‘‘weak’’ learners together
2. to infer combined ER/ HER2 status. 
   -  a prediction across 4 different classes: (1) ER-positive, HER2-positive, (2) ER-positive, HER2-negative, (3) ER-negative, HER2- positive and (4) ER-negative, HER2-negative.
   -  **LADTree** (multi-class counterpart of ADTree): pro- duces a multi-class, prediction-alternating, decision tree and uses the LogitBoost strategy, which performs additive logistic regression

To predict the molecular markers based on tumor features, the tree is traversed, adding a differ- ent score at each decision stump. The final prediction is made by choosing the molecular marker with the high- est score after reaching the lowest leaf. 

2 separate classifiers corresponding to low- and high-grade tumors,

*SEER Survival in Pre-screening Period Defines the Baseline Survival Curves*

To ultimately determine ER/HER2-specific survival for patients in the absence of screening mammography and adjuvant treatment, we leveraged SEER survival curves for cases detected between 1975 and 1981, when the use of screening and adjuvant treatment were not widespread in the general population.

*Molecular-specific Parameter Estimation*
The premise behind our estimation is 

- to generate the **individual-level characteristics** of each breast cancer case during its pre-clinical course and at the time of clinical detection in the absence of screening,
- use this information to **determine the ER/HER2 status by apply- ing the molecular-subtype classifiers**.

adjuvant treatment or survival outcomes are not used to assign molecular subtype at detection

the features of the breast tumor at clinical detection are used to deter- mine survival in the absence of screening and treatment by sampling from the 1975 to 1981 SEER survival curves

we can compute (or ‘‘**back- calculate’’) survival curves for each molecular subtype**. In essence, by following this procedure, we are **generat- ing a ‘‘virtual population’’ of women, conditioning each of the tumor features and outcomes to form a sample from a distribution**

a **modeling approach to estimate ER- and HER2-specific breast cancer features and survival in the absence of screening and treatment**

data-driven as it leverages data from 2 large sets, SEER and BCSC, and builds a link between them by **making use of a model of the natural history of breast cancer**. 

In general, **ER-negative status and HER2- positive status are associated with higher tumor aggres- siveness, are harder to detect by mammography, and are more frequent (percent-wise) among younger women.** 

The underlying survival estimated through this approach suggests that **both ER and HER2-statuses are strong predictors of long-term prognosis**, **even in the absence of screening and treatment**. 

Our analysis also reveals a **crossover between ER+ and ER2 annual hazards of breast cancer death, where the latter exhibits a much higher risk of death in the first 5 y but has lower risk afterwards.**

Our **ER/HER2 estimates can serve as input for simu-lation models aimed at recreating incidence and mortality trends by molecular subtype**. All models in the CISNET Breast Cancer Working Groups have already implemen- ted these estimates to model molecular-specific incidence and mortality trends.2

our method is computationally intensive, it has several advantages. 

1. it is data-driven 
2. it maintains the observed correlations among patient age, stage, survival, and ER status. 
3. it enables us to **associate molecular subtype with the probability of being screen detected**, correcting for the length bias induced by over-sampling ofER+ cases by screen detec- tion. 
4.  is flexible enough to account for new evidence that may demonstrate a different rela- tionship among age, stage, survival and molecularly spe- cific subgroups.

Our approach does not account for ductal carcinoma in situ (DCIS). DCIS was not considered in Model S due to issues of nonidentifia- bility concerning the estimation of natural history para- meters that describe the progression of DCIS to invasive disease.

### Structure, Function, and Applications of the Georgetown–Einstein (GE) Breast Cancer Simulation Model - Schechter et al. 2018



The model is a **discrete events microsi- mulation of single-life histories of women from multiple birth cohorts**. 

Events are simulated **in the absence of screen- ing and treatment**, and **interventions are then applied to assess their impact on population breast cancer trends**

The model accommodates differences in natural history associated with estrogen receptor (ER) and human epidermal growth factor receptor 2 (HER2) biomarkers, as well as conventional breast cancer risk factors

The approach for simulating breast cancer natural history is **phenomenological**, **relying on dates, stage, and age of clinical and screen detection for a tumor molecular subtype without explicitly modeling tumor growth**.

The inputs to the model are reg- ularly updated to reflect current practice

The model has been used in collaboration with other CISNET models to assess cancer control policies and will be applied to evaluate clinical trial design, recurrence risk, and poly- genic risk-based screening.

Simulation model- ing can be a useful research method to synthesize existing data and compare a broader range of alternatives than can be feasibly included in clinical trials or other stud- ies.5,6,7

(CISNET) was launched by the National Cancer Institute in 2000 to promote open collaboration to advance modeling science. **The aim of CISNET** is to provide a range of decision makers with **tools for synthe- sizing evidence to determine the impact of alternative cancer control strategies on US population incidence and mortality**

he Georgetown University–Albert Einstein College of Medicine breast cancer simulation model (Model GE). Model GE was 1 of 7 original breast cancer simulation models.

GE pre-C, however, modeled only a single birth cohort, did not model any secular trend in age-specific breast cancer incidence, and could simulate only simple, strictly periodic screening pro- grams. 

examples of the use of GE pre-C related to evaluating strategies to improve breast cancer outcomes in African-American women8 and to determine upper age limits for mammography screening

The **initial goal of CISNET** was to have investigators with existing breast cancer **simulation models adapt them to examine the relative contributions of the disse- mination of mammography screening and the utilization of adjuvant treatment to the decline in breast cancer mortality** observed from 1975 to 2000 in the US.

The model is a continuous-time, event-driven microsimulation utilizing a parallel universes approach

The parallel universes approach starts with the **genera- tion of a basic life history for each simulated woman in the absence of any screening or adjuvant treatment**

The **effects of each screening and adjuvant treatment strategy under study are then simulated starting using the same basic life history**. In this manner, the outputs for the dif- ferent screening and adjuvant treatment strategies are **matched pairs (tuples).** 

Breast cancer **incidence depends on age, time period, and birth cohort,** and can be modi- fied based on risk. 

The incidence **includes a subset of ductal carcinoma in situ (DCIS) tumors that never sur- face clinically and eventually regress**. 

Breast cancers include **4 molecular subtypes based on estrogen receptor (ER) and human epidermal growth factor receptor 2 (HER2) status.** 

The approach to simulating breast cancer natural history is **phenomenological**, relying on **dates, stage, and age of clinical and screen detection for a tumor molecular subtype without explicitly modeling tumor growth**. 

Model GE uses the common CISNET input parameters either directly or as calibra- tion targets to conduct this modeling;17 other parameters are model-specific, 

*Model Top-Level Iteration Cycle*

Several life history event counts are provided as model output, including in each age (single years) and calendar year, the number of women alive at the start of the year, incident breast cancers (disaggre- gated by stage, ER, and HER2), mammograms (disag- gregated into true positive, false positive, true negative and false negative), deaths from breast cancer, and deaths from other causes.

Basic Life History in the Absence of Screening and Treatment. The generation of simulated individual life histories is done by random sampling from specified probability distributions. 

The basic life history object contains a date of birth, breast density at ages 40, 50, and 65 years, a date of onset of clinically diagnosed breast cancer (which may be never), and a date of death from breast cancer in the absence of adjuvant treatments (which may be never, even if breast cancer is diagnosed in the simulated woman’s lifetime). 

*Modeling Effects of Screening*
There are 2 separate aspects to modeling the effects of screening: 

1. a simulated screening schedule; through calling 2 functions.
   1. uses the simulated woman’s birth year to sample a date of first mammogram.
   2. another func- tion that samples the next mammogram date and updates the running information about the woman’s screening history
2. screen detection.

Stage Shift Due to Screen Detection. **Model GE works backwards** in time from the date at which the lesion would have been diagnosed clinically and its stage at that time, to the time at which a true-positive mammogram is obtained. 

*Simulated Treatment*
When simulating a strategy that includes adjuvant treat- ment, the simulation must identify a treatment to apply, and a possible modification of underlying survival according to the treatment. 

We have simulated 3 treat- ment approaches: 

1. no adjuvant treatment
2. dissemination of adjuvant treatment (intended to reflect adjuvant ther- apy as actually used in the US since its introduction in the 1980s), 
3.  ‘‘optimal’’ treatment (intended to repre- sent the most effective therapy available for the woman at the time she is diagnosed). 

The dissemination and opti- mal treatment strategies may result in the application of no adjuvant therapy to some women. 

The probability of the application of each adjuvant treatment to a woman is conditional on calendar year, age, stage, and ER and HER2 status, based on CISNET common inputs.

The optimal treatment strategy is also a common input, cal- culated by selecting the age, stage, ER-, and HER2- specific treatment associated with the greatest prolonga- tion of survival among therapies available in the year of diagnosis.

### Reflecting on 20 years of breast cancer modeling in CISNET: Recommendations for future cancer systems modeling efforts - Trentham-Dietz et al. 2021

Since 2000, the National Cancer Institute’s Cancer Intervention and Surveillance Modeling Network (CISNET) modeling teams have developed and applied microsimulation and statis- tical models of breast cancer

. The 6 CISNET breast cancer models embody the key features of systems modeling by incorporating numerous data sources and reflecting tumor, person, and health system factors that change over time and interact to affect the burden of breast cancer. 

Our 6 CISNET breast cancer models embody the key features ofsystems modeling by incorporating numerous data sources and reflecting tumor, person, and health system factors that change over time and interact to represent the burden ofbreast cancer

investigate questions related to breast cancer biology, compare strategies to improve the balance ofbenefits and harms ofscreening mammography, and support insights into the delivery ofcare by modeling outcomes following clinical decisions about breast cancer treatment

(CISNET) modeling teams have developed and applied microsimulation and statistical models of several types ofcan- cer, including breast cancer [4].

. CISNET breast cancer models incorporate data on distribu- tions oftumor characteristics, women’s risk factors, and healthcare use of breast cancer control interventions

used to evaluate the impact of various screening and treatment interventions on multiple health outcomes in the overall United States population and population subgroups that differ by race, risk, and/or breast density

Each modeling group begins with a common set ofinputs. 

*Tumour level*

**The models estimate age-specific incidence for first diagnosis of breast cancer overall and by molecular subtype over time and by birth cohort in the absence ofscreening**

Second breast cancers including recurrences and new primary breast cancers are not yet modeled but are currently being added for forthcoming studies. 

Five models use common **estimates of breast cancer incidence rates in the absence of screening (“background incidence rate”) for each calendar year and single year ofage derived from an age–period–cohort model**

One model (**Model M) assumes a linear model for the annual incidence rates during the years 1975 to 2012** under a hypothetical scenario of no screening **then adds the effect of screening dissemination patterns** to the linear model. A Bayesian approach is applied to adjust the linear model parameters so that the model output matches the SEER rates in 1975 to 2012

All 6 models consider **4 breast cancer molecular subtypes** based on age-specific proportions of **ER and HER2** positive and negative breast cancers in the population

**Stage** of breast cancer at detection with and without screening is based on data from the Breast Cancer Surveillance Consortium (BCSC) [8]. 

Tumors are assigned **sojourn times** (defined as the time from when tumors are detectable by screening until they are detectable by clinical symptoms) or mean tumor **doubling times** specific to subtype conditioning on age group (�40, 40 to 50, and �50) [22]. 

**Tumor growth** for model WH is based on a Gompertz-type function with a lag to account for the difference in timing ofdetection by screening or symptoms [21]. 

Tumors can be detected at a younger age and earlier stage (or smaller size) if they are screen detected than if they are clinically detected

In all models, some tumors could be considered overdiagnoses, where **overdiagnosis** is defined as screen-detected cancers that would not have been diagnosed within the woman’s lifetime in the absence ofscreening. These tumors have no effect on breast cancer–specific mortality

*Individual person level*

The models **simulate the life history of each individual woman until death**. 

At the start ofthe simulation, a woman is assigned a date at birth and date ofdeath from other-cause mortality. Based on the incidence parameters described above, **some women develop breast cancer** and are assigned a date (and age) ofsymptomatic clinical detection for breast cancer in the absence ofscreening. If the date ofbreast cancer is before the date ofother-cause death, the cancer can be screen or clinically detected. Women can die ofnon-breast cancer causes at any time; non- breast cancer mortality rates by age and calendar year are derived from national data [29]. 

**Breast density** is a radiographic feature observed on mammogram images that reflects the degree to which fibroglandular tissue is radio-lucent (white on the image) or fatty and radio- opaque (dark on the image). Women are assigned one of **4 breast density levels** [Breast Imag- ing Reporting and Data System (BI-RADS) categories: almost entirely fat, scattered fibrogland- ular density, heterogeneously dense, or extremely dense] at age 40 [33]. Women are assigned to either the same breast density category or to the next lower category at ages 50 and 65 based on observed age-specific prevalence in the BCSC [34,35]. We assume density does not change after after age 65.

*Health system level*
Although women are assigned **receipt ofscreening or therapy** at the individual person level in the models, these interventions are delivered via the health system. 

**assigned an age at first mammogram and screening frequency** based on the distribution observed for their birth cohort using data from the BCSC, the National Health Interview Survey, and the US Food and Drug Administration’s Mammography Quality Stan- dards Act and Program [24].

Model inputs for mammography performance are based on data from the BCSC and depend on a woman’s age (25 to 39, 40 to 49, 50 to 64, and �65), density level, screening interval (annual, biennial, and triennial/infrequent), and whether the mammo- gram was the first screening or not [8,24]. 

All women diagnosed with breast cancer are **assumed to receive initial therapy** with mastectomy or lumpectomy with radiation, but local therapy is not explicitly modeled. Subtype-spe- cific **adjuvant treatment** (chemotherapy, endocrine therapy, and trastuzumab) is assigned according to a dissemination model based on SEER patterns ofcare special studies (1980 to 1996) and the National Comprehensive Cancer Network Outcomes Database (1997 to 2012) [24]. 

Breast cancer **survival depends on age group** (<40, 40 to 49, 50 to 59, 60 to 69, and 70 to 84), and American Joint Committee on Cancer (AJCC)/SEER **stage or tumor size** in the absence ofscreening and treatment (background survival) as estimated from our prior research [8,22]. 

Systemic treatment reduces the hazards ofbreast cancer death (Models D, GE, M, and S) or results in cure for some cases (Models E, WH) based on age-specific data from the most current Oxford Overview of clinical trials [32]. Since the Overview did not find age differences in efficacy, hazards reductions are applied to all age group

*Model output analysis*

 The population ofUS women is modeled **starting in the year 1975 until the most recent year of data** available in the SEER Program. 

- Either **the entire population can be modeled (all ages and birth cohorts) or **
- **a single recent birth cohort can be selected to simulate expected outcomes of a contemporary cohort ofwomen receiving the current standard ofcare for breast cancer screening and treatment**. 

The models generate a wide range ofbenefit and harm outcomes

**Mortality reductions can be attributable to screening alone, treatment alone, or the combination in a given calendar year** by calculating the difference between the mortality rates pre- dicted with an intervention and the background mortality rate in the absence ofscreening and treatment. 

There are several methods for these estimates that **consider the potential for nega- tive synergy** between the contributions ofscreening and systemic treatment on mortality reductions (i.e., as treatment becomes more effective at later stages, the contribution ofscreen- ing to mortality reductions decreases)

### Experimentally-driven mathematical modeling to improve combination targeted and cytotoxic therapy for HER2+ breast cancer - Jarrett et al. 2019

experimentally and computationally investigate combination trastuzumab- paclitaxel therapies and identify potential synergistic effects due to sequencing of the therapies with in vitro imaging and mathematical modeling

Longitudinal alterations in cell confluence are reported for an in vitro model of BT474 HER2+ breast cancer cells following various dosages and timings of paclitaxel and trastuzumab combination regimens.

Results of combination drug regimens are evaluated for drug interaction relationships based on order, timing, and quantity of dose of the drugs.

. **Two mathematical models** are introduced that are **constrained by the in vitro data to simulate the tumor cell response to the individual therapies.**

A **collective model** merging the two individual drug response models was designed to investigate **the potential mechanisms of synergy** for paclitaxel-trastuzumab combinations.

The synergy derived from the model is found to be **in agreement with the combination index**, where both indicate a spectrum of additive and synergistic interactions between the two drugs dependent on their dose order.

breast cancer patients receive systemic delivery of targeted and cytotoxic drugs (in parallel or sequentially) for various subgroups of receptor positive breast cancers.

The additive effects of introducing targeted therapy to cytotoxic treatment have been explored experimentally; however, the synergistic (or even antagonistic) effects of sequencing these regimens in breast cancer have not been systematically investigated. This is likely due to the enormity of such a study as there are nearly limitless combinations of dosing, timing, and ordering of treatments indicated for any given subtype of breast cancer

, **experiment-driven, mathematical modeling** could alleviate these challenges by **investigating a myriad of combination therapy strategies in silico to identify potential treatment regimens for focused in vivo and in vitro investigations**

standard combination therapy for the treatment of breast cancers that overexpress the human epidermal growth factor receptor 2 (HER2) is the simultaneous administration of paclitaxel and trastuzumab. HER2

d 25–30% of all breast cancer cases are considered HER2+4,5, which is associated with poorer overall prognoses with more aggressive

Paclitaxel, a chemotherapy, causes cell death by stabilizing microtubules during mitosis—impeding normal cytokinesis and equal cellular divisions and proliferation6.

Trastuzumab binding blocks HER2/neu, preventing receptor dimerization and interfering with intracellular signaling8. This obstruction induces cell cycle arrest, inhibits cell proliferation and migration9–11, and causes HER2 internalization and subsequent degradation

targeted therapy in combination with systemic cytotoxic therapy has been shown to increase overall survival rates, nearly 26% of patients have recurring disease within 10 years14. 

Since the clinical initiation of trastuzumab, cytotoxic drugs in combination with targeted anti-HER2 therapy, has remained the standard-of-care practice for treatment of primary HER2+ breast cancer. 

 in vivo preclinical animal data suggests the order, dosing, and timing of combination therapy has yet to be optimized1

time-resolved **microscopy data that captures changes in in vitro cell confluence in response to combination paclitaxel and trastuzumab therapy**.

evaluated for **synergism versus additive or antagonistic relationships** based on order, timing, and quantity of dose of the two drugs. 

**Mathematical models** motivated from and **constrained by single treatment in vitro data** to simulate and predict the tumor cell response are then derived. 

a **collective model merging the two individual drug response models** is subsequently designed to **reveal potential synergistic and non-synergistic (additive and antagonistic) effects** between trastuzumab and paclitaxel due to dosage and timing of the therapies.

**Mathematical models  to describe the response of tumor cells to each of the single agent therapeutic regimens**. 

The models are **ordinary differential equations (ODEs),** which describe **how the tumor cell number changes as a function of time and drug concentration.**

matlab model https://github.com/OncologyModelingGroup/InVitroPacandTmAbSynergy

All of the parameters of the models were calibrated to the time-resolved microcopy data (apart from HER2exp, which is assigned based on the cell line17,20,21).

Synergy for combination regimens in the mathematical model. 

First, the combined mathematical model was calibrated for growth and carrying capacity values using data acquired during the initial 24 hours.

Then the model was simulated for the different combination regimens using the corresponding values derived from the single drug dose results to determine a baseline performance of the model where synergy was not considered. 

Then, to explore possible interaction effects between the two drugs, the parameter S was allowed to calibrate (from its nominal value of 1, which represents no synergy) for the combination regimens.

The sequen- tial calibration strategy, as described above for the single agent regimens was used;

Additionally, the data from day 1 to day 2 was used to perform a limited recalibration of the parameters associated with the first drug applied; s

the **combination index** (using the Loewe additivity null reference principle23) was calculated for the different regimens at day 4 to determine if the combined effects of the two drugs are synergistic, additive, or antagonistic:

The **concordance correlation coefficient** (CCC) was calculated to compare the overall agreement in temporally
dynamic trends of the model simulations to the in vitro da

Calculation of the combination index shows the change in synergism quantitatively with paclitaxel first resulting in only a trend toward an additive effect. 

While the two synergy measures agree (the combination index and the calibrated values of S), only the results
from the mechanistic mathematical modeling are hypothesis generating in regard to the biological phenom- ena behind this effect.

Many studies have focused on the **quantification and definition of synergistic effects** of drugs1,23,27–30, and the results of these studies are useful for identifying the presence and measuring the amount of synergy between drugs. However, the results **cannot be used to help understand the molecular/cellular drivers of synergy because expressions for drug mechanisms have not been incorporated**. 

Our work utilizes **in vitro cellular scale data to describe the synergism between the two drugs as well as elicit information about the biological cause of this synergy**—**biologically identifying mechanisms of drug synergy is difficult and data-driven**. 

Mathematical models can be used to recognize potential cellular mechanisms from where synergy might be derived. Other **mechanistic modeling efforts for synergy of anti-cancer therapeutics have directly considered signaling networks with kinetic, partial differential equation, agent-based, and even multi-scale models**–34. However, these studies are limited by their **computational expensive and data demanding nature** (requiring large amounts of time course proteomics and/or transcriptomics) and **do not consider synergism** dependent on the sequencing or timing of drugs.

### Development and Validation of a Simulation Model-Based Clinical Decision Tool: Identifying Patients Where 21-Gene Recurrence Score Testing May Change Decisions - Jayasekera et al. 2021

a need for industry-independent **decision tools** that integrate clinicopathologic features, comorbidities, and genomic information for women with node-negative, invasive, hormone receptor–positive, human epidermal growth factor receptor-2–negative (early-stage) breast cancer

The goal of this study was to develop and validate an independent **clinical decision tool** that **integrates individual, clinical, and pathologic characteristics with 21-gene recurrence score assay information** to **guide adjuvant chemotherapy treatment decisions** in women diagnosed with node-negative, invasive, hormone receptor–positive, human epidermal growth factor receptor-2–negative (early-stage) breast cancer

**simulation model–based web clinical decision tool** (BTxChoice) provides the 10-year risk of distant recurrence and life-years gained with chemoendocrine versus endocrine therapy considering a woman’s age, tumor size, tumor grade, and comorbidities with and without 21-gene recurrence score test results.

, we **adapted and validated a CISNET a breast cancer simulation model to power a new clinical decision tool** (BTxChoice, Washington, DC).15-

This independent tool provides estimates of 

- 10-year risk of distant recurrence, 
- breast cancer–specific mortality, 
- other-cause mortality, and 
- life-years gained with chemoendocrine versus endocrine therapy 

for individual women based on 

- age, 
- tumor size, 
- grade, and 
- comorbidity level with and without 21-gene RS results. This tool is flexible, includes the range of uncer- tainty for each outcome, fills current clinical needs, and facilitates personalized treatment decisions in HR1, HER2–, node-negative breast cancer.

*Model Inputs*

Individual and clinicopathologic characteristics (eg, age, grade, tumor size, estrogen re- ceptor or progesterone receptor status, and RS)

Time-to-distant recurrence was defined as time from di- agnosis to tumor recurrence at a distant site according to the Standardized Definitions for Efficacy End Points in Adjuvant Breast Cancer Trials (STEEP) criteria

Contralateral disease, second primary cancers, and death without distant recurrence were considered cen- soring events.

*Model Overview* 

The simulation model uses an **empiric Bayesian analytical approach that captures uncertainty in all predictors’ effects on outcomes and the sampling variation**

modeled women from the point of having newly diagnosed breast cancer to death.

The population were women with HR1, HER2–, invasive, node-negative breast cancer with tumor size <= 5.0 cm that have received lumpectomy (with ra- diotherapy) or mastectomy.

First, we **simulated individual characteristics for each vir- tual woman** based on the joint distributions of age, tumor size, grade, hormonal status (estrogen receptor or pro- gesterone receptor), comorbidity level, and RS.16

Each woman could remain event-free, experience a distant re- currence, die of breast cancer or other causes conditional on her treatment (5 years of endocrine therapy or che- moendocrine therapy), age, RS, tumor size, grade, and comorbidity level. A woman with distant recurrence could die of breast cancer or other causes. 

Time-to-events for each virtual woman was identified by applying the sub- hazard ratios of the predictive attributes to the baseline cumulative incidence functions estimated using the competing-risk models described above. Distant recur- rence and breast cancer death were modeled separately because of limited data that were available to model the distribution of time from distant recurrence to breast cancer death (Data Supplement). We simulated the effects of endocrine versus chemo- endocrine therapy for 1,512 subgroups defined by com- binations of age (within 5-year bands), tumor size (# 2cm or . 2 cm), grade (low; intermediate; or high), and comorbidity level (no; mild; moderate; or severe), all with RS (0-100 in increments of 5) or without RS.16

We chose 5-year age bands based on current computational capacity; using individual years of age and RS would have required more than 87,000 subgroups.

[RS = recurrence score]

First, we calculated the proportion of women belonging to each RS category within each subgroup defined by age (5- year bands), tumor size, and grade. 

Next, we estimated the 10-year distant recurrence rates and 95% CIs for che- moendocrine versus endocrine therapy using Kaplan- Meier curves for each subgroup with clinical-pathologic features alone (age, tumor size, and grade), and then re- peated the analysis with clinical-pathologic features com- bined with RS results. 

Since chemotherapy generally prevents most recurrences within 5 years after diagnosis,6 only recurrences up to 10 years postdiagnosis were estimated

Average life-years were calculated from the date of diag- nosis to date of breast cancer death or age- or comorbidity- specific other-cause death with chemoendocrine or en- docrine therapy. Life-years gained were estimated by calculating the difference in average life-years for che- moendocrine therapy versus endocrine therapy for each subgroup.

We quantified the uncertainty related to sampling variability for any given input parameter value by replicating the simulation for each of the 1,512 subgroups up to 1,000 times

Independent validation of results was performed to confirm model accuracy.

This novel simulation model–based clinical decision tool can be useful in clinical practice in several ways. 

1. can be used to guide treatment decisions even if RS testing is not available
2. provides additional information about RS-category– specific benefits, given a woman’s age, comorbidities, and clinicopathologic characteristics. 
3. provides treatment out- comes based on a woman’s age- and comorbidity-specific life expectancy
4. pro-vides individualized estimates of absolute chemotherapy benefit 

