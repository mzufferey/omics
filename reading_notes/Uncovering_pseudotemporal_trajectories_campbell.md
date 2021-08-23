### Uncovering pseudotemporal trajectories with covariates from single cell and bulk expression data (Campbell & Yau 2018)



**Pseudotime algorithms**: employed to extract latent temporal information from cross-sectional data sets allowing dynamic biological processes to be studied in situations where the collection of time series data is challenging or prohibitive



**PhenoPath**: a novel statistical framework (**hybrid regression-latent variable model**) that **learns how pseudotime trajectories can be modulated through covariates that encode such factors**.



* **longitudinal studies** are often challenging to conduct and cohort sizes limited by logistical and resource availability
* **cross-sectional surveys** of a population are relatively easier to conduct in large numbers and more prevalent for molecular ‘omics based studies. 
  * do not directly capture the changes in disease characteristics in patients but it may be possible to recapitulate aspects of temporal variation by applying **“pseudotime” computational analysis**



The objective of **pseudotime analysis** is to take a collection of high-dimensional molecular data from a cross-sectional cohort of individuals and to map these on to a series of one-dimensional quantities, called **pseudotimes**



These **pseudotimes** **measure the relative progression of each of the individuals along the biological process of interest**, e.g., disease progression, cellular development, etc., allowing us to **understand the (pseudo)temporal behaviour of measured features without explicit time series data**

- possible when individuals in the cross-sectional cohort behave asynchronously and each is at a different stage of progression
- by creating a relative ordering of the individuals, we can define a series of molecular states that constitute a trajectory for the process of interest



Pseudotime methods generally rely on the **assumption that any two individuals with similar observations should carry correspondingly similar pseudotimes** and algorithms will attempt to **find some ordering of the individuals that satisfies some overall global measure** that best adheres to this assumption

- differ in the way “similarity” is defined
- when applied to molecular data, typically capture some dominant mode of variation that corresponds to the continuous (de)activation of a set of biological pathways



gained particular popularity in the domain of **single-cell** gene expression analysis (where **each “individual” is now a single cell**) e.g. to model the differentiation (cf. https://github.com/agitter/single-cell-pseudotime)

- use advanced machine learning techniques (e.g. can characterise cell cycle,  model branching behaviours)



these single-cell applications were **predated by more general applications** in modeling disease progression

- provided early inspiration for single-cell pseudotime methods



To date, little cross-over between these distinct application domains (different contexts of application)

- interesting possibilities by translating recent advances in single-cell pseudotime modelling to disease progression modelling

  



recent single-cell pseudotime approaches for branching pseudotime trajectories, these **can only be retrospectively examined for their association with prior factors of interest****

$\rightarrow$ develop **a statistical model in which these factors could be explicitly incorporated** into pseudotime analysis

- would provide **a mechanism to account for known genetic, phenotypic or environmental factors allowing gene expression variability to be decomposed into different contributory factors**
- would allow us to answer questions related to the **interaction between heterogeneity in these external factors and biological progression**



a novel **Bayesian statistical framework for pseudotime trajectory modelling that allows explicit inclusion of prior factors of interest**

- **allows to incorporate information in the form of covariates** that can modulate the pseudo- temporal progression allowing sub-groups within the cross- sectional population to each develop their own trajectory
- **combines linear regression and latent variable modelling and allows for interactions between the covariates and temporally driven components** of the model
- first approach to allow for **modelling pseudotime trajectories on heterogeneous backgrounds** allowing its **utility in both single and non-single cell** applications



**PhenoPath** provides a **probabilistic ordering of high-dimensional gene expression measurements across objects** (e.g., cells, tumours, patients, etc)

- achieved by **compressing** the information contained within the data on to a **unidimensional axis**
  - construct an axis such that **relative positions along the axis correspond to some meaningful biological or disease progression**
- novelty: introduce the notion that **objects may have different labels (covariates)** attached to them corresponding to different innate properties or exposure to external stimuli
  - these factors might cause the objects to evolve over (pseudo)time differently
-  **simultaneously learns a pseudotemporal axis** that is common to the different object labels, **while decomposing gene expression variability into static and dynamic components**

$\Rightarrow$ a **Bayesian** statistical framework that integrates **linear regression and latent variable modelling**

- the observed data ($y_n$) for the nth individual is a linear function of both measured covariates ($x_n$) and an unobserved latent variable ($z_n$) corresponding to latent progression that we will term pseudotime



the model involves **three components**:

1. **gene expression**:  a **static** component based on your covariate status ($Ax^T_n $)

2. a **dynamic** component related to **how far** along the biological process you are ($λz_n$)

3. (main novelty) an **interaction component** which allows your **covariate status to change the direction of the dynamic component** of the gene expression ($Bx^T_nz_n$)

   

* if only 1. used =  linear regression based differential expression analysis
* if only 2. used = factor analysis 



the **covariates** in this study are binary quantities, any arbitrary design matrix that can be used for standard regression may be used for $x$ (

**sparse Bayesian prior probability distributions** are used to constrain the parameters ($A$, $B$, $\lambda$) so that covariates only drive the emergence of distinct trajectories if there is sufficient information within the data to do so



**fast and highly scalable** variational Bayesian inference framework that can handle thousands of features and samples in minutes using a standard personal computer 



variational inference of such hierarchical Bayesian models can be sensitive to **hyperparameters values and parameter initialisation** we found PhenoPath to be **robust** to such choices 

