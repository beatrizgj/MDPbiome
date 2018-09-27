# MDPbiome
*Authors: Tomás de la Rosa, Beatriz García-Jiménez, Mark D. Wilkinson*


### Description
**MDPbiome** is a software developed in R that uses Markov Decision Processes to create "policy prescriptions" for microbiome engineering. MDPbiome performs a variety of analysis describing the robustness of the prescription, as well as a variety of visualizations to assist in manual interpretation of the state transitions and biological understanding of a microbiome's dynamics.

**Citation**: Please, cite us when you use MDPbiome in your work:  

**Beatriz García-Jiménez, Tomás de la Rosa, Mark D. Wilkinson; MDPbiome: microbiome engineering through prescriptive perturbations, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i838–i847, [doi: 10.1093/bioinformatics/bty562](https://doi.org/10.1093/bioinformatics/bty562)**


### Pre-requisites 
MDPbiome uses the **phyloseq** package to manage data from OTU tables, and the **MDPtoolbox** package to perform MDP algorithms.  In addition, the system needs commonly used R package for data analysis, such as reshape, ggplot2 and Hmisc.  For a complete list, check the imported libraries in the initMDPBiome.R file.

### Dataset Setup
Given the heterogeneity of longitudinal microbiome datasets, MDPbiome needs to pre-load a configuration script file where user specifies the elements and features of the dataset. The input to MDPbiome is a phyloseq object obtained from importing an OTU table and a mapping file containing the meta-data of the microbiome samples. 
This mapping file, in CSV or TSV format, should contain the following columns:

* a variable to identify subjects (e.g., "Subject.ID")
* a variable to identify time steps for sample time series (e.g., "day")
* a variable to indicate the cluster of each sample (named as "cluster")
* variables to indicate individual or combined perturbations (e.g. "antibiotic", "probiotic")
* (optional) a variable to represent the "goodness" of samples. The per-cluster average of this feature could serve as state utility function for defining reward functions.

If the mapping file does not contain one of these elements, we suggest to append them within the configuration file.  Next, we will show a step-by-step description for creating the configuration file and for running high-level functions of MDPbiome.

#### Loading Data
First, you set the location of MDPbiome as your working directory, load the sources from the initialization file.

```{r eval=FALSE}
source("initMDPBiome.R")
```

To load the OTU table and the mapping file, use the appropriate phyloseq functions

```{r eval=FALSE}
biom.otu<-import_biom("../Data/Example/otu_table.biom")
map<-import_qiime_sample_data("../Data/Example/mapping.tsv")
data.raw <- merge_phyloseq(biom.otu,map)
```

#### Normalization
We propose to normalize microbiome samples following the approach of David et. al. 2014, included as part of the software.

```{r eval=FALSE}
data.norm <- MDPpreprocess(data.raw, processDir = "../Data/Example/")

```

#### Clustering for Creating MDP States
This step involves creating MDP states as a form of abstracting groups of similar microbiome samples.
The common way of doing this abstraction is by clustering techniques.  

Note that, if you do not have a cluster assignment in your mapping file, you could use our *Robust Clustering* procedure, also provided with MDPbiome sources. This will compute a good performing set of clusters based on several metrics.  If the process doesn't find clusters it returns -1. If the process succeeds it will append a "cluster" attribute to the sample data of the phyloseq object. 

```{r eval=FALSE}
source("robust.clustering.metagenomics.functions.r")
clus.results <- robust.clustering(data.norm, "age", "Ballou2016")
data.norm <- robust.clustering.decision(data.norm,clus.results,"age","Ballou2016")
```

#### Additional Attribute Generation and Final Setup
For output generation we recommend to have short names for MDP actions (distinct values of perturbation). In some cases this implies additional pre-processing of the mapping data, for instance, renaming as "yes"/"no" the occurrence of a perturbation in a sample. For simplicity we omit these steps, since they are not strictly necessary. This is also the point where you may want to create a combined perturbation, for instances appending the labels of two single perturbations.

Then, you should specify the set of perturbations and the required variables

```{r eval=FALSE}
# the set of perturbations
Perturbations <- c('Combined','Salmonella','Probiotic')

# sample attribute for time points
stepVar <- "age"

# sample attribute for subject 
subjectVar <- "Subject.Id"

# default utility function
goalVar <- "diversity"

```

Finally, you need to set a path to the location where output files will be generated. The *createTreeDir* function prepares a directory structure considering perturbations and high-level functions

```{r eval=FALSE}
dirdata <- "../Data/Example/"
createTreeDir(dirdata,Perturbations)
```

### Running MDPbiome
We describe here the high-level functions of the system. To see examples of the output you can expect, please check the pre-computed results presented in [MDPbiome results](https://tomdelarosa.shinyapps.io/mdpbiome/)

<!-- #### 1. Analyzing Sample Time Series -->
<!-- This generates set of plots for analyzing sample time series. You need to provide the phyloseq object with the required elements described above, and optionally the phyloseq object without the normalization step.  This second object is used to compute the alpha-diversity of samples using original abundances. -->

<!-- ```{r eval=FALSE} -->
<!-- mdpBiomePreAnalysis(data.norm,data.raw) -->
<!-- ``` -->

#### 1. Computing Optimal Policy & Policy Reliability
This computes the optimal policy for the given data, and performs a "Monte Carlo"-like sampling experiment to determine the stability ratio for the policy and by each individual action. By default, the alpha-diversity criteria is used to determine the reward function.

```{r eval=FALSE}
mdpBiomeBase()
```

Alternatively, you can provide an attribute indicating the goodness of the sample.  The per-cluster average of this attribute can serve as the utility vector U(s), which is the basis for reward schemas.

```{r eval=FALSE}
mdpBiomeBase(utilityVar="altDiversity")
```

#### 2. Computing Policy Generality
This computes a *leave-one-out* cross validation to determine how general policies are when they are prescribed
to subjects not included in the modelling data.  The performance is measured comparing the outcome when subject follow or not the proposed optimal policy

```{r eval=FALSE}
mdpBiomeLoocv("preferGood",goalVar = "diversity")
```
In this example, we select the "preferGood" reward schema and use the alpha diversity as a criteria for assigning utility to the clusters.








