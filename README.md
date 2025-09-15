# Rethinking Item Fairness Using Single World Intervention Graphs

Youmi Suk^1 and Weicong Lyu^2

1: Teachers College, Columbia University
2: University of Macau

## Overview

Since the 1960s, the testing community has striven to ensure fairness in tests and test items. Differential item functioning (DIF) is a widely used statistical notion for items that may unfairly disadvantage specific subgroups of test-takers. However, traditional DIF analyses focus only on statistical relationships in observed data and cannot explain why such unfairness occurs. To fill this gap, we introduce a novel causal framework for defining and detecting unfair items using single world intervention graphs (SWIGs). By leveraging SWIGs and potential outcomes, we define causal DIF (CDIF) as the difference in item functioning between two hypothetical worlds: one where individuals were assigned to one group and another where they were assigned to a different group. We also connect CDIF to related fairness concepts, including group versus individual fairness and item impact. In particular, we use SWIGs to graphically distinguish between item fairness at the individual level and the population level. Additionally, we discuss causal identification strategies using SWIGs and demonstrate how our approach differs from traditional DIF methods through a simulation study. We further illustrate its application to a controversial item from New York's Regents math exam and a real dataset from the Trends in International Mathematics and Science Study, concluding with broader implications of promoting causal fairness in testing practices.

For more details of our proposed framework, see our paper: [(Preprint)](https://osf.io/preprints/psyarxiv/bue62_v2). 
Here, we provide `R` codes to reproduce our simulation study and replicate our data analysis using the the Trends in International Mathematics and Science Study (TIMSS) 2011 data. 

## Simulation Study

* `aux.R`  

   This `R` file includes data generating codes and an outcome regression function that detects CDIF.

* `simu_D1.R`, `simu_D2.R`, `simu_D2.R`
 
   These `R` file includes simulation codes under simulation designs 1, 2, and 3. For more information on designs and evaluation, see our paper: [(Preprint)](https://osf.io/preprints/psyarxiv/bue62_v2). 


## TIMSS Data Study

* `dat.RData`

  [ADD Explanation].

* `analysis.R` 
 
  [ADD Explanation].
  
