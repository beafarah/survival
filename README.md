# Survival Analysis
This repository gathers survival analysis codes developed in the context of an M2 internship on Sample size calculation based on differences of quantiles from censored data. This work was developed under joint supervision of Aurélien Latouche (U900, Institut Curie) and Olivier Bouaziz (MAP5, Université Paris Cité).

This repository is organised as follows:
- _logrank_: contains the codes for the implementation of logrank test using counting process approach
- _kosorok_: contains the implementation of the method introduced in the paper "Two-sample quantile tests under general conditions" (Kosorok, 1999). This approach aims at testing the null hypothesis that the quantiles between two independent samples are equal, in the setting of right-censored survival data.

# References
Kosorok, Michael R. "Two-sample quantile tests under general conditions." Biometrika 86.4 (1999): 909-921.

SAS implementation of the Two-Sample Quantile Test for Right Censored Survival Data (available on https://mkosorok.web.unc.edu/qtest/)
