Introduction to Lethals Project
============

This document was created as a repository of the scripts used in the analyses of the article: XXXXXX (doi: XXXXX).
Files used in the pipeline are described as follow
- *sim.slim* Forward in time simulations were conducted using the software [SLiM 3](https://github.com/MesserLab/SLiM)<sup>1</sup>.
- *demographic_inference.py* This python script contains the code used for the inference demographic parameters using the software [dadi](https://dadi.readthedocs.io/en/latest/#welcome-to-dadi)<sup>2</sup>.
- *selection_inference.py* This python script contains the code used for the inference of the distribution of fitness effects (DFE) using the software [fitdadi](https://github.com/LohmuellerLab/fitdadi)<sup>3</sup>.

Scripts for each figure are:
- *figures/script/exp_v_inferred.R* : Figure 2B, 5B 
- *figures/script/dot_plot.R* : Figure 2A, 5A 
- *figures/script/log_profile.R* : Figure 6, S5 
- *figures/script/compare_ss.R* : Figure 3, S4
- *figures/script/compare_sfs.R* : Figure 1, 4, S3
- *figures/script/compare_inf_sim_sfs.R* : Figure S2

Data produced to generate figures are:
- *figures/data/all_inferences.csv* 
- *figures/data/all_log_inferences.csv* 
- *figures/data/exp_v_inferred.csv* 
- *figures/data/sim_avg_sfs10.csv* 
- *figures/data/sim_avg_sfs100.csv*  
- *figures/data/sim_avg_sfs1000.csv* 
- *figures/data/computed_avg_sfs10.csv* 
- *figures/data/computed_avg_sfs100.csv*  
- *figures/data/computed_avg_sfs1000.csv* 

Scripts for mutation-selection-drift balance results are available at: https://github.com/ckyriazis/lethals_scripts

**References**
- <sup>1</sup> Haller, Benjamin C., and Philipp W. Messer. 2019. “SLiM 3: Forward Genetic Simulations Beyond the Wright–Fisher Model.” *Molecular Biology and Evolution* 36 (3): 632–37.
- <sup>2</sup> Gutenkunst, Ryan, Ryan Hernandez, Scott Williamson, and Carlos Bustamante. 2010. “Diffusion Approximations for Demographic Inference: DaDi.” *Nature Precedings*, June, 1–1.
- <sup>3</sup> Kim, Bernard Y., Christian D. Huber, and Kirk E. Lohmueller. 2017. “Inference of the Distribution of Selection Coefficients for New Nonsynonymous Mutations Using Large Samples.” Genetics 206 (1): 345–61.
