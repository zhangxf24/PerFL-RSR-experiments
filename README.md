# PerFL-RSR-experiments 

Code to reproduce simulation results, figures, and real data analysis results from the paper "Robust Personalized Federated Learning with Sparse Penalization" by Weidong Liu, Xiaojun Mao, Xiaofei Zhang and Xin Zhang. 

## Organization

### simulation-code  

Simulation_6.1.1.R evaluates effect of grouping and produces the results displayed in Table 1. 

Simulation_6.1.2_p.R, Simulation_6.1.2_n.R, and Simulation_6.1.2_M.R evaluates effect of heavy-tail noise and produce the results displayed  in Figure 3.  

Simulation_6.1.3.R evaluates the effect of tuning parameters and produces the results displayed in Figure 4.

Simulation_6.1.4.R evaluates effect of partial participation and produces the results displayed in Figure 5.

Functions_PerFL-RSR.R contains usr defined functions for implementing of the the proposed method.

Functions_other_models.R contains functions for implementing the compared methods in the paper.

Functions_regularizer.R contains functions for the regularization.

Functions_evaluation.R contains evaluation functions.


### simulation-results  

Contains all the results from running the simulation codes as described above. 

### figures  

Produces Table 1-3, Table S.1, and Figure 3-6 based on the simulation results.

### figures-code  

Contains codes for generating numbers in Table 1-3 and Table S.1, and plotting Figure 3-5.

### data

Contains the Communities and Crime dataset (crimedata.csv) and county information for communities (county.csv).

### real-data-code  

crime_county.R runs the real data analysis and produces results displayed in Table 3 and Figure 6.

county_data_processing.R loads data, combines county information and select useful data.

rea_func.R contains functions of implementing the real data analysis. 


