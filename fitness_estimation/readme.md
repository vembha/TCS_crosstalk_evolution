# Fitness estimation

Files in the subfolders of this folder helps us estimate the fitness for all the possible phenotypes for a bacterium with _N_ TCSs. Workflow is the following:
 - First, choose the environment in which fitness values are to be estimated, and access the `bossfile.m` file in the folder. This calls the file `*evaluator.m` file which performs the calculations.
 - Here, input the value of _N_, the number of TCSs, `time_diff`, the time gap between the two signals, and `decay_factor`, the time for signal to decay 1000-fold of the initial value in an exponential fashion.
 - Then, choose the number of values of γ using `num_gamma` for which you wish to estimate the fitness.
 - Run the file. The final file will be saved with the name `Non_path_fitness_for_N_*.mat` or `Fitness_for_N_*.mat` depending on whether the estimations are for random or programmed environments, respectively. It contains a matrix of `fitness` over values of γ and indexed phenotypes, the vector `max_index` which denotes the phenotype index which has the highest fitness at that value of γ, and the vector `max_fit` denoting the maximum fitness value achieved by that phenotype.

The file `K_matrix_assignment.m` gives the interaction matrix as the output for a specific value of `N`, the number of TCSs, `K_ratio`, the value of γ, and `count`, the phenotype's index. Such a matrix can be used as the input, `KC`, along with `N` and `K_ratio`, to get the phenotype identity number. These identification numbers are unique.

Finally, the files `data_set.m` and `initialization.m` are the files that have the ODE model and the initial condition vectors respectively. We use the mid-point numerical integration algorithm coded in the `mid_point_int.m` file. 
