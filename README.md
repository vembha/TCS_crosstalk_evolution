# **CODE REPOSITORY FOR THE PUBLICATION**
## An evolutionary paradigm favoring crosstalk between bacterial two-component signaling systems
### Bharadwaj Vemparala, Arjun Valiya Parambathu, Deepak Kumar Saini, Narendra M Dixit
<br/>
<br/>

There are three folders named `fitness_estimation`, `wright_fisher`, and `genetic_analyses`, wherein each of them contains the files to:
1. `fitness_estimation`: estimate the fitness for various values of γ and for all phenotypes possible.
2. `wright_fisher`: perform the Wright-Fisher evolution simulations for a given mutation rate, using the fitness values estimated in the previous case, to identify the dominating phenotype.
3. `genetic_analyses`: construct the phylogenetic trees for histidine kinases (HKs) and response regulators (RRs) as well as to estimate the ratio of non-synonymous to synonymous mutations (K<sub>A</sub>/K<sub>S</sub>) between the selected neighbors.

Both these folders have two subfolders namely `random_environment` and `programmed_environment`, and as the names identify, the files are separately stored for random environment case, where the signals are elicited in random fashion, and for programmed environment case, where we considered the strict sequence in signals given by _1_, _2_,... _N_, with _N_ being the total number of two-component signaling systems (TCSs). Contents in these folders are described below:

#### 1. `fitness_estimation`
The file `data_set.m` is coded with the differential equations of the model, where `evaluator.m` is the file that is called for various values of _N_, γ, and phenotype for estimation of fitness values. `K_matrix_assignment.m` gives out the interaction matrix for every phenotype (refer to Fig. 2A for the example of _N_ = _2_ case of the manuscript) which is identified by a unique integer. Finally, `bossfile.m` is the main file which can be run to generate fitness data for a given _N_ and γ. Data for the cases _N_ = _2_, _3_, and _4_ are provided in the folder.

The above explained files are for both kinds of environments, but additionally, the random environment case contains an extra `non_path_signal_sequence.m` file to generate all the possible signal sequences, and the fitness of a bacterium in that random environment is the average of fitness values estimated for all possible sequences.

Also, although not used during calculations, the file `K_matrix_inverse.m`, for a given value of _N_, γ, and the matrix structure (refer to Fig. 5B of the manuscript) gives out the unique phenotype ID.

#### 2. `wright_fisher`
The fitness values estimated previously shall be used in the evolution simulations, procedure mentioned in the manuscript. For each kind of environment, as segregated by subfolders, the main file for homogeneous and mixed population initial conditions are separately given as `uniform_evolution.m` and `distributed_evolution.m` respectively. The baseline fitness values are the control fitnesses which decide whether the bacteria, in each generation, die or replicate.

#### 3. `genetic_analyses`
There are four `.txt` documents which contain the amino acid and nucleotide sequences of HKs and RRs in _M. tuberculosis_. The amino acid sequences were aligned using Clustal Omega (`fullHKalignment.fasta` and `fullRRalignment.fasta`), and then used these alignments to construct phylogenetic trees for HKs and RRs separately using the software package MEGA (version 7) and the resulting trees are `HK_full_tree.nwk` and `RR_full_tree.nwk` respectively.

The Excel file `domains_of_interest.xlsx` contains the starting and ending positions of kinase domains of HKs and receiver domains of RRs (refer to the manuscript on how they were identified) in nucleotide and amino acid sequences. Using this information, the alignments are spliced to extract the kinase domain alignment for HK an dreceiver domain alignment for RR respectively from the full protein alignments mentioned previously. The resulting files are `fullHKkinasedomainalignment.fasta` and `fullRRreceiverdomainalignment.fasta` respectively. These domain alignments, which contain amino acids, are converted into respective nucleotide alignments using the nucleotide sequences, and finally the information is used to estimate the K<sub>A</sub>/K<sub>S</sub> ratios for the domains of interest
