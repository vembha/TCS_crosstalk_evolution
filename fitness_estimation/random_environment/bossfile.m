%% |BOSSFILE| |FITNESS OF ALL PHENOTYPES & MAXIMUM FITNESS WITH gamma|
%% SETTING-UP THE VALUES
N = 2;                      % Number of TCSs
decay_factor = 500;         % Time for input signal to decay by 1000 fold
time_diff = 500;            % Time gap between input signals
sim_time = 2*N*time_diff;   % Total simulation time
fitness = zeros(20,2^(N*(N-1)));
%% FITNESS EVALUATION OF FOR ALL POSSIBLE PHENOTYPES
for i=1:20
    gamma = i/100;
    fitness(i,:) = non_path_evaluator(N,gamma,decay_factor,time_diff,sim_time);
    i
end
delete(gcp);
%% ESTIMATING THE PHENOTYPE WITH MAXIMUM FITNESS
for i=1:20
    [max_fit, max_index] = max(fitness');
end
%% SAVING THE REQUIRED VARIABLES
% The matrix 'fintess' contains fitness values estimated for each phenotype
% in columns and for each value of gamma in rows
save(['Non_path_fitness_for_N_',num2str(N),'.mat'],'fitness','max_fit','max_index');