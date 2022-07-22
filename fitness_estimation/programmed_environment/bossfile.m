%% |BOSSFILE| |FITNESS OF ALL PHENOTYPES & MAXIMUM FITNESS WITH gamma|
%% SETTING-UP THE VALUES
N = 2;                      % Number of TCSs
decay_factor = 0;           % For square pulse signals
% decay_factor = 500;         % Time for input signal to decay by 1000 fold
time_diff = 500;            % Time gap between input signals
sim_time = N*time_diff;     % Total simulation time
fitness = zeros(50,2^(N*(N-1)));
num_gamma = 50;             % Number of values of gamma for which the estimation needs to be carried
%% FITNESS EVALUATION OF FOR ALL POSSIBLE PHENOTYPES
% tic
for i=1:num_gamma
    gamma = i/100;
    fitness(i,:) = evaluator(N,gamma,decay_factor,time_diff,sim_time);
end
% toc
%% ESTIMATING THE PHENOTYPE WITH MAXIMUM FITNESS
for i=1:num_gamma
    [max_fit, max_index] = max(fitness');
end
%% SAVING THE REQUIRED VARIABLES
% The matrix 'fitness' has fitness data for each gamma value in rows and
% each phenotype in columns
save(['Fitness_for_N_',num2str(N),'.mat'],'fitness','max_fit','max_index');
