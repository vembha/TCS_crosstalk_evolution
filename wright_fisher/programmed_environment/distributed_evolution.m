function [] = distributed_evolution(N,KR,mut,majority_percent)
%% SETTING UP VARIABLES
M = N*(N-1);
P = 2^M;                            % Total number of phenotypes possible
prob_mut = 10^(-mut);               % Mutation rate
K_ratio = KR/100;                   % The value of gamma
majority = majority_percent/100;    % Cut-off percentage of dominant phenotype in the population to stop the simulation

load(['Fitness_for_N_',num2str(N),'.mat']);
vals = zeros(P,1);
locs = zeros(P,1);
for i=1:P
    [vals(i,1), locs(i,1)] = max(fitness(:,i));
end
load(['Baseline_fitness_N_',num2str(N),'.mat']);

npp = 100;
total_bact = P*npp;
gen = 0;
%% ASSIGN FITNESSES FOR THE FIRST GENERATION
temp = zeros(npp,P);
for i=1:P
    temp(:,i) = i;
end
bact_assign(:,1) = temp(:);
%% LOOPING FOR EVOLUTION CALCULATIONS
while (gen<=10000)
    %Assigning fitnesses for the population%
    gen = gen + 1;
    bact_fitness = zeros(total_bact,1);
    bact_fitness(1:total_bact,1) = vals(bact_assign(1:total_bact),1);
    
    %Copying for computations%
    temp = bact_assign;
    
    %Counting the number of phenotypes of each kind%
    for j=1:P
        numb(j,gen) = sum(bact_assign(1:total_bact)==j);
    end
    avg_fitness_gen(gen) = mean(bact_fitness);
    
    %Killing loop%
    for k=1:total_bact
        kill = baseline_fitness + (1 - baseline_fitness)*rand();
        if (bact_fitness(k)<kill)
            temp(k) = 0;
        end
    end
    
    %Updating for the next generation%
    bact_assign = generation_calculator(KR,temp,total_bact,N,prob_mut);
    fixation = mode(bact_assign);
    
    %Printing%
    if(mod(gen,100)==0)
        gen
        fixation
    end
end
%% FIXATED PHENOTYPE CROSS-TALK PATTERN
KC_fix = K_matrix_assignment(N,K_ratio,fixation);
%% SAVING THE RESULTS
save(['Distributed_evolution_N_',num2str(N),'_mut_',num2str(mut),'_maj_',num2str(majority_percent),'.mat']);
end