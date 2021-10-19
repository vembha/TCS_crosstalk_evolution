function [] = Non_path_uniform_evolution(N,mut,majority_percent)
%% SETTING UP VARIABLES
M = N*(N-1);
P = 2^M;                            % Total number of phenotypes possible
prob_mut = 10^(-mut);               % Mutation rate
K_ratio = KR/100;                   % The value of gamma
majority = majority_percent/100;    % Cut-off percentage of dominant phenotype in the population to stop the simulation

load(['Non_path_fitness_for_N_',num2str(N),'.mat']);
vals = zeros(P,1);
locs = zeros(P,1);
for i=1:P
    [vals(i,1), locs(i,1)] = max(fitness(:,i));
end
load(['Non_path_baseline_fitness_N_',num2str(N),'.mat']);

total_bact = 10000;
gen = 0;
%% ASSIGNMENT FOR FIRST GENERATION
bact_assign = ones(total_bact,1);
count = 1;
fixation(count) = 1;
%% LOOPING FOR EVOLUTION
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
    
    %Varying input signal sequence%
    input_sequence = randperm(power(N,N),1);
    
    %Killing loop%
    for k=1:total_bact
        kill = baseline_fitness(input_sequence) + (1 - baseline_fitness(input_sequence))*rand();
        if (bact_fitness(k)<kill)
            temp(k) = 0;
        end
    end
    
    %Updating for the next generation%
    bact_assign = generation_calculator(KR,temp,total_bact,N,prob_mut);
    fix = mode(bact_assign);
    
    %Saving for the evolution pattern%
    if ~ismember(fix,fixation)
        count = count + 1;
        fixation(count) = fix;
    end
    if(numb(max_index(1,KR),gen)>=majority*total_bact)
        break;
    end
    
    %Printing%
    if(mod(gen,100)==0)
        gen
        fixation
    end
end
%% SAVING THE RESULTS
save(['Non_path_uniform_evolution_N_',num2str(N),'_mut_',num2str(mut),'_maj_',num2str(majority_percent),'.mat']);
end