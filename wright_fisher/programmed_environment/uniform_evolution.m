function [] = uniform_evolution(N,KR,mut,majority_percent)
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
    
    %Killing loop%
    for k=1:total_bact
        kill = baseline_fitness + (1 - baseline_fitness)*rand();
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
%     if(numb(max_index(1,KR),gen)>=majority*total_bact)
%         break;
%     end
    
    %Printing%
    if(mod(gen,100)==0)
        gen
        fixation
    end
end
%% SAVING THE RESULTS
save(['Uniform_evolution_N_',num2str(N),'_mut_',num2str(mut),'_maj_',num2str(majority_percent),'.mat']);
%% PLOTTING THE GRAPHS
h(1) = figure;
plot(avg_fitness_gen);
title('Average fitness v/s Generation');
xlabel('Generation');
ylabel('Average fitness');

h(2) = figure;
for i=1:count
    plot(1:gen,numb(fixation(i),1:gen)/total_bact);
    legendCell{i} = [num2str(fixation(i),'phenotype = %d')];
    hold on;
end
legend(legendCell,'Location','bestoutside');
title('Fixated phenotypes population fraction v/s generation');
xlabel('Generation');
ylabel('Population fraction');

savefig(h,['Uniform_evolution_N_',num2str(N),'_mut_',num2str(mut),'_maj_',num2str(majority_percent),'.fig']);
close(h);
end