function [output] = steady_state_values(N,gamma,phenotype,k_cat_2_matrix)
decay_factor = 500;
% if (decay_factor==0)
%     k_deg_input = 0;
% else
%     k_deg_input = 2.3026/decay_factor;
% end
k_deg_input = 0;
sim_time = 2000000;
%% INITIAL CONDITIONS
HK0 = 100;
RR0 = 1000;
% O0 = 100;
I0 = 0;
% Im = 10000;
% K1 = 5e+05;
%% SOLVING FOR DIFFERENT PHENOTYPES
KC = K_matrix_assignment(N,gamma,phenotype);
init_conds = initialization(N,HK0,RR0,I0);
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[T,Y] = ode15s(@data_set,[0 sim_time],init_conds,options,N,KC,k_deg_input,k_cat_2_matrix);
len = length(T);
output = Y(len,:);
output = output';
end