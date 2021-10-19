function [fitness] = non_path_evaluator(N,gamma,decay_factor,time_diff,sim_time)
%% SETTING UP THE INPUT TIMES AND k_deg_input
if (decay_factor==0)
    k_deg_input = 0;
else
    k_deg_input = 2.3026/decay_factor;
end

% toi = [0:(N-1)]*time_diff;
% toi(N+1) = sim_time;
toi = 0:time_diff:sim_time;
%% INITIAL CONDITIONS
HK0 = 100;
RR0 = 1000;
O0 = 100;
I0 = 10000;
Im = I0;
K1 = 5e+05;
% K1 = 10;

avg_fitness = zeros(1,2^(N*(N-1)));
%% LOOPING FOR CALCULATIONS
time = [];
Yfunc = [];
output_RRO = [];
input = [];

%For each phenotype...%
for count = 1:2^(N*(N-1))
    KC = K_matrix_assignment(N,gamma,count);
    phen_fitness = zeros(1,power(N,N));
    k_cat_2_matrix = k_cat_2_generator(N,count);
    
    %For each input sequence...%
    for in_seq = 1:power(N,N)
        sequence = non_path_signal_sequence(N,in_seq);  % Function to produce the various signal sequences
        len = 0;
        init_conds = steady_state_values(N,gamma,count,k_cat_2_matrix);
        init_conds((6*N + 4*N*N + sequence(1,1)),1) = I0;
        %For each signal in the sequence...%
        for k=1:N
            %Running for ODE solver%
            options = odeset('RelTol',1e-8,'AbsTol',1e-8);
            [T,Y] = ode15s(@data_set,[toi(k) toi(k+1)],init_conds,options,N,KC,k_deg_input,k_cat_2_matrix);
            for i=1:length(T)
                time(i + len) = T(i);
                Yfunc(i + len,:) = Y(i,:);
            end
            for l=1:N
                RR(l,len + 1:len + length(T)) = Yfunc(len + 1:len + length(T),5*N + l); % RR-P for each TCS
                input(l,len + 1:len + length(T)) = Yfunc(len + 1:len + length(T),6*N + 4*N*N + l);  % Input signal dynamics
            end
            
            %Varying initial conditions to add the next signal in sequence%
            if (k<N)
                init_conds(1:(7*N + 4*N*N),1) = Y(length(T),1:(7*N + 4*N*N));
                for q = 1:1:(N-1)
                    init_conds(6*N + 4*N*N + q + 1) = 0;
                end
                init_conds(6*N + 4*N*N + sequence(k + 1),1) = I0;
            end
            len = len + length(T);
        end
        
        %Fitness estimation%
        FVs = zeros(N,len);
        for l=1:N
            FVs(l,1:len) = exp((K1./(K1 + (RR(l,1:len)).^2)).*(-input(l,1:len)/Im));
        end
        FV = prod(FVs);
        phen_fitness(1,in_seq) = mid_point_int(time(1:len),FV(1:len))/sim_time;
        clearvars time Yfunc input;
    end
    avg_fitness(1,count) = mean(phen_fitness);
%     count
end
%% RETURNING THE FITNESS
fitness = avg_fitness;
end