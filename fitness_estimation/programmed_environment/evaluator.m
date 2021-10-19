function [fitness] = evaluator(N,gamma,decay_factor,time_diff,sim_time)
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

avg_fitness = zeros(1,2^(N*(N-1)));
%% LOOPING FOR CALCULATIONS
for count = 1:2^(N*(N-1))
    KC = K_matrix_assignment(N,gamma,count);        % Interaction matrix
    k_cat_2_matrix = k_cat_2_generator(N,count);    % Phosphatase rate matrix
    init_conds = steady_state_values(N,gamma,count,k_cat_2_matrix);
    init_conds(6*N + 4*N*N + 1,1) = I0;
    len = 0;
    for k=1:N
        %Running for ODE solver%
        options = odeset('RelTol',1e-8,'AbsTol',1e-8);
        [T,Y] = ode15s(@data_set,[toi(k) toi(k+1)],init_conds,options,N,KC,k_deg_input,k_cat_2_matrix);
        for i=1:length(T)
            time(i + len) = T(i);
            Yfunc(i + len,:) = Y(i,:);
        end
        for l=1:N
            RR(l,len + 1:len + length(T)) = Yfunc(len + 1:len + length(T),5*N + l); % RR-P profile for each TCS
            input(l,len + 1:len + length(T)) = Yfunc(len + 1:len + length(T),6*N + 4*N*N + l);  % Input signal profile stimulated TCS
        end
        
        %Varying initial conditions to add for subsequent signal%
        if (k<N)
            init_conds(1:(7*N + 4*N*N),1) = Y(length(T),1:(7*N + 4*N*N));
            for q = 1:1:(N-1)
                init_conds(6*N + 4*N*N + q + 1,1) = 0;
            end
            init_conds(6*N + 4*N*N + k + 1,1) = I0;
        end
        len = len + length(T);
    end
    
    %Fitness evaluation%
    FVs = zeros(N,len);
    for l=1:N
        FVs(l,1:len) = exp((K1./(K1 + (RR(l,1:len)).^2)).*(-input(l,1:len)/Im));    % 'phi' for each TCS
    end
    FV = prod(FVs);
    avg_fitness(1,count) = mid_point_int(time(1:len),FV(1:len))/sim_time;   % '<phi>' value
    clearvars time Yfunc input RR;
end
%% RETURNING THE FITNESS
fitness = avg_fitness;
end