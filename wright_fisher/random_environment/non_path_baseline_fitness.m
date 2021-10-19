function [] = non_path_baseline_fitness(N,decay_factor)
%% SETTING UP THE INPUT TIMES AND k_deg_input
time_diff = 500;
sim_time = 2*N*time_diff;

if (decay_factor==0)
    k_deg_input = 0;
else
    k_deg_input = 2.3026/decay_factor;
end

toi = [0:(N-1)]*time_diff;
toi(N+1) = sim_time;
%% INITIAL CONDITIONS
HK0 = 100;
RR0 = 1000;
O0 = 100;
I0 = 10000;
Im = I0;

baseline_fitness = zeros(1,power(N,N));
%% LOOPING FOR CALCULATIONS
time = [];
Yfunc = [];
input = [];

for u=1:power(N,N)
    %For each signal sequence...%
    in_seq = non_path_signal_sequence(N,u);
    init_conds = zeros(N,1);
    init_conds(in_seq(1),1) = I0;
    len = 0;
    
    %Running for ODE solver%
    for k=1:N
        [T,Y] = ode15s(@data_set_baseline,[toi(k) toi(k+1)],init_conds,[],N,k_deg_input);
        for i=1:length(T)
            time(i + len) = T(i);
            Yfunc(i + len,:) = Y(i,:);
        end
        
        %Varying initial conditions%
        if (k<N)
            init_conds(1:N,1) = Y(length(T),1:N);
            init_conds(in_seq(k+1),1) = I0;
        end
        len = len + length(T);
    end
    %Fitness evaluation%
    FVs = zeros(N,len);
    for l=1:N
        input(l,1:len) = Yfunc(1:len,l);
        FVs(l,1:len) = exp(-input(l,1:len)/Im);
    end
    FV = prod(FVs);
    baseline_fitness(1,u) = trapz(time(1:len),FV(1:len))/sim_time;
end
%% SAVING THE RESULTS
save(['Non_path_baseline_fitness_N_',num2str(N),'.mat'],'baseline_fitness');
end