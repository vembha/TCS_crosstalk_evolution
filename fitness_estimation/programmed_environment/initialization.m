function [output] = initialization(N,HK0,RR0,I0)
    init_conds = zeros((7*N + 4*N*N),1);
    init_conds(1:N,1) = HK0;
    init_conds((4*N + 1):5*N,1) = RR0;
    init_conds((6*N + 4*N*N + 1),1) = I0;
    output = init_conds;
end