function [dy] = data_set(t,y,N,KC,k_deg_input,k_cat_2_matrix)
dy = zeros((7*N + 4*N*N),1);
%% PARAMETERS
kfi = 1e-10;
k_prim_fi = 0.1;
kbi = 0.1;
k_prim_bi = 0.1;
k_plus_i = 5;
k_prim_plus_i = 5;
k_minus_i = 1000;
k_prim_minus_i = 1e-06;
k_cat_1ij = 0.055;
k_cat_2ij = k_cat_2_matrix;
k_plus_1 = KC;
k_minus_1ij = 1;
k_plus_2 = KC;
k_minus_2ij = 1;
% k_p_bind = 0.01;
% k_p_unbind = 0.1;
k_deg = 6e-05;
alpha = 10;
beta = 6e-03;
lambda = 0.1;
K1 = 5e+05;
% K1 = 10;
P_T = 100;
%% EQUATIONS
sum_HK = zeros(N,1);
sum_HKstar = zeros(N,1);
sum_IHK = zeros(N,1);
sum_IHKstar = zeros(N,1);
sum_RR = zeros(N,1);
sum_RRstar = zeros(N,1);

%For HK%
for i=1:N
    sum_HK(i) = 0;
    for j=1:N
        sum_HK(i) = sum_HK(i) - (k_plus_2(i,j)*y(i)*y(5*N + j)) + (k_minus_2ij*y(6*N + N*N + N*(i-1) + j) + k_cat_1ij*y(6*N + N*(i-1) + j) + k_cat_2ij(i,j)*y(6*N + N*N + N*(i-1) + j));
    end
    dy(i) = -(kfi*y(i) + k_plus_i*y(i)*y(6*N + 4*N*N + i)) + (kbi*y(N + i) + k_minus_i*y(2*N + i)) + lambda*beta*P_T*(1 + alpha*y(5*N + i)*y(5*N + i)/K1)/(1 + y(5*N + i)*y(5*N + i)/K1) - k_deg*y(i) + sum_HK(i);
end

%For HK*%
for i=1:N
    sum_HKstar(i) = 0;
    for j=1:N
        sum_HKstar(i) = sum_HKstar(i) - (k_plus_1(i,j)*y(N + i)*y(4*N + j)) + (k_minus_1ij*y(6*N + N*(i-1) +j));
    end
    dy(N + i) = -(kbi*y(N + i) + k_prim_plus_i*y(6*N + 4*N*N + i)*y(N + i)) + (kfi*y(i) + k_prim_minus_i*y(3*N + i)) - k_deg*y(N + i) + sum_HKstar(i);
end

%For IHK%
for i=1:N
    sum_IHK(i) = 0;
    for j=1:N
        sum_IHK(i) = sum_IHK(i) - (k_plus_2(i,j)*y(2*N + i)*y(5*N + j)) + (k_minus_2ij*y(6*N + 3*N*N + N*(i-1) + j) + k_cat_1ij*y(6*N + 2*N*N + N*(i-1) +j) + k_cat_2ij(i,j)*y(6*N + 3*N*N + N*(i-1) + j));
    end
    dy(2*N + i) = -(k_minus_i*y(2*N + i) + k_prim_fi*y(2*N + i)) + (k_plus_i*y(6*N + 4*N*N + i)*y(i) + k_prim_bi*y(3*N + i)) - k_deg*y(2*N + i) + sum_IHK(i);
end

%For IHK*%
for i=1:N
    sum_IHKstar(i) = 0;
    for j=1:N
        sum_IHKstar(i) = sum_IHKstar(i) - (k_plus_1(i,j)*y(3*N + i)*y(4*N + j)) + (k_minus_1ij*y(6*N + 2*N*N + N*(i-1) +j));
    end
    dy(3*N + i) = -(k_prim_minus_i*y(3*N + i) + k_prim_bi*y(3*N + i)) + (k_prim_plus_i*y(6*N + 4*N*N + i)*y(N + i) + k_prim_fi*y(2*N + i)) - k_deg*y(3*N + i) + sum_IHKstar(i);
end

%For RR%
for j=1:N
    sum_RR(j) = 0;
    for i=1:N
        sum_RR(j) = sum_RR(j) - (k_plus_1(i,j)*y(N + i)*y(4*N + j) + k_plus_1(i,j)*y(3*N + i)*y(4*N + j)) + (k_minus_1ij*y(6*N + j + (i-1)*N) + k_minus_1ij*y(6*N + 2*N*N + j + (i-1)*N) + k_cat_2ij(i,j)*y(6*N + N*N + j + (i-1)*N) + k_cat_2ij(i,j)*y(6*N + 3*N*N + j +(i-1)*N));
    end
    dy(4*N + j) = beta*P_T*(1 + alpha*y(5*N + j)*y(5*N + j)/K1)/(1 + y(5*N + j)*y(5*N + j)/K1) - k_deg*y(4*N + j) + sum_RR(j);
end

%For RR*%
for j=1:N
    sum_RRstar(j) = 0;
    for i=1:N
        sum_RRstar(j) = sum_RRstar(j) - (k_plus_2(i,j)*y(i)*y(5*N + j) + k_plus_2(i,j)*y(2*N + i)*y(5*N + j)) + (k_cat_1ij*y(6*N + j + (i-1)*N) + k_cat_1ij*y(6*N + 2*N*N + j + (i-1)*N) + k_minus_2ij*y(6*N + N*N + j + (i-1)*N) + k_minus_2ij*y(6*N + 3*N*N + j + (i-1)*N));
    end
    dy(5*N + j) = - k_deg*y(5*N + j) + sum_RRstar(j);
end

%For HK*RR%
for i=1:N
    for j=1:N
        dy(6*N + N*(i-1) + j) = -(k_minus_1ij*y(6*N + N*(i-1) + j) + k_cat_1ij*y(6*N + N*(i-1) + j) + k_prim_plus_i*y(6*N + 4*N*N + i)*y(6*N + N*(i-1) + j)) + (k_plus_1(i,j)*y(N + i)*y(5*N + j) + k_prim_minus_i*y(6*N + 2*N*N + N*(i-1) + j)) - k_deg*y(6*N + N*(i-1) + j);
    end
end

%For HKRR*%
for i=1:N
    for j=1:N
        dy(6*N + N*N + N*(i-1) + j) = -(k_minus_2ij*y(6*N + N*N + N*(i-1) + j) + k_cat_2ij(i,j)*y(6*N + N*N + N*(i-1) + j) + k_plus_i*y(6*N + 4*N*N + i)*y(6*N + N*N + N*(i-1) + j)) + (k_plus_2(i,j)*y(i)*y(5*N + j) + k_minus_i*y(6*N + 3*N*N + N*(i-1) + j)) - k_deg*y(6*N + N*N + N*(i-1) + j);
    end
end

%For IHK*RR%
for i=1:N
    for j=1:N
        dy(6*N + 2*N*N + N*(i-1) + j) = -(k_minus_1ij*y(6*N + 2*N*N + N*(i-1) + j) + k_cat_1ij*y(6*N + 2*N*N + N*(i-1) + j) + k_prim_minus_i*y(6*N + 2*N*N + N*(i-1) + j)) + (k_plus_1(i,j)*y(3*N + i)*y(4*N + j) + k_prim_plus_i*y(6*N + 4*N*N + i)*y(6*N + N*(i-1) + j)) - k_deg*y(6*N + 2*N*N + N*(i-1) + j);
    end
end

%For IHKRR*%
for i=1:N
    for j=1:N
        dy(6*N + 3*N*N + N*(i-1) + j) = -(k_minus_2ij*y(6*N + 3*N*N + N*(i-1) + j) + k_cat_2ij(i,j)*y(6*N + 3*N*N + N*(i-1) + j) + k_minus_i*y(6*N + 3*N*N + N*(i-1) + j)) + (k_plus_2(i,j)*y(2*N + i)*y(5*N + j) + k_plus_i*y(6*N + 4*N*N + i)*y(6*N + N*N + N*(i-1) + j)) - k_deg*y(6*N + 3*N*N + N*(i-1) + j);
    end
end

%For I%
for i=1:N
    dy(6*N + 4*N*N + i) = -k_deg_input*y(6*N + 4*N*N + i);
end
end