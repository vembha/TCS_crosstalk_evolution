function [dy] = data_set_baseline(t,y,N,k_deg_input)
dy = zeros(N,1);
for i=1:N
    dy(i) = -k_deg_input*y(i);
end
end