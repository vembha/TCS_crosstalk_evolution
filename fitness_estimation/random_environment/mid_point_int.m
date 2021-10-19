function [output] = mid_point_int(ar_x,ar_y)
M = 0;
for i=2:length(ar_x)
    M = M + 0.5*(ar_y(i) + ar_y(i-1))*(ar_x(i) - ar_x(i-1));
end
output = M;
end