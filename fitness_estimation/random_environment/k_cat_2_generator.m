function [output] = k_cat_2_generator(N,count)
temp_matrix = K_matrix_assignment(N,1,count);
temp_matrix = temp_matrix.*(10^3);
for i = 1:1:N
    for j = 1:1:N
        if (i == j)
            temp_matrix(i,j) = 0.01;
        else
            if (temp_matrix(i,j) ~= 0)
                temp_matrix(i,j) = 1/600;
            end
        end
    end
end
output = temp_matrix;
end