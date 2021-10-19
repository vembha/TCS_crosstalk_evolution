function [output] = K_matrix_assignment(N,K_ratio,count)
K_matrix = zeros(N,N);
for i=1:N
    K_matrix(i,i) = 0.001;
end

array = zeros(1,N*(N-1));

for i=1:N*(N-1)
    unit = 2^(i-1);
    nos = 1;
    diff = count - unit;
    while(diff>0)
        diff = diff - unit;
        nos = nos + 1;
    end
    
    if (mod(nos,2)==0)
        array(1,i) = K_ratio;
    end
end
array = array.*0.001;

record = 1;
for i=1:N
    for j=1:N
        if (i~=j)
            K_matrix(i,j) = array(1,record);
            record = record + 1;
        end
    end
end

output = K_matrix;
end