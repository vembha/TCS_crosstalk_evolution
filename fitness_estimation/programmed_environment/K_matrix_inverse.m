function [output] = K_matrix_inverse(N,K_ratio,KC)
array = zeros(1,(N*(N-1)));
count = 1;
for m=1:N
    for n=1:N
        if (m~=n)
            array(count) = KC(m,n);
            count = count + 1;
        end
    end
end

for i=1:(N*(N-1))
    if (array(i)~=0)
        array(i) = 1;
    end
end

len = 1;
for i=1:(N*(N-1))
    len = len + (2^(i-1))*(array(i));
end

output = len;
end