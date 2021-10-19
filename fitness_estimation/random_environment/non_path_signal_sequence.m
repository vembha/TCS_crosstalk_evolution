function [output] = non_path_signal_sequence(N,count)
sequence = zeros(1,N);
temp = count;

for i=1:N
    point = power(N,(N-i));   
    for j=1:N
        if ((temp - j*point)<=0)
            sequence(1,i) = j;
            temp = temp - (j-1)*point;
            break;
        end
        
    end
    
end

output = sequence;
end