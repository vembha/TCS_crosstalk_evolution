function [next_generation] = generation_calculator(KR,temp,total_bact,N,prob_mut)
%% REMOVING THE KILLED BACTERIA
clear new_bact_assign;
temp = sort(temp,'descend');
count = length(temp) - sum(temp(:)==0);
temp = temp(1:count);
%% RECPLICATION OF BACTERIA ALIVE
d = zeros(2,count);
d(1,:)=temp(:);
d(2,:)=temp(:);
daughter=d(:);
%% FILLING EMPTY SPOTS/REMOVING EXTRA ONES
if (length(daughter)>total_bact)
    X_len = total_bact - length(d(1,:));
    daughter_array = d(1,:);
    rand_shuff = daughter_array(randperm(length(daughter_array)));
    new_bact_assign(1:count) = d(1,:);
    new_bact_assign(count + 1:total_bact) = rand_shuff(1:X_len);
else
    new_bact_assign(1:length(daughter)) = daughter(1:length(daughter));
    new_bact_assign(length(daughter) + 1:total_bact) = randperm(2^(N*(N-1)),1);
end
%% ESTIMATING FOR MUTATIONS
for i=1:total_bact
    KC_temp = K_matrix_assignment(N,KR/100,new_bact_assign(i));
    array = zeros(N*(N-1),1);
    count = 1;
    for m=1:N
        for n=1:N
            if (m~=n)
                array(count) = KC_temp(m,n);
                count = count + 1;
            end
        end
    end
    
    count = 1;
    for j=1:N*(N-1)
        R1 = rand();
        if (R1<prob_mut)
            if (array(count)==0)
                array(count) = KR/100;
            elseif (array(count~=0))
                array(count) = 0;
            end
        end
        count = count + 1;
    end
    
    KC_mutated = ones(N,N);
    count = 1;
    for m=1:N
        for n=1:N
            if(m~=n)
                KC_mutated(m,n) = array(count);
                count = count + 1;
            end
        end
    end
    KC_final = KC_mutated.*0.001;
    
    new_bact_assign(i) = K_matrix_inverse(N,KR/100,KC_final);
end
%% RETURNING THE NEXT GENERATION ARRAY
next_generation = new_bact_assign;
end