function [b, index] = CPWL_function(a, aVect, mu0, mu1, etaVect)

    b = mu0 + mu1*a + etaVect*abs(a-aVect)'; % + nuVect'*sign(a_old-aVect);
    
    % Calculation needed to find the new slope for adaptation condition
    for j = 1:length(aVect)-1
        if min(aVect(j+1), aVect(j)) <= a && a <= max(aVect(j+1), aVect(j))
            index = j;
            break
        elseif a <= aVect(1)
            index = 1;
            break
        elseif a >= aVect(end)
            index = length(aVect)-1;
            break
        end
    end
    
end

