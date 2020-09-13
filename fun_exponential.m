function [y] = fun_exponential(xx)

    % The function is usually evaluated on the hypercube xi ∈ [-1.28, 1.28], for all i = 1, …, d.
  
    
    d = length(xx);
    sum = 0;
    for ii = 1:d
        xi = xx(ii);
        sum = sum + xi^2;
    end
    
    y = exp(sum/2)-1;
    
end