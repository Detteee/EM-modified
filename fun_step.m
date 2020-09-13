function [y] = fun_step(xx)

    % The function is usually evaluated on the hypercube xi ∈ [-100, 100], for all i = 1, …, d.
  
    
    d = length(xx);
    sum = 0;
    for ii = 1:d
        xi = xx(ii);
        sum = sum + (floor(xi+0.5))^2;
    end
    
    y = sum;
    
    end