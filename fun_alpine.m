function [y] = fun_alpine(xx)

    % The function is usually evaluated on the hypercube xi ∈ [-10, 10], for all i = 1, …, d.
  
    
    d = length(xx);
    sum = 0;
    for ii = 1:d
        xi = xx(ii);
        sum = sum + abs(xi*sin(xi)+0.1*xi);
    end
    
    y = sum;
    
    end