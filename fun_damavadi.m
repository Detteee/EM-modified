function y = fun_damavadi(x)
    x1 = x(1);x2 = x(2);
    y1 = pi*(x1-2);y2 = pi*(x2-2);
    y = (1-abs(sin(y1)*sin(y2)/y1/y2)^5)*(2+(x1-7)^2+2*(x2-7)^2);
end