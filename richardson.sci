function [xkp1] = richardson (A, b)
    vp = spec(A);
    alpha = 2/(min(vp) + max(vp));
    
    xk = 0;
    xkp1 = xk + alpha*(b-A*xk);
    iter = 0;
    while (abs(xkp1 - xk) < 10e-3 & iter < 1000) :
        xk = xkp1;
        xkp1 = xk + alpha*(b-A*xk);
        iter = iter + 1;
    end
endfunction
