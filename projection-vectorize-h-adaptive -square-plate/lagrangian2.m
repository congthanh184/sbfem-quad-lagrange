function [NN] = lagrangian2(n)
    NN = [  -0.5*n*(1-n), ...
            (1+n)*(1-n), ...
            0.5*n*(1+n) ];
end