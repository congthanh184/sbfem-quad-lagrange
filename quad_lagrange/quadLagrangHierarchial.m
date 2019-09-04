function [N, NN, dN, dNN] = quadLagrangHierarchial(n, lb, k)
    NN = [  -0.5*n*(1-n), ...
	        (1+n)*(1-n), ...
	        0.5*n*(1+n) ];
	dNN = [ -0.5+n, ...
            -2*n, ...
			0.5+n ];
    
    newNN = [];
    newdNN = [];
    ab = [];
    for ilam = 1:lb
        for ik = 1:2^ilam
            ank = -1 + (ik-1)/2^(ilam-1);
            bnk = ank + 1/2^(ilam-1);
            if n >= ank && n <= bnk
                t = -1+2^(ilam)*(n - ank);
                newNN = [newNN 1-t^2];
                newdNN = [newdNN -2*t];
                tt = [ank; bnk];
                ab = [ ab tt];
            else
                newNN = [newNN 0];
                newdNN = [newdNN 0];
            end
        end
    end
    temp = NN;
    NN([2 4]) = newNN;
    NN([1 3 5]) = temp;

    temp = dNN;
    dNN([2 4]) = newdNN;
    dNN([1 3 5]) = temp;

    N = (convertSingleArrayToEyeMatrix(NN, 2))';
    dN = (convertSingleArrayToEyeMatrix(dNN, 2))';
end