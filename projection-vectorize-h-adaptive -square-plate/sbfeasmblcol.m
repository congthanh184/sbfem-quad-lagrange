function [kk]= sbfeasmbl(kk,k,index)
%----------------------------------------------------------
% Purpose
%   Assembly of element matrices into the system matrix
%
% Synopsis
%   [kk]= sbfeasmbl(kk,k,index)
% Variable description
%   kk - system matrix, divided into 4 squares each of size sdof*sdof,
%   and with order [1 3
%                   2 4]
%   k - element matrix
%   index - dofs vector associated with an element
%-----------------------------------------------------------

edof= length(index);

for i=1:edof
    ii= index(i);
    
    
                
        kk(ii)= kk(ii)+k(i);
        
    
end
