function [localfc] = project_local(pf, localf, supelf, pc, localc, supelc)
% this function will project the localf to localc 
% to find the localfc (local fine coarse)

% degree of freedom of each node
nnel = supelf.nnel;

% each element of localc is supported by 2 elements of localf
% c and lambda will be a combination of both elements
% c and lambda of finecoarse actually the same with the fine, therefore comment those lines
% supelfc.c = supelf.c;
% supelfc.values = supelf.values;

ngl = 16;
[p64,w64] = sbfeglqd1(ngl);

dNj = [];
dNhalf1 = [];
dNhalf2 = [];

for i = 1:ngl
    [N, nn, dnn] = lobatto( p64(i), 2);
    dNj = [dNj, dnn(:)];
    [N, nn, dnn] = lobatto( 0.5 * p64(i) - 0.5, 2);
    dNhalf1 = [dNhalf1; dnn(2)];
    [N, nn, dnn] = lobatto( 0.5 * p64(i) + 0.5, 2);
    dNhalf2 = [dNhalf2; dnn(2)];
end

% k1 = dNj * (dNhalf1 .* w64(:)) / 2;
% k2 = dNj * (dNhalf2 .* w64(:)) / 2;
k1 = dNj * (dNhalf1 .* w64(:));
k2 = dNj * (dNhalf2 .* w64(:));

K = [k1', k2'];
temp = zeros(1, numel(K) * nnel);
temp(1:2:end) = K;
K = [temp; [0 temp(1:end-1)]];

% localfc has the same element as localc so iterate through 
% all elements of localc
for iel = 1:numel(localc)		   
    ielf1 = iel * 2 - 1;
    ielf2 = iel * 2;
	localfc(iel).phi = zeros(size(localf(ielf1).phi));
    
    % project uf1 on uc1 by project phi, the 1st nnel rows
    localfc(iel).phi(1:nnel, :) = localf(ielf1).phi(1:nnel, :);
    % project uf2 on uc2 by project phi
    localfc(iel).phi(end - nnel + 1:end, :) = localf(ielf2).phi(end - nnel + 1:end, :);
    % project middle node of localf to localc
    if pc(iel) == 2
        middle_nodes = 3:4;
        localphi = [localf(ielf1).phi; localf(ielf2).phi];
        localfc(iel).phi(middle_nodes, :) = K * localphi;
    end
end   

end