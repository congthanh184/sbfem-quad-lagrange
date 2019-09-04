clear;
% tempCoord = linspace(1000, 0, 5);
% temp1 = tempCoord(2:end);
% horiElem = [fliplr(tempCoord); 1000*ones(size(tempCoord))];
% vertElem = [1000*ones(size(temp1)); temp1];
% orimesh = [horiElem'; vertElem'];
% orimesh = [[1:1:size(orimesh, 1)]', orimesh];
% centrecoord=[0 0];
% p = ones( size(orimesh, 1)-1, 1);

% tempCoord = linspace(1000, 0, 5);
% temp1 = tempCoord(2:end);
% xElem = [1000*ones(size(tempCoord)), temp1];
% yElem = [fliplr(tempCoord), 1000*ones(size(temp1))];
% orimesh = [xElem(:), yElem(:)];
% orimesh = [[1:1:size(orimesh, 1)]', orimesh];
% centrecoord=[0 0];

orimesh = [ 1 1000 2000;...
    6 500 2000; ...
    2 0 2000;...
    9 0 1500; ...
    7 0 1000; ...
    10 0 500; ...
    3 0 0;...
    11 500 0; ...
    8 1000 0; ...
    12 1500 0; ...
    4 2000 0;...
%     9 2000 500; ...
    5 2000 1000];
centrecoord=[1000 1000];
p = 2*ones( size(orimesh, 1)-1, 1);

% tempCoord = linspace(1000, 0, 5);
% temp1 = tempCoord(1:end);
% xElem = [tempCoord, 0, 1000*ones(size(temp1))];
% yElem = [1000*ones(size(tempCoord)), 0, fliplr(temp1)];
% orimesh = [xElem', yElem'];
% orimesh = [[1:1:size(orimesh, 1)]', orimesh]
% centrecoord=[500 500];
% p = ones( size(orimesh, 1)-1, 1);

% tempCoord = linspace(1000, 0, 5);
% temp1 = tempCoord(2:end);
% xElem = [(tempCoord)'; 0; 1000];
% yElem = [1000*ones(size(tempCoord))'; 0; 0];
% orimesh = [xElem, yElem];
% orimesh = [[1:1:size(orimesh, 1)]', orimesh];
% p = 2*ones( size(orimesh, 1)-1, 1);

% centrecoord=[1000 500];

% distributed force between node 1 & 2
% constrain x on node1, constrain y on node 9
[supel] = initialize_problem(orimesh, centrecoord, p);
[local,supel,globnodesdisp] = today3_modify(p, supel);
% [sigma_star] = nodalStressCalculate3(p, local, supel, 1:5, 1);
[sigma_star] = nodalStressCalculate_patch_p2(p, orimesh, local, supel, 1:5, 2);
% [recStress] = recoverFromModalStress(sigma_star, supel, local, p);
% [recStress] = recoverFromModalStress_p2(sigma_star, supel, local, p);
[recStress] = recoverFromModalStress_p2_vectorize(sigma_star, supel, local, p);
plot_stress(local, recStress);
% [deltaerror, energyStress, nitaError]= error_calculate3(sigma_star, p, local, supel, globnodesdisp)
% plot_disp(local);
% test_recovery_patch
% plot4_stress(local, local)