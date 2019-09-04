function [local, orimesh, supel, error_estimator, computation_time, h_mesh_patch_save] = test_h_refine_projection()
clear;

tempCoord = linspace(1000, 0, 5);
temp1 = tempCoord(2:end);
xElem = [1000*ones(size(tempCoord)), temp1];
yElem = [fliplr(tempCoord), 1000*ones(size(temp1))];
orimesh = [xElem(:), yElem(:)];
orimesh = [[1:1:size(orimesh, 1)]', orimesh];
% p = ones( size(orimesh, 1)-1, 1);
centrecoord=[0 0];
nitaErrorMean = [];
error_estimator = [];
computation_time = [];

h_mesh_projection_save = {};

for iRepeat = 1:10
    tic();

    p = 2*ones( size(orimesh, 1) - 1, 1);
    [supel1] = initialize_problem(orimesh, centrecoord, p);

    finemesh = refinemesh(orimesh, 1);
    pf = 2*ones( size(finemesh, 1) - 1, 1);
    [supel2] = initialize_problem(finemesh, centrecoord, pf);

%     [localcoarse,refineIdx,deltaerror,rhs65,error1,rhs66,errest, ET, cc]= main21_modify(p, supel1, supel2);
    % [localcoarse,refineIdx,deltaerror,rhs65,error1,rhs66,errest, ET, cc]= new_error_estimator_vectorize(p, supel1, supel2);
    [localcoarse,refineIdx,deltaerror,rhs65,error1,rhs66,errest, ET, cc] = h_project_error_estimator_vectorize(p, supel1, pf, supel2);
%     [localcoarse1,refineIdx1,deltaerror1,rhs651,error11,rhs661,errest1, ET1]= new_error_estimator_vectorize(p, supel1, supel2);
    T = toc();
    error_estimator = [error_estimator; errest]    
    computation_time = [computation_time; T]

	h_mesh_projection_save{iRepeat} = orimesh;
    save('h_projection_uniform_p2.mat', 'h_mesh_projection_save', 'error_estimator');

%     if max(refineIdx) == 0
%         break;
%     end
%     newNodeIdx = size(orimesh, 1);
%     for iRefineIdx = numel(refineIdx):-1:1
%         if refineIdx(iRefineIdx) == 1
%             newNodeIdx = newNodeIdx + 1;
%             newNode = [ newNodeIdx, mean(orimesh(iRefineIdx:iRefineIdx+1, 2:3)) ];
%             orimesh = [ orimesh(1:iRefineIdx, :); newNode; orimesh(iRefineIdx+1:end, :)];
%         end
%     end
    orimesh = finemesh;
end

end