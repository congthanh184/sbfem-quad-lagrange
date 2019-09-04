function [local, orimesh, supel, error_estimator, computation_time, h_mesh_patch_save] = test_h_projection()
clear;

orimesh = [ 1 1000 2000;...
    6 500 2000; ...
    2 0 2000;...
    7 0 1000;...
    3 0 0;...
    8 1000 0; ...
    4 2000 0;...
    9 2000 500; ...
    5 2000 1000];
centrecoord=[1000 1000];
nitaErrorMean = [];
error_estimator = [];
computation_time = [];

h_mesh_projection_save = {};

for iRepeat = 1:10
    tic();
    p = ones( size(orimesh, 1)-1, 1);
%     [supel] = initialize_problem(orimesh, centrecoord, p);
    [supel1] = initialize_problem(orimesh, centrecoord, p);
    [supel2] = initialize_problem(orimesh, centrecoord, p+ones(size(p)));

%     [localcoarse,refineIdx,deltaerror,rhs65,error1,rhs66,errest, ET, cc]= main21_modify(p, supel1, supel2);
    [localcoarse,refineIdx,deltaerror,rhs65,error1,rhs66,errest, ET, cc]= new_error_estimator_vectorize(p, supel1, supel2);
%     [localcoarse1,refineIdx1,deltaerror1,rhs651,error11,rhs661,errest1, ET1]= new_error_estimator_vectorize(p, supel1, supel2);
    T = toc();
    error_estimator = [error_estimator; errest]    
    computation_time = [computation_time; T]

	h_mesh_projection_save{iRepeat} = orimesh;
    save('h_projection_result.mat', 'h_mesh_projection_save', 'error_estimator');

    if max(refineIdx) == 0
        break;
    end
    newNodeIdx = size(orimesh, 1);
    for iRefineIdx = numel(refineIdx):-1:1
        if refineIdx(iRefineIdx) == 1
            newNodeIdx = newNodeIdx + 1;
            newNode = [ newNodeIdx, mean(orimesh(iRefineIdx:iRefineIdx+1, 2:3)) ];
            orimesh = [ orimesh(1:iRefineIdx, :); newNode; orimesh(iRefineIdx+1:end, :)];
        end
    end
    orimesh
end

end