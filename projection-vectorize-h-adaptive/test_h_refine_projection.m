function [local, orimesh, supel, error_estimator, computation_time, h_mesh_patch_save] = test_h_refine_projection()
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
local_save = {};
globdisp_save = {};

for iRepeat = 1:5
    tic();

    p = 2*ones( size(orimesh, 1) - 1, 1);
    [supel1] = initialize_problem(orimesh, centrecoord, p);

    finemesh = refinemesh(orimesh, 1);
    pf = 2*ones( size(finemesh, 1) - 1, 1);
    [supel2] = initialize_problem(finemesh, centrecoord, pf);

    [localcoarse,globnodedisp,refineIdx,deltaerror,rhs65,error1,rhs66,errest, ET, cc] = h_project_error_estimator_vectorize(p, supel1, pf, supel2);
    T = toc();
    error_estimator = [error_estimator; errest]    
    computation_time = [computation_time; T]

	h_mesh_projection_save{iRepeat} = orimesh;
    local_save{iRepeat} = localcoarse;
    globdisp_save{iRepeat} = globnodedisp;
    save('h_projection_test_p2.mat', 'h_mesh_projection_save', 'error_estimator', 'local_save', 'globdisp_save');

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
    size(orimesh)
end

end