function [local, orimesh, supel, error_estimator, computation_time, h_mesh_uniform_save] = main_h_adaptive()
clear;

tempCoord = linspace(1000, 0, 5);
temp1 = tempCoord(2:end);
xElem = [1000*ones(size(tempCoord)), temp1];
yElem = [fliplr(tempCoord), 1000*ones(size(temp1))];
orimesh = [xElem(:), yElem(:)];
orimesh = [[1:1:size(orimesh, 1)]', orimesh];
% p = 2 * ones( size(orimesh, 1)-1, 1);
centrecoord=[0 0];

nitaErrorMean = [];
error_estimator = [];
computation_time = [];

h_mesh_uniform_save = {};
local_save = {};
globdisp_save = {};

for iRepeat = 1:4
    tic();
    p = 2*ones( size(orimesh, 1)-1, 1);
    [supel] = initialize_problem(orimesh, centrecoord, p);

    [local,supel,globnodesdisp] = today3_modify(p, supel);
    [sigma_star] = nodalStressCalculate_patch_p2(p, orimesh, local, supel,[1,5,9], 2);
    [deltaerror_nosqrt, energyStress, nitaError]= generate_phi(sigma_star, p, local, supel);

    nitaError;
    meanOfNitaErr = mean(nitaError);
    maxOfNitaErr = max(nitaError);
    nitaErrorMean = [nitaErrorMean, sum(nitaError)*100];
    errTotal = sqrt(sum(deltaerror_nosqrt));       
    sigma_smooth = sqrt(sum(energyStress));
    error_estimator = [error_estimator; (errTotal / sigma_smooth) * 100]
        
    refineIdx = 100*nitaError >= 0.2;
%     refineIdx = nitaError > (1/3)*maxOfNitaErr;
    
    T = toc(); % in minute
    computation_time = [computation_time; T]
    h_mesh_quad_lob_save{iRepeat} = orimesh;
    local_save{iRepeat} = local;
    globdisp_save{iRepeat} = globnodesdisp;
    save('h_quad_lobatto.mat', 'h_mesh_quad_lob_save', 'error_estimator', 'local_save', 'globdisp_save');

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
%     orimesh = refinemesh(orimesh, 1);
    size(orimesh)
end

end