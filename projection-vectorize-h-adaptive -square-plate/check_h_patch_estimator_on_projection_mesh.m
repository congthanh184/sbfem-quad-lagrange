error_estimator_patch_projection = [];

for idxMesh = 1:numel(h_mesh_projection_save)
    orimesh = h_mesh_projection_save{idxMesh};

    centrecoord=[1000 1000];
    p = ones( size(orimesh, 1)-1, 1);
    [supel] = initialize_problem(orimesh, centrecoord, p);
    [local,supel,globnodesdisp] = today3_modify(p, supel);
    [sigma_star] = nodalStressCalculate3(p, local, supel, 1:5, 1);
%     [sigma_star] = nodalStressCalculate_patch_p2(p, orimesh, local, supel, 1:5, 1);
    % [deltaerror_nosqrt, energyStress, nitaError]= error_calculate3(sigma_star, p, local, supel, globnodesdisp);
    [deltaerror_nosqrt, energyStress, nitaError]= generate_phi(sigma_star, p, local, supel);

    errTotal = sqrt(sum(deltaerror_nosqrt));       
    sigma_smooth = sqrt(sum(energyStress));
    error_estimator_patch_projection = [error_estimator_patch_projection; errTotal / sigma_smooth]
    
end

error_estimator_patch_projection

%%
dof_projection_mesh = [];
dof_patch_mesh = [];
for idxMesh = 1:numel(h_mesh_projection_save)
    dof_projection_mesh = [dof_projection_mesh; size(h_mesh_projection_save{idxMesh}, 1)];   
end

for idxMesh = 1:numel(h_mesh_patch_save)
    dof_patch_mesh = [dof_patch_mesh; size(h_mesh_patch_save{idxMesh}, 1)];
end

loglog(dof_projection_mesh*2, error_estimator_patch_projection*100, 'o-');
hold on
loglog(dof_patch_mesh*2, error_estimator*100, '*-');