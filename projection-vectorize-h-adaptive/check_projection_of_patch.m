load_patch = load('h_quad_lobatto.mat');
load_projection = load('h_projection_test_p2.mat');
load_uniform = load('h_projection_uniform_p2.mat');

h_patch_error = [];
for idxMesh = 1:numel(load_patch.h_mesh_quad_lob_save)
    orimesh = load_patch.h_mesh_quad_lob_save{idxMesh};
    centrecoord=[1000 1000];
    
    p = 2*ones( size(orimesh, 1) - 1, 1);
    [supel1] = initialize_problem(orimesh, centrecoord, p);

    finemesh = refinemesh(orimesh, 1);
    pf = 2*ones( size(finemesh, 1) - 1, 1);
    [supel2] = initialize_problem(finemesh, centrecoord, pf);
    

%     p = ones( size(orimesh, 1)-1, 1);
%     [supel] = initialize_problem(orimesh, centrecoord, p);
%     [local,supel,globnodesdisp] = today3_modify(p, supel);
%     [sigma_star] = nodalStressCalculate3(p, local, supel, 1:5, 1);
%     [sigma_star] = nodalStressCalculate_patch_p2(p, orimesh, local, supel, 1:5, 1);
    % [deltaerror_nosqrt, energyStress, nitaError]= error_calculate3(sigma_star, p, local, supel, globnodesdisp);
%     [deltaerror_nosqrt, energyStress, nitaError]= generate_phi(sigma_star, p, local, supel);
    
    [localcoarse,refineIdx,deltaerror,rhs65,error1,rhs66,errest, ET, cc] = h_project_error_estimator_vectorize(p, supel1, pf, supel2);
    h_patch_error = [h_patch_error; errest]    
    
%     errTotal = sqrt(sum(deltaerror_nosqrt));       
%     sigma_smooth = sqrt(sum(energyStress));
%     error_estimator_patch_projection = [error_estimator_patch_projection; errTotal / sigma_smooth]
%     
end

h_patch_error

%%
dof_projection_mesh = [];
dof_patch_mesh = [];
dof_uni = [];
for idxMesh = 1:numel(load_projection.h_mesh_projection_save)
    dof_projection_mesh = [dof_projection_mesh; size(load_projection.h_mesh_projection_save{idxMesh}, 1)];   
end

for idxMesh = 1:numel(load_patch.h_mesh_quad_lob_save)
    dof_patch_mesh = [dof_patch_mesh; size(load_patch.h_mesh_quad_lob_save{idxMesh}, 1)];
end

for idxMesh = 1:numel(load_uniform.h_mesh_projection_save)
    dof_uni = [dof_uni; size(load_uniform.h_mesh_projection_save{idxMesh}, 1)];
end


loglog(dof_projection_mesh*2, load_projection.error_estimator, 'o-');
hold on
loglog(dof_patch_mesh*2, h_patch_error, '*-');
loglog(dof_uni * 2, load_uniform.error_estimator, '*-');