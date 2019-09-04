% load_patch = load('h_quad_lobatto_test.mat');
load_patch = load('h_quad_lobatto_p2.mat');
load_projection = load('h_projection_fixed_p2_D_test.mat');
% load_uniform = load('h_projection_uniform_p2.mat');

h_patch_error = [];
for idxMesh = 1:numel(load_patch.h_mesh_quad_lob_save)
    orimesh = load_patch.h_mesh_quad_lob_save{idxMesh};
    centrecoord=[0 0];
    
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

% load_projection.error_estimator(5) = 0.6651;
for idxMesh = 1:4
    dof_projection_mesh = [dof_projection_mesh; size(load_projection.h_mesh_projection_save{idxMesh}, 1)];   
end

for idxMesh = 1:4
    dof_patch_mesh = [dof_patch_mesh; size(load_patch.h_mesh_quad_lob_save{idxMesh}, 1)];
end

% for idxMesh = 1:numel(load_uniform.h_mesh_projection_save)
%     dof_uni = [dof_uni; size(load_uniform.h_mesh_projection_save{idxMesh}, 1)];
% end


semilogx(dof_projection_mesh*2, load_projection.error_estimator(1:4), 'o-');
hold on
semilogx(dof_patch_mesh*2, load_patch.error_estimator(1:4)*100, '*-');
semilogx(dof_patch_mesh(1:3) * 2, h_patch_error(1:3), '*-');