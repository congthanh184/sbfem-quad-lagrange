close all;
mesh_patch = h_mesh_patch_save{6};
mesh_projection = h_mesh_projection_save{5};

plot(mesh_patch(:, 2), mesh_patch(:, 3), 'x')
hold on
plot(mesh_projection(:, 2), mesh_projection(:, 3), 'o')