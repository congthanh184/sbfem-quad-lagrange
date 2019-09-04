for i = 1:6
    orimesh = h_mesh_projection_save{i};
    plot(orimesh(:, 2), orimesh(:, 3), 'rx');
    axis([0 1200 0 1200]);
%     hold on
    pause
end

%%
for i = 1:5
    orimesh = h_mesh_projection_save{i};
    x = zeros(size(orimesh(:,2), 1), 1);
    y = zeros(size(orimesh(:,3), 1), 1);
    
    x = [x, orimesh(:, 2)]';
    y = [y, orimesh(:, 3)]';
    plot(orimesh(:, 2), orimesh(:, 3), 'r-', 'LineWidth', 2);
    hold on
    plot(x, y, 'bo-', 'LineWidth',2)
    axis([0 1200 0 1200]);
    hold off
    pause
end

%%
for i = 2:2:5
    orimesh = h_mesh_projection_save{i};
    centrecoord=[0 0];
    p = 2 * ones( size(orimesh, 1) - 1, 1);
    [supel] = initialize_problem(orimesh, centrecoord, p);
    [local,supel,globnodesdisp] = today3_modify(p, supel);
    plot_stress(local, local);
end