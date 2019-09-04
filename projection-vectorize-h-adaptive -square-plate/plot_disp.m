function plot_disp(local)
figure;
hold on;

for iSi = 1:11
    dispx = [];
    dispy = [];
    for jEle = 1:size(local,2)
        x = local(jEle).x(iSi, :) + 20000*local(jEle).dispx(iSi, :);
        y = local(jEle).y(iSi, :) + 20000*local(jEle).dispy(iSi, :);
%         x =  local(jEle).x(iSi, :);
%         y =  local(jEle).y(iSi, :);

        dispx = [dispx; x(:)];
        dispy = [dispy; y(:)];
        
    end
    plot( dispx, dispy);
end
end