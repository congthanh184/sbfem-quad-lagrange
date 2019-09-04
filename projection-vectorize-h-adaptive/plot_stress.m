function plot_stress(local, rec_str)
    figure; 
    x = [];
    y = [];
    strx = [];
    recx = [];
    for iEle = 1:size(local, 2)
        x = [x; local(iEle).x];    
        y = [y; local(iEle).y];
        strx = [strx; local(iEle).strx];
        recx = [recx; rec_str(iEle).strx];
    end
    contourf(x, y, recx, 8);
    colorbar;
    figure;
    contourf(x, y, strx, 8);
    colorbar;
end