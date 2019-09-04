function plot_stress(local, rec_str)
    figure; 
    x = [];
    y = [];
    strx = [];
    recx = [];
    for iEle = 1:size(local, 2)
        x = [x; local(iEle).x];    
        y = [y; local(iEle).y];
        strx = [strx; local(iEle).stry];
        recx = [recx; rec_str(iEle).stry];
    end    
%     figure
%     subplot(1, 2, 1);
    contourf(x, y, strx, 8);
    colorbar;
    figure
%     subplot(1, 2, 2);
    contourf(x, y, recx, 8);
    colorbar;      
end