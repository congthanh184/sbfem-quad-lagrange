figure;
hold on;

for iSi = 1:11
    dispx = [];
    dispy = [];
    for jEle = 1:size(local,2)
        x = local(jEle).x(iSi, :) + local(jEle).dispx(iSi, :);
        y = local(jEle).y(iSi, :) + local(jEle).dispy(iSi, :);
%         x =  local(jEle).x(iSi, :);
%         y =  local(jEle).y(iSi, :);

        dispx = [dispx; x(:)];
        dispy = [dispy; y(:)];
        
    end
    plot( dispx, dispy);
end

% figure;
% hold on;
% [X, Y] = meshgrid(0:1:2000);
% Z = size(X);
% for iSi = 1:11
%     dispx = [];
%     dispy = [];
%     stressxy = [];
%     for jEle = 1:4
%         x = local(jEle).x(iSi, :);
%         y = local(jEle).y(iSi, :);
% %         x =  local(jEle).strx(iSi, :);
% %         y =  local(jEle).stry(iSi, :);
%         stress_temp = local(jEle).strx(iSi, :);
%         dispx = [dispx; x(:)];
%         dispy = [dispy; y(:)];
%         stressxy = [stressxy; stress_temp(:)];
%     end
%     for k = 1:numel(dispx)
%         Z( dispx(k), dispy(y) ) = stressxy(k);
%     end
% %     strssxy = diag(stressxy);
% %     contour( dispx, dispy, strssxy);
%     plot(dispx, dispy);
% end
% 
% [point1, wt1] = sbfeglqd1(64);
% numel(point1)
% res = [];
% subplot(121);
% plot(point1);
% subplot(122);
% plot(wt1);
% % xsys = [1000; 2000; 0; 2000];
% xsys = [0; 2000; 0; 0];
% xy = [];

% for i = 1:numel(point1)
%     [N, NN, dN] = lobatto(point1(i), 2);
%     
%     res = [res; NN];  
%     
% %     temp = N * xsys;
% %     xy= [xy; temp];
% end
% % 
% plot(res(:, 2))
% hold on;
% plot(res(:, 1), 'x')
% plot(res(:, 3))
% 
% x = xy(1:2:numel(xy));
% y = xy(2:2:numel(xy));
% 
% plot(x, y, 'x');
% 
% p = 2;
% % nita = sbfeglqd1(64);
% nita = -1:0.1:1;
% 
% x = [0; 0; 0];
% y = [2000; 0; 0];
% 
% cenCoordx = 1000;
% cenCoordy = 1000;
% 
% lob_save = [];
% u_nita = []
% for i = 1:numel(nita)
%     [N, NN] = lobatto(nita(i), p);
%     lob_save = [lob_save; NN];
%     tx = cenCoordx + NN * (x-cenCoordx);
%     ty = cenCoordy + NN * (y - cenCoordy);
%     u_nita = [u_nita; [tx, ty]];
% end
% 
% figure;
% hold on;
% 
% for i = 1:p+1
%     plot(lob_save(:, i));
% end
% 
% figure;
% plot(u_nita(:, 1), u_nita(:, 2), 'x')