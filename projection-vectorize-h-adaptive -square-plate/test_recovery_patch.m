
mainNodes = [1; 2; 3; 4; 5];

[matmtrx,G] = sbfematisolq(1,supel.emodule,supel.poisson);
D = G*matmtrx;

nodalStress = [];
for iMainNodes = 1:numel(mainNodes)-1
    [~, idInSupelNodes] = ismember(mainNodes(iMainNodes:iMainNodes+1), supel.nodes);
    startIdx = idInSupelNodes(1);
    endIdx = idInSupelNodes(2);

    sampleX = []; 
    sampleY = [];
    sampleTest = [];
    for iel=startIdx:endIdx-1
        if (startIdx - endIdx) > 1
            [x, y, gaussStr] = calculateModalStress(D, iel, supel, local, p, 0);
        else
            [x, y, gaussStr] = calculateModalStress(D, iel, supel, local, p, -0.5774);
            sampleTest = [sampleTest; gaussStr(:)'];
            sampleY = [sampleY; y];
            sampleX = [sampleX; x];
            [x, y, gaussStr] = calculateModalStress(D, iel, supel, local, p, 0.5774);            
        end
        sampleTest = [sampleTest; gaussStr(:)'];
        sampleY = [sampleY; y];
        sampleX = [sampleX; x];
    end
    nodalX = supel.gcoord(startIdx:endIdx, 1);
    nodalY = supel.gcoord(startIdx:endIdx, 2);
    % check which value is the same, then switch with X
    if numel(unique(sampleX))==1
        sampleX = sampleY;
        nodalX = nodalY;
    end

    orderInterp = 1;
    if numel(sampleX) > 2
        orderInterp = 2;
    end
    stressTempPatch = interpolateStress( sampleX, sampleTest, nodalX, orderInterp);
    
    if isempty(nodalStress)
        nodalStress = stressTempPatch;
    else
        nodalStress = [ nodalStress(1:end-1, :); ... 
                        (nodalStress(end, :) + stressTempPatch(1, :))./2; ...
                        stressTempPatch(2:end, :)];
    end

end

unpackStress = [];
numRow = 3;
numCol = size(nodalStress, 2) / numRow;
for iSRow = 1:size(nodalStress, 1)
    unTemp = reshape(nodalStress(iSRow, :), 3, numCol);
    unpackStress = [unpackStress; unTemp];
end

% [stress_field] = recoverFromModalStress(unpackStress, supel, local, p);


% [locMainNodes, ~] = ismember(supel.nodes, 1:5);
% nodalY = supel.gcoord(locMainNodes, 2);
% stressPatch1 = interpolateStress( sampleY, sampleMidStrY, nodalY, 2);
% stressPatchTest = interpolateStress( sampleY, sampleTest, nodalY, 2);

% sampleMidStrX = [];
% sampleMidStrY = [];
% sampleMidStrXY = [];
% for iel=5:8
%     [x, y, gaussStr] = calculateModalStress(D, iel, supel, local, p, 0);
%     sampleMidStrX = [sampleMidStrX; gaussStr(1,:)];
%     sampleMidStrY = [sampleMidStrY; gaussStr(2,:)];
%     sampleMidStrXY = [sampleMidStrXY; gaussStr(3,:)];
%     sampleX = [sampleX; x];
% end
% [locMainNodes, ~] = ismember(supel.nodes, 5:9);
% nodalX = supel.gcoord(locMainNodes, 1);
% stressPatch2 = interpolateStress( sampleX, sampleMidStrY, nodalX, 2);

% [x, y, strNode1] = calculateModalStress(D, 1, supel, local, p, -1);
% [x, y, strNode9] = calculateModalStress(D, 8, supel, local, p, 1);
% %%
% % ssY = [	(strNode1(2, :) + stressPatch1(1, :))./2; ...
% % 		stressPatch1(2:end-1, :); ...
% % 		(stressPatch1(end, :) + stressPatch2(1, :))./2; ...
% % 		stressPatch2(2:end-1, :); ...
% % 		(strNode9(2, :) +stressPatch2(end, :))./2 ];
% ssY = [	
% 		stressPatch1(1:end-1, :); ...
% 		(stressPatch1(end, :) + stressPatch2(1, :))./2; ...
% 		stressPatch2(2:end, :); ...
% 		];
% % ssY = [	stressPatch1(1:end, :); ...
% % 		
% % 		stressPatch2(2:end, :) ];
    
% sdoff = supel.nnode*2 + 2*(sum(p)-supel.nel);
% sigma_star = zeros(supel.nnode*3, sdoff);
% sigma_star(2:3:supel.nnode*3, :) = ssY;

% [stress_field] = recoverFromModalStress(sigma_star, supel, local, p);
