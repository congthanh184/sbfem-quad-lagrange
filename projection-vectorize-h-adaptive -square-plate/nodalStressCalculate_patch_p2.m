function [unpackStress] = nodalStressCalculate_patch_p2(p, orimesh, local, supel, mainNodes, ngl)
%

    % mainNodes = [1; 2; 3; 4; 5];

    [matmtrx,G] = sbfematisolq(1,supel.emodule,supel.poisson);
    D = G*matmtrx;

    daNodes = orimesh(:, 1);

    [point1, ~] = sbfeglqd1(ngl);
    nodalStress = [];
    for iMainNodes = 1:numel(mainNodes)-1
        [~, idInSupelNodes] = ismember(mainNodes(iMainNodes:iMainNodes+1), daNodes);
        startIdx = idInSupelNodes(1);
        endIdx = idInSupelNodes(2);

        sampleX = []; 
        sampleY = [];
        sampleTest = [];
        for iel=startIdx:endIdx-1
            for iP = 1:numel(point1)
                [x, y, gaussStr] = calculateModalStress(D, iel, supel, local, p, point1(iP));
                sampleTest = [sampleTest; gaussStr(:)'];
                sampleY = [sampleY; y];
                sampleX = [sampleX; x];
            end
        end

        % nodalX = supel.gcoord(startIdx:endIdx, 1);
        % nodalY = supel.gcoord(startIdx:endIdx, 2);
        nodalX = orimesh(startIdx:endIdx, 2);
        nodalY = orimesh(startIdx:endIdx, 3);
        % check which value is the same, then switch with X
        if numel(unique(round(sampleX,5))) == 1
            sampleX = sampleY;
            nodalX = nodalY;
        end


        % stressTempPatch = interpolateStress( sampleX, sampleTest, nodalX, 2);
        
        % something wrong here, should produce sigma_star with the same order p
        % for example, if p = 1, that element should have 2 sigma_star
        % if p = 2, then have 3 sigma_star
        if (endIdx - startIdx == 1)
%         if size(sampleTest, 1) < 3
            nodalX = [nodalX(1); mean(nodalX); nodalX(end)];
            stressTempPatch = interpolateStress( sampleX, sampleTest, nodalX, 1);            
        else
            stressTempPatch = [];
            for iel = 1:numel(nodalX)-2
                startIdxTemp = (iel-1)*ngl + 1;
                endIdxTemp = startIdxTemp + 2*ngl - 1;

                nodalXTemp = nodalX(iel:iel+2);
                for idx = numel(nodalXTemp):-1:2
                    nodalXTemp = [nodalXTemp(1:idx-1); mean([nodalXTemp(idx), nodalXTemp(idx-1)]); nodalXTemp(idx:end)];
                end

                stressTemp = interpolateStress( sampleX(startIdxTemp:endIdxTemp, :), sampleTest(startIdxTemp:endIdxTemp, :), ... 
                                                nodalXTemp, 2);            
                if isempty(stressTempPatch)
                    stressTempPatch = stressTemp;
                else
                    stressTempPatch = [ stressTempPatch(1:end-2, :); ...                                       
                                        mean([stressTempPatch(end-1, :); stressTemp(2, :)]); ...
                                        stressTemp(3:end, :)];
                end
            end
        end

        
        if isempty(nodalStress)
            nodalStress = stressTempPatch;
        else
            nodalStress = [ nodalStress(1:end-1, :); ... 
                            (nodalStress(end, :) + stressTempPatch(1, :))./2; ...
                            stressTempPatch(2:end, :)];
        end

    end

    [x, y, gaussStr] = calculateModalStress(D, 1, supel, local, p, -1);
    gaussTemp1st = gaussStr(:)';

    [x, y, gaussStr] = calculateModalStress(D, numel(p), supel, local, p, 1);
    gaussTemplast = gaussStr(:)';

    nodalStress = [ gaussTemp1st; ...
                    nodalStress(2:end-1, :); ...
                    gaussTemplast];
                    
    unpackStress = [];
    numRow = 3;
    numCol = size(nodalStress, 2) / numRow;
    for iSRow = 1:size(nodalStress, 1)
        unTemp = reshape(nodalStress(iSRow, :)', 3, numCol);
        unpackStress = [unpackStress; unTemp];
    end

% [stress_field] = recoverFromModalStress(unpackStress, supel, local, p);
end