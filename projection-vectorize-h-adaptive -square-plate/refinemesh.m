function [newMesh] = refinemesh(oriMesh, number_refine)
    newMesh = oriMesh;

    for id = 1:number_refine
        newNodeIdx = max(newMesh(:, 1));

        for idxNewMesh = size(newMesh, 1):-1:2
            newNodeIdx = newNodeIdx + 1;
            newNode = [newNodeIdx, mean(newMesh(idxNewMesh-1:idxNewMesh, 2:3))];
            newMesh = [newMesh(1:idxNewMesh-1, :); newNode; newMesh(idxNewMesh:end, :)];
        end
    end
end