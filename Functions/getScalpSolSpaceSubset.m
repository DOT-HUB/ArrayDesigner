function [scalpSubset, subsetInd] = getScalpSolSpaceSubset(gmSurfaceMesh,scalpSolSpace,ROIgmNodeList,maxRho)

posGM = gmSurfaceMesh.node;
posROI = posGM(ROIgmNodeList,1:3);
nposROI = length(posROI);
nScalpPos = length(scalpSolSpace);

count = 1;
for i = 1:nScalpPos
    
    %First find gmdepth
    [~, ~, gmDepth] = DOTHUB_nearestNode(scalpSolSpace(i,:),posGM(:,1:3));
    %Include a scalp pos if it is < hyp away from a ROI GM node where hyp
    %is the hypoteneuse of a triangle formed by the gm-scalp depth and half
    %the maxRho. This approach should include scalp pos in the solution
    %space out far enough that the gm nodes at the edge of an ROI can still be
    %mid way between a maxRho channel. This might cause problems with computation time
    %as the solution space becomes unecessarily large.
    hyp = sqrt(gmDepth^2 + (maxRho/2)^2);
    
    dists = sqrt(sum((posROI - repmat(scalpSolSpace(i,:),nposROI,1)).^2,2));
    if min(dists) < hyp
        scalpSubset(count,:) = scalpSolSpace(i,:);
        subsetInd(count,1) = i;
        count = count+1;
    end
end
        