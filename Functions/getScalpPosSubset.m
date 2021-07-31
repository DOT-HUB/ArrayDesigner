function [scalpSubset, subsetInd] = getScalpPosSubset(gmSurfaceMesh,scalpSurfaceMesh,ROIgmNodeList,maxRho)

posGM = gmSurfaceMesh.node;
posROI = posGM(ROIgmNodeList,:);
nposROI = length(posROI);
posScalp = scalpSurfaceMesh.node;
nScalpPos = length(posScalp);

count = 1;
for i = 1:nScalpPos
    
    %First find gmdepth
    [~, ~, gmDepth] = DOTHUB_nearestNode(posScalp(i,:),posGM);
    hyp = sqrt(gmDepth^2 + (maxRho/2)^2);
    
    dists = sqrt(sum((posROI - repmat(posScalp(i,:),nposROI,1)).^2,2));
    if min(dists) < hyp
        scalpSubset(count,:) = posScalp(i,:);
        subsetInd(count,1) = i;
        count = count+1;
    end
end
        