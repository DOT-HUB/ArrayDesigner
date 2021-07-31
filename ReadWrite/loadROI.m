function ROI = loadROI(pathnameROI)

tmp = importdata(pathnameROI);
ROI.gmNodeList = tmp(1,:)';
ROI.solutionSpaceNodeList = tmp(2,~isnan(tmp(2,:)))';




