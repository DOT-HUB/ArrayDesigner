function [hAxis, viableChannelDists] = plotArray(AD,A,varargin)

% Manage Inputs ############################### 
varInputs = inputParser;
addParameter(varInputs,'hAxes','',@ishandle);
addParameter(varInputs,'view',[],@isnumeric);
parse(varInputs,varargin{:});
varInputs = varInputs.Results;
if isempty(varInputs.hAxes)
    hAxes = gca;
else
    hAxes = varInputs.hAxes;
end
if isempty(varInputs.view)
    viewAng = [0 90];
else
    viewAng = varInputs.view;
end

%Load HeadModel elements
load([AD.inputs.pathnameHeadModel '/GMSurfaceMesh.mat']);
load([AD.inputs.pathnameHeadModel '/ScalpSurfaceMesh.mat']);
scalpPos = importdata([AD.inputs.pathnameHeadModel '/scalpPos.txt']);

%Unpack
array = [AD.results.source AD.results.detector];
nS = AD.inputs.nD;
nD = AD.inputs.nD;
nodes = GMSurfaceMesh.node;
face  = GMSurfaceMesh.face;
ROI = loadROI(AD.inputs.pathnameROI);

if ~exist('A','var')
    AROIoverlap = zeros(length(nodes),1);
elseif isempty(A)
    AROIoverlap = zeros(length(nodes),1);
else %create ROI overlap map if A exists
    colMap(1,:) = [0.7    0.7    0.7];%grey
    colMap(2,:) = [1 0 1];%magenta - 1 = uncovered ROI
    colMap(3,:) = [0 1 1];%cyan - 2 = superfluous coverage
    colMap(4,:) = [0    1    0];%green - 3 = covered ROI
    Athresh = A>AD.inputs.coverageThresh;
    ROImat = false(size(Athresh));
    ROImat(ROI.gmNodeList) = true;
    AROIoverlap = zeros(size(Athresh));
    AROIoverlap(ROImat) = 1;
    AROIoverlap(Athresh) = 2;
    AROIoverlap(ROImat & Athresh) = 3;
end



% Plot ############################### 
hAxis = gca;
hPatch = trisurf(face(:,1:3), nodes(:,1), nodes(:,2), nodes(:,3),AROIoverlap,'EdgeColor',[0.8 0.8 0.8],'EdgeAlpha',1);
shading('flat');
set(hPatch,'diffusestrength',.7,'specularstrength',.2,'ambientstrength',.2);
set(hPatch,'Facelighting','phong');
view(viewAng);
camlight(viewAng(1),viewAng(2));
camlight(viewAng(1)+90,0);
camlight(viewAng(1)+180,0);
camlight(viewAng(1)+270,0);
axis equal;axis off;
colormap(hAxis,colMap);
caxis([-0.5 3.5])
colorbar off;
hold on;


nchan = 0;
%Add channels to diagram
for i = 1:nS
    dist = sqrt(sum((scalpPos(array(nS+1:end),:) - repmat(scalpPos(array(i),:),nD,1)).^2,2));
    for j = 1:length(dist)
        if dist(j) >= AD.inputs.minRho && dist(j) < AD.inputs.maxRho
            Hl = line([scalpPos(array(i),1) scalpPos(array(nS+j),1)],[scalpPos(array(i),2) scalpPos(array(nS+j),2)],[scalpPos(array(i),3) scalpPos(array(nS+j),3)],'color','m','LineWidth',2);
            nchan = nchan+1;
            viableChannelDists(nchan) = dist(j);
        end
        hold on;
    end
end

plotmesh(ScalpSurfaceMesh.node,ScalpSurfaceMesh.face,'FaceColor','none','FaceAlpha',0,'EdgeColor','k','EdgeAlpha',0.1);

optodeMarkerSize = 100;
Hs = scatter3(scalpPos(array(1:nS),1),scalpPos(array(1:nS),2),scalpPos(array(1:nS),3),optodeMarkerSize,'MarkerEdgeColor','k','MarkerFaceColor','r');
Hd = scatter3(scalpPos(array(nS+1:end),1),scalpPos(array(nS+1:end),2),scalpPos(array(nS+1:end),3),optodeMarkerSize,'MarkerEdgeColor','k','MarkerFaceColor','b');
view([0 90]);
hold off;


% if ~exist('scalpFlag','var');
%     scalpFlag = 1;
% end
% if scalpFlag == 1;
%     plotmesh(ScalpSurfaceMesh.node,ScalpSurfaceMesh.face,'FaceColor','none','FaceAlpha',0,'EdgeColor','k','EdgeAlpha',0.3);
% end

%Add ROI centre
%scatter(roiCent(1),roiCent(2),100,'go','filled');
%axis equal tight
%Hleg = legend([Hs Hd Hl],{'Sources','Detectors','Channels'},'FontSize',14);
%set(Hleg,'Position',[0.8605    0.1564    0.0797    0.0763]);
