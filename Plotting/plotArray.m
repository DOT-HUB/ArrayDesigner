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

if ~exist('A','var')
    A = zeros(length(nodes),1);
elseif isempty(A)
    A = zeros(length(nodes),1);
end

load('CMAPgreyGreen.mat');

% Plot ############################### 
hAxis = gca;
hPatch = trisurf(face(:,1:3), nodes(:,1), nodes(:,2), nodes(:,3),A,'EdgeColor',[0.8 0.8 0.8],'EdgeAlpha',1);
shading('interp');
set(hPatch,'diffusestrength',.7,'specularstrength',.2,'ambientstrength',.2);
set(hPatch,'Facelighting','phong');
view(viewAng);
camlight(viewAng(1),viewAng(2));
camlight(viewAng(1)+90,0);
camlight(viewAng(1)+180,0);
camlight(viewAng(1)+270,0);
axis equal;axis off;
caxis([0 1]);
colormap(greyGreen);
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

Hs = scatter3(scalpPos(array(1:nS),1),scalpPos(array(1:nS),2),scalpPos(array(1:nS),3),70,'MarkerEdgeColor','k','MarkerFaceColor','r');
Hd = scatter3(scalpPos(array(nS+1:end),1),scalpPos(array(nS+1:end),2),scalpPos(array(nS+1:end),3),70,'MarkerEdgeColor','k','MarkerFaceColor','b');
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
