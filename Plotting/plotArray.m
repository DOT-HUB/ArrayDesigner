function [hAxes, viableChannelDists] = plotArray(AD,varargin)

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

pathnameMSHS = getMSHSpath(AD.inputs.pathnameHeadModel);
%Load HeadModel elements
load(pathnameMSHS,'gmSurfaceMesh','-mat');
load(pathnameMSHS,'scalpSurfaceMesh','-mat');
scalpPos = importdata([AD.inputs.pathnameHeadModel '/Utils/scalpSolutionSpace_10_2p5.txt']);

%Unpack
array = [AD.results.source AD.results.detector];
nS = AD.inputs.nD;
nD = AD.inputs.nD;
nodes = gmSurfaceMesh.node;
face  = gmSurfaceMesh.face;
ROI = loadROI(AD.inputs.pathnameROI);
sources = AD.results.source;
detectors = AD.results.detector;
nChannels = AD.results.nChannels;
measList = AD.results.measList;

A = AD.results.sensitivityMap;
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
hPatch = trisurf(face(:,1:3), nodes(:,1), nodes(:,2), nodes(:,3),AROIoverlap,'EdgeColor',[0.8 0.8 0.8],'EdgeAlpha',1,'Parent',hAxes);
shading(hAxes,'flat');
set(hPatch,'diffusestrength',.7,'specularstrength',.2,'ambientstrength',.2);
set(hPatch,'Facelighting','phong');

camlight(hAxes,viewAng(1),viewAng(2));
camlight(hAxes,viewAng(1)+90,0);
camlight(hAxes,viewAng(1)+180,0);
camlight(hAxes,viewAng(1)+270,0);
hAxes.DataAspectRatio = [1 1 1];
colormap(hAxes,colMap);
caxis(hAxes,[-0.5 3.5]);
view(hAxes,viewAng);
hAxes.Visible = 'off';
hold(hAxes,'on');

%Add channels to diagram
for i = 1:nChannels
    Hl = line(hAxes,[scalpPos(sources(measList(i,1)),1) scalpPos(detectors(measList(i,2)),1)],[scalpPos(sources(measList(i,1)),2) scalpPos(detectors(measList(i,2)),2)],[scalpPos(sources(measList(i,1)),3) scalpPos(detectors(measList(i,2)),3)],'color','m','LineWidth',2);
end

%Add scalp
hPatch = trisurf(scalpSurfaceMesh.face(:,1:3), scalpSurfaceMesh.node(:,1), scalpSurfaceMesh.node(:,2), scalpSurfaceMesh.node(:,3),'FaceColor','none','FaceAlpha',0,'EdgeColor','k','EdgeAlpha',0.05,'Parent',hAxes);

optodeMarkerSize = 100;
Hs = scatter3(hAxes,scalpPos(sources,1),scalpPos(sources,2),scalpPos(sources,3),optodeMarkerSize,'MarkerEdgeColor','k','MarkerFaceColor','r');
Hd = scatter3(hAxes,scalpPos(detectors,1),scalpPos(detectors,2),scalpPos(detectors,3),optodeMarkerSize,'MarkerEdgeColor','k','MarkerFaceColor','b');
view(hAxes,viewAng);
hold(hAxes,'off');

hAxes.XLimMode = 'auto';
hAxes.YLimMode = 'auto';
hAxes.ZLimMode = 'auto';
set(hAxes, 'XLimSpec', 'Tight','YLimSpec','Tight');

%Add colorbar
cb = colorbar(hAxes);
width = 0.015;
height = 0.1;
cb.Position = [0.860, 0.01, width, height];
cb.Ticks = [1 2 3];
cb.TickLabels = {'Uncovered','Additional','Covered'};
cb.FontSize = 10;
