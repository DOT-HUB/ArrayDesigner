function hAxis = plotROI(pathnameHeadModel,pathnameROI,varargin)

%This function plots an Array Designer ROI on the GM surface mesh
%
%############################### INPUTS ###################################
%
% pathnameHeadModel     =   the path of the HeadModel (contained in
%                           ArrayDesigner/HeadModels) OR the GMSurfaceMesh
%                           structure
%
% pathnameROI           =   the path to the ROI.txt (which should also be 
%                           saved within the ROIs folder in the selected 
%                           HeadModel directoryfile OR the ROI structure
%                           that is output by loadROI
%
% viewAng               =   Optional. Viewing angle [a,b] Defaults to [

%############################# Dependencies ###############################
% Iso2Mesh
%
% #########################################################################
% RJC, UCL, Dec 2020
% 
% ############################## Updates ##################################
% #########################################################################
%
% ############################### TO DO ###################################
% #########################################################################

% Manage Inputs ############################### 
varInputs = inputParser;
addParameter(varInputs,'hAxes','',@ishandle);
addParameter(varInputs,'view',[],@isnumeric);
parse(varInputs,varargin{:});
varInputs = varInputs.Results;
if isempty(varInputs.hAxes)
    hAxes = gca;
end
if isempty(varInputs.view)
    viewAng = [0 90];
end

if isstr(pathnameROI)
    ROI = loadROI(pathnameROI);
elseif isstruct(pathnameROI)
    ROI = pathnameROI;
end

if ~isempty(pathnameHeadModel)
    if isstr(pathnameHeadModel)
        load([pathnameHeadModel '/GMSurfaceMesh.mat']);
    else
    GMSurfaceMesh = pathnameHeadModel;
    end
else
    error('Please specify HeadModel (2nd input)');
end

load('CMAPgreyGreen.mat')

% Reorg
nodes = GMSurfaceMesh.node;
face  = GMSurfaceMesh.face;
ROImask = zeros(length(nodes),1);
ROImask(ROI.gmNodeList) = 1;

% Plot ############################### 

hAxis = gca;
hPatch = trisurf(face(:,1:3), nodes(:,1), nodes(:,2), nodes(:,3),ROImask,'EdgeColor',[0.8 0.8 0.8],'EdgeAlpha',1);
shading('interp');
set(hPatch,'diffusestrength',.7,'specularstrength',.2,'ambientstrength',.2);
set(hPatch,'Facelighting','phong');
view(viewAng);
camlight(viewAng(1),viewAng(2));
camlight(viewAng(1)+90,0);
camlight(viewAng(1)+180,0);
camlight(viewAng(1)+270,0);
axis equal;axis off;

colormap(greyGreen);
caxis([0 1]);
colorbar off
axis off
view(viewAng);
