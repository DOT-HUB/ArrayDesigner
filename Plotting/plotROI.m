function hAxis = plotROI(gmSurfaceMesh,pathnameROI,varargin)

%This function plots an Array Designer ROI on the GM surface mesh
%
%############################### INPUTS ###################################
%
% gmSurfaceMesh         =   the gmSurfaceMesh
%                           structure directly
%
% pathnameROI           =   the path to the ROI.txt (which should also be 
%                           saved within the ROIs folder in the selected 
%                           HeadModel directoryfile OR the ROI structure
%                           that is output by loadROI)
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
else
    hAxes = varInputs.hAxes;
end
if isempty(varInputs.view)
    viewAng = [0 90];
end

if isstr(pathnameROI)
    ROI = loadROI(pathnameROI);
elseif isstruct(pathnameROI)
    ROI = pathnameROI;
end

if ~isempty(gmSurfaceMesh)
    if isstr(gmSurfaceMesh)
        load(gmSurfaceMesh,'gmSurfaceMesh','-mat');
    else
        gmSurfaceMesh = gmSurfaceMesh;
    end
else
    error('Please specify HeadModel (2nd input)');
end

colMap(1,:) = [0.7    0.7    0.7];%grey
colMap(2,:) = [1    0    1];%magenta

% Reorg
nodes = gmSurfaceMesh.node;
face  = gmSurfaceMesh.face;
ROImask = zeros(length(nodes),1);
ROImask(ROI.gmNodeList) = 1;

% Plot ############################### 

hPatch = trisurf(face(:,1:3), nodes(:,1), nodes(:,2), nodes(:,3),ROImask,'EdgeColor',[0.8 0.8 0.8],'EdgeAlpha',1,'Parent',hAxes);
shading(hAxes,'flat');
set(hPatch,'diffusestrength',.7,'specularstrength',.2,'ambientstrength',.2);
set(hPatch,'Facelighting','phong');
camlight(hAxes,viewAng(1),viewAng(2));
camlight(hAxes,viewAng(1)+90,0);
camlight(hAxes,viewAng(1)+180,0);
camlight(hAxes,viewAng(1)+270,0);
hAxes.DataAspectRatio = [1 1 1];
colormap(hAxes,colMap);
caxis(hAxes,[-0.5 1.5]);
view(hAxes,viewAng);
hAxes.Visible = 'off';

