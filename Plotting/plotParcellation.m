function hAxis = plotParcellation(pathnameHeadModel,parcellationNumber,varargin)

%This function plots a parcellation on the GM surface mesh
%
%############################### INPUTS ###################################
%
% pathnameHeadModel     =   the path of the HeadModel (contained in
%                           ArrayDesigner/HeadModels) OR the GMSurfaceMesh
%                           structure
%
% parcellationIndex     =   index of parcellation withing headmodel to use.
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
else
    viewAng = varInputs.view;
end

if ~isempty(pathnameHeadModel)
    if isstr(pathnameHeadModel)
        pathnameMSHS = getMSHSpath(pathnameHeadModel);
        load(pathnameMSHS);
        GMSurfaceMesh = gmSurfaceMesh;
    else
        GMSurfaceMesh = pathnameHeadModel;
    end
else
    error('Please specify HeadModel');
end

if ~isfield(GMSurfaceMesh,'parcellation')
    return
elseif isempty(GMSurfaceMesh.parcellation)
        return
end

% Reorg
nodes = GMSurfaceMesh.node;
face  = GMSurfaceMesh.face;
parcellation = GMSurfaceMesh.parcellation(parcellationNumber).index;
nParcels = length(unique(nonzeros(parcellation)));
cmap = 0.7*ones(nParcels+1,3);
cmap(2:end,:) = jet(nParcels);
    
% Plot ############################### 

hPatch = trisurf(face(:,1:3), nodes(:,1), nodes(:,2), nodes(:,3),parcellation,'EdgeColor',[0.8 0.8 0.8],'EdgeAlpha',1,'Parent',hAxes);
shading(hAxes,'flat');
set(hPatch,'diffusestrength',.7,'specularstrength',.2,'ambientstrength',.2);
set(hPatch,'Facelighting','phong');
camlight(hAxes,viewAng(1),viewAng(2));
camlight(hAxes,viewAng(1)+90,0);
camlight(hAxes,viewAng(1)+180,0);
camlight(hAxes,viewAng(1)+270,0);
hAxes.DataAspectRatio = [1 1 1];
colormap(hAxes,cmap);
caxis(hAxes,[0 nParcels]);
view(hAxes,viewAng);
hAxes.Visible = 'off';

