
clear all
mshs = load('MNI152_Symmetric.mshs','-mat');
mshs = rmfield(mshs,'fileName');
mshs = rmfield(mshs,'logData');

gmSurfaceMesh = mshs.gmSurfaceMesh;

%Index of zero is a node not assigned to a parcel.
%Reset label list
parcelIndTmp = unique(nonzeros(gmSurfaceMesh.parcellation.index));
gmSurfaceMesh.parcellation.labels = gmSurfaceMesh.parcellation.labels(parcelIndTmp);

%Fix holes in parcels?
nNode = length(gmSurfaceMesh.node);
nParcels = length(parcelIndTmp);
newIndex = gmSurfaceMesh.parcellation.index;
h = waitbar(0,'Running through parcels');
for i = 1:length(parcelIndTmp)
    ind = parcelIndTmp(i);
    nodeMaskTmp = gmSurfaceMesh.parcellation.index == ind;
    nodeListTmp = find(nodeMaskTmp);
    
    %Add existing parcel nodes to new index
    newIndex(nodeListTmp) = i;
    %
    %figure;
    %subplot(1,2,1);
    %DOTHUB_plotSurfaceImage(gmSurfaceMesh,double(nodeMaskTmp),[],'flat');
    
    %if the nodes associated with a face are parts of faces that themselves have >=2 nodes in a parcel, those nodes should also be in the parcel
    faceList = ismember(gmSurfaceMesh.face(:,1:3),nodeListTmp); %index of if a face has nodes in the parcel
    faceListParcel = any(double(faceList'))';
    nodeList = gmSurfaceMesh.face(faceListParcel,1:3);
    nodeList = unique(nodeList(:));
    tmp = false(nNode,1);
    tmp(nodeList) = 1;
    newIndex(tmp & gmSurfaceMesh.parcellation.index==0) = i; %This should be i now
    %nodeMaskTmp = index == i;
    %subplot(1,2,2);
    %DOTHUB_plotSurfaceImage(gmSurfaceMesh,double(nodeMaskTmp),[],'flat');
    waitbar(i/nParcels,h);
end
delete(h)

tmpMask = zeros(nNode,1);
tmpMask(newIndex==60) = 1;
DOTHUB_plotSurfaceImage(gmSurfaceMesh,tmpMask,[],'flat');

figure
subplot(1,2,1);
DOTHUB_plotSurfaceImage(gmSurfaceMesh,gmSurfaceMesh.parcellation.index,[],'flat');
colormap('colorcube');
subplot(1,2,2);
DOTHUB_plotSurfaceImage(gmSurfaceMesh,newIndex,[],'flat');
colormap('colorcube');

%Write out changes
gmSurfaceMesh.parcellation.index = newIndex;
mshs.gmSurfaceMesh = gmSurfaceMesh;

save('MNI152_Symmetric.mshs','-struct','mshs');