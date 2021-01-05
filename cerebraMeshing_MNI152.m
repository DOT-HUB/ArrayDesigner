
a = loadnifti('mni_icbm152_CerebrA_tal_nlin_sym_09c.nii');
img = a.NIFTIData;
dim = size(img);
mask = img;
%mask = double(img>0);
mask_filled = zeros(dim);
for i_slice = 1:dim(3)
    mask_filled(:,:,i_slice) = imfill(mask(:,:,i_slice),'holes');
end

%Need to figure out how to get regions into the surface mesh as the 4th
%node dimension. Could make volume mesh first?
for i = 1:101
    opt(i).maxnode = 5000;%500000
    opt(i).radbound = 1;
    opt(i).maxVol = 1;
    opt(i).downSampling = 0.3;
    opt(i).regions = i; 
end
    [node_surf,face_surf] = createSurfMesh(mask_filled,opt);

function [node_s,face_r] = createSurfMesh(mask_tissue,opt)

dim = size(mask_tissue);

% Fill the obtained tissue mask; in this way the surface mesh
% created will have only the external surface of the tissue and not
% the internal boundary one
mask_filled = zeros(dim);
for i_slice = 1:dim(3)
    mask_filled(:,:,i_slice) = imfill(mask_tissue(:,:,i_slice),'holes');
end

% Create surface mesh
[node_surf,face_surf] = vol2surf(mask_filled,1:dim(1),1:dim(2),1:dim(3),opt,1);

% Repairing the mesh, downsampling it and smoothing it
[node_r,face_r] = meshcheckrepair(node_surf,face_surf);
[newnode,newface] = meshresample(node_r,face_r,opt.downSampling); % downsample
[node_r,face_r] = meshcheckrepair(newnode,newface);
conn = meshconn(face_r,size(node_r,1));
node_s = smoothsurf(node_r,[],conn,10,0.7,'lowpass'); % smoothing
[node_s,face_r]=removeisolatednode(node_s,face_r);
[node_s,face_r] = meshcheckrepair(node_s,face_r);
[node_s,face_r] = surfreorient(node_s,face_r);

end