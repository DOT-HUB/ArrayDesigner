% This function computes the fluence in the 10-2.5 locations on the
% volumetric head mesh with Toast
% 
% INPUTS: 
%
% headVolumeMesh          The multi-layer volume mesh structure. 
%                         Contains fields: node, face, elem, labels
%
% landmarks               A matrix containing the landmarks coordinate
%
% ADSolutionSpace         ten-two.five locations for the mesh (.positions
%                         (nx3) and .labels {n})
%
% optProp                 Struct with field prop, one entry for each tissue
%                         type, following order of headVolumeMesh.labels.
%                         It should contain, for each tissue, its optical
%                         properties, order as mus, anisotropy, mua, refractive index
%                         (considered as homogeneous, we expect the same for each tissue type)
%
% pathnameSave            pathname to the Phis folder in the AD folder
%                         within the specific head model folder
%
%
% Dependencies: TOAST
%

function computeFluenceToast(pathnameMSHS,ADSolutionSpacePath,pathnameSave)

disp('COMPUTING FLUENCE DISTRIBUTIONS...');

if ~exist(pathnameSave,'dir')
    mkdir(pathnameSave)
end

% Load head volume
load(pathnameMSHS,'-mat');

% get AD solution space
fid = fopen(fullfile(ADSolutionSpacePath,'Utils','scalpSolutionSpace_10_2p5.txt'),'r');
tmp = textscan(fid,'%.6f %.6f %.6f');
fclose(fid);  
ADSolutionSpace.positions = [tmp{1} tmp{2} tmp{3}];

% Create toast mesh
eltp = ones(length(headVolumeMesh.elem),1)*3;
headVolumeMesh.elem(:,1:4) = [headVolumeMesh.elem(:,4), headVolumeMesh.elem(:,1), headVolumeMesh.elem(:,2), headVolumeMesh.elem(:,3)];
hMesh = toastMesh(headVolumeMesh.node(:,1:3),headVolumeMesh.elem(:,1:4),eltp);

nNodes = hMesh.NodeCount;
n_tissue = length(unique(headVolumeMesh.elem(:,5)));

% Optical properties 
c0 = 0.3;   % speed of light in vacuum
nWavs = 1;

% First assign
mua = ones(nNodes,1);
mus = ones(nNodes,1);
ref = ones(nNodes,1);

% Determine tissue optical properties and populate vectors
for tiss = 1:n_tissue
    tmpInd = find(strcmpi({'scalp','skull','ECT','CSF','GM','WM'},headVolumeMesh.labels{tiss}), 1);
    if isempty(tmpInd)
        error('Unknown tissue label, please correct rmap.headVolumeMesh.labels and try again');
    end
    tissueNodeList = headVolumeMesh.node(:,4)==tiss;
    for wav = 1:nWavs
        [muaT, musPrimeT, refIndT] = DOTHUB_getTissueCoeffs(headVolumeMesh.labels{tiss},800);%Potentially update this to be age-specific?
        mua(tissueNodeList,1) = muaT;
        mus(tissueNodeList,1) = musPrimeT;
        ref(tissueNodeList,1) = refIndT;
    end
end

c_medium = c0./ref;

% Create mus and mua vectors
% % for i_t = 1:n_tissue
% %     mua(headVolumeMesh.node(:,4)==i_t,1) = optProp(i_t).prop(3);
% %     mus(headVolumeMesh.node(:,4)==i_t,1) = optProp(i_t).prop(1).*(1-optProp(i_t).prop(2));
% % end

% Calcualte boundary mismatch factors
theta_c = asin(1./ref(1));
R0 = ((ref(1)-1).^2)/((ref(1)+1).^2);
A = (2./(1-R0) - 1 + abs(cos(theta_c)).^3)./(1-abs(cos(theta_c)).^2);
zeta = (1./(2.0.*A)).*ones(nNodes,1);

% Build the system matrix assuming mus independent of mua (baseline and
% perturbation).
K = hMesh.SysmatComponent('PDD', c_medium(1).*1./(3*(mua+mus)));

F = hMesh.SysmatComponent('BNDPFF', c_medium(1).*ones(nNodes,1).*zeta);

M = hMesh.SysmatComponent('PFF', mua.*c_medium(1));

S = K+M+F;

clear K F M

%Loop around ADSolutionSpace positions as sourcepos and save every 100
measpos = ADSolutionSpace.positions(1,:); % Use Nasion as measpos, this shouldn't matter; SB: changed to first pos
count = 1;
batchsize = 100;
for i = 1:batchsize:size(ADSolutionSpace.positions,1)
    
    fprintf('Running batch %d of %d\n',count,floor(size(ADSolutionSpace.positions,1)/batchsize));
    
    if i+batchsize-1>size(ADSolutionSpace.positions,1)
        finind = size(ADSolutionSpace.positions,1);
    else
        finind = i+batchsize-1;
    end
    
    sourcepos = ADSolutionSpace.positions(i:finind,:);
    
    % Get the source vectors, which in toast are also not multiplied by the
    % speed of light.
    hMesh.SetQM(sourcepos, measpos);
    qvec = hMesh.Qvec('Neumann', 'Gaussian', 2);
    
    % Forward and adjoint field for all sources and detectors, baseline mua
    phi = S\qvec;
        
    fprintf('Saving batch %d\n',count);
    fname = ['Phi_Batch' num2str(count) '.mat'];
    save(fullfile(pathnameSave,fname),'phi','-v7.3');
    
    count = count+1;
      
end