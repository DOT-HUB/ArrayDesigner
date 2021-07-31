% This function computes the PMDFs on the GM surface for all possible 10-2.5 combinations
% within the SDrange and save them as PMDFs.data and PMDFs.idx in the PMDFs folder
% within the specified head model folder. They are correctly scaled. 
% It also saves the avNodalVol.txt file in the specified head model folder. 
% 
% INPUTS: 
%
% mshs                      :  A file containing a structure of:
%
%                           % headVolumeMesh     :   The multi-layer volume mesh structure. 
%                                                    Contains fields: node, face, elem, labels
%
%                           % gmSurfaceMesh      :   The gm surface mesh structure. 
%                                                    Contains fields: node, face.
%
%                           % scalpSurfaceMesh   :   The scalp surface mesh structure.
%                                                    Contains fields: node, face.
%
%                           % vol2gm             :   The sparse matrix mapping from head volume mesh
%                                                    space to GM surface mesh space
%
%                           % landmarks          :   A matrix containing the
%                                                    landmarks coordinate
%
%                           % tenFive            :   ten-five locations for
%                                                    the mesh (.positions
%                                                    (nx3) and .labels {n})
%
%                           % ADSolutionSpace    :   ten-two.five locations for
%                                                    the mesh (.positions
%                                                    (nx3) and .labels {n})
%
%                           % logData            :   optional
%
%                           % fileName           :   The path of the saved mshs file
%
% SDrange                 range of distance for meaningful PMDFs (e.g.,[0
%                         60])
%
% pmdfThresh              threhsold to set to 0 smaller values of PMDF
%                         (e.g., 1e-6)
%
% pathnamePhis            pathname to the PMDFs folder
%
% OUTPUT:
%
% mshs                      :  A file containing a structure of:
%
%                           % headVolumeMesh     :   The multi-layer volume mesh structure. 
%                                                    Contains fields: node, face, elem, labels
%
%                           % gmSurfaceMesh      :   The gm surface mesh structure. 
%                                                    Contains fields: node, face.
%
%                           % gmSurfaceMeshDownsampled: The gm surface mesh structure downsampled by a factor of 2.5. 
%                                                    Contains fields: node, face.
%
%                           % scalpSurfaceMesh   :   The scalp surface mesh structure.
%                                                    Contains fields: node, face.
%
%                           % vol2gm             :   The sparse matrix mapping from head volume mesh
%                                                    space to GM surface mesh space
%
%                           % vol2gmDownsampled  :   The sparse matrix mapping from head volume mesh
%                                                    space to GM downsampled surface mesh space
%
%                           % landmarks          :   A matrix containing the
%                                                    landmarks coordinate
%
%                           % tenFive            :   ten-five locations for
%                                                    the mesh (.positions
%                                                    (nx3) and .labels {n})
%
%                           % ADSolutionSpace    :   ten-two.five locations for
%                                                    the mesh (.positions
%                                                    (nx3) and .labels {n})
%
%                           % PMDFNormFactor     :   normalization factor
%                                                    required for the PMDFs weight computation for SNR 
%
%                           % logData            :   optional
%
%                           % fileName           :   The path of the saved mshs file
%
% Dependencies: TOAST, iso2mesh
%

function [mshs] = computePMDFs(pathnamePhis,mshs,SDrange,pmdfThresh)

files = dir(fullfile(pathnamePhis,'*.mat'));

phi = [];
for i = 1:length(files)
    inname = fullfile(pathnamePhis,['Phi_Batch' num2str(i) '.mat']);
    phis_tmp = load(inname);
    phi = [phi phis_tmp.phi];
end
clear phis_tmp

posScalp = mshs.ADSolutionSpace.positions;
nposScalp = size(posScalp,1); % number of scalp lattice points

% Downsample GMSurfaceMesh with a factor of 2.5 and save it also as
% GMpos.txt
if ~isfield(mshs,'gmSurfaceMeshDownsampled')
    [mshs.gmSurfaceMeshDownsampled.node,mshs.gmSurfaceMeshDownsampled.face] = remeshsurf(mshs.gmSurfaceMesh.node,mshs.gmSurfaceMesh.face,2.5);
    
    % Save as GMpos.txt to be fed to the algorithm
    fid = fopen(fullfile(pathnamePhis(1:end-6),'GMpos.txt'),'w');
    for iV = 1:size(mshs.gmSurfaceMeshDownsampled.node,1)
        fprintf(fid,'%.6f %.6f %.6f \n',mshs.gmSurfaceMeshDownsampled.node(iV,1:3));
    end
    fclose(fid); 
end

% Compute new vol2gm matrix and nodal volume vector
if ~isfield(mshs,'vol2gmDownsampled')
    [mshs.vol2gmDownsampled,GMnodalVol] = vol2gmDirectVol(mshs.headVolumeMesh,mshs.gmSurfaceMeshDownsampled.node);
end

% fluence distribution = 1 x nposGM.
phiSurf = mshs.vol2gmDownsampled*phi;
% Output the PMDF normalization factors for use later (in SNR estimation)
PMDFnormFactor = zeros(nposScalp,nposScalp); %Save scaling factor (max val of src distribution at detector location or vice versa)
% Dist mat
distMat = zeros(nposScalp,nposScalp);

% Compute distance matrix between ADSolutionSpace points
noskipmat = ones(nposScalp,nposScalp);
noskipmat(1:nposScalp+1:end) = 0;

h1 = waitbar(0,'Calculating Dist Mat...');
for ii = 1:nposScalp
    waitbar(ii/nposScalp,h1)
    
    for jj = 1:nposScalp
    
        if noskipmat(ii,jj)
            tmpDist = sqrt( sum(  (posScalp(ii,:) - posScalp(jj,:)).^2,2));
            distMat(ii,jj) = tmpDist;
            distMat(jj,ii) = tmpDist;
            noskipmat(ii,jj) = 0;
            noskipmat(jj,ii) = 0;
        end
        
    end
end
delete(h1);

% Compute normalization factor
noskipmat = ones(nposScalp,nposScalp);
noskipmat(1:nposScalp+1:end) = 0;

h1 = waitbar(0,'Calculating PMDF norms...');
for ii = 1:nposScalp
    waitbar(ii/nposScalp,h1)
    
    for jj = 1:nposScalp
    
        if distMat(ii,jj) < SDrange(1) || distMat(ii,jj) > SDrange(2) % keep only if within SD range
            noskipmat(ii,jj) = 0;
            noskipmat(jj,ii) = 0;
        end
        
        if noskipmat(ii,jj)
            
            [~,mxIndq] = max(abs(phi(:,ii)));
            [~,mxIndm] = max(abs(phi(:,jj)));
            PMDFnormFactor(ii,jj) = min([phi(mxIndm,ii); phi(mxIndq,jj)]);
            PMDFnormFactor(jj,ii) = PMDFnormFactor(ii,jj);
            
            noskipmat(ii,jj) = 0;
            noskipmat(jj,ii) = 0;       
        
        end
    end
end
delete(h1)
clear phi; %save some memory

% Compute PMDFs in the surface, correctly scaled.  Cell of sparse arrays
ALLPMDFs = cell(nposScalp,nposScalp);

noskipmat = ones(nposScalp,nposScalp);
noskipmat(1:nposScalp+1:end) = 0;

tic;
h1 = waitbar(0,'Calculating PMDFs...');
for ii = 1:nposScalp
    waitbar(ii/nposScalp,h1)
    
    for jj = 1:nposScalp
             
        if distMat(ii,jj)  < SDrange(1) || distMat(ii,jj)  > SDrange(2) % keep only if within SD range
            noskipmat(ii,jj) = 0;
            noskipmat(jj,ii) = 0;
        end
        
        if noskipmat(ii,jj)

            pmdf = (phiSurf(:,ii).*phiSurf(:,jj).*GMnodalVol)./PMDFnormFactor(ii,jj);
            % Correctly scaled PMDF, but positive (not point being negative)
                        
            % Get rid of values below threshold
            pmdfmask = zeros(size(pmdf));
            maskInd = pmdf > max(pmdf)*pmdfThresh;
            pmdfmask(maskInd) = pmdf(maskInd);
            
            ALLPMDFs{ii,jj} = sparse(pmdfmask);
            
            noskipmat(ii,jj) = 0;
            noskipmat(jj,ii) = 0;
        end
    end
    
    t = toc;
    remaining = sum(noskipmat(:));
    done = nposScalp.^2 - remaining;
    
    fprintf(['Estimated time remaining = ' num2str((t/done)*remaining/60,'%0.1f') ' mins\n'])
end
delete(h1);

mshs.PMDFnormFactor = PMDFnormFactor;

% Save the PMDFs file in the .data and .idx format
fid = fopen(fullfile(pathnamePhis,'PMDFs.idx'),'w');
for iR = 1:size(ALLPMDFs,2)
    %fprintf('%d \n',iR)
    for iC = 1:size(ALLPMDFs,1)
        if iC > iR
            if ~isempty(ALLPMDFs{iR,iC})
                nTot = sum(full(ALLPMDFs{iR,iC} >= 1e-9));
                fprintf(fid,'%d ',nTot);
            else
                fprintf(fid,'%d ',0);
            end
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);

h1 = waitbar(0,'Saving PMDFs...');
fid = fopen(fullfile(pathnamePhis,'PMDFs.data'),'w');
for iR = 1:size(ALLPMDFs,2)
    waitbar(iR/size(ALLPMDFs,2),h1)
    for iC = 1:size(ALLPMDFs,1)
        if ~isempty(ALLPMDFs{iR,iC})
            for iG = 1:length(ALLPMDFs{iR,iC})
                if ALLPMDFs{iR,iC}(iG) >= 1e-9
                    fwrite(fid,iG-1,'int32');
                    fwrite(fid,full(ALLPMDFs{iR,iC}(iG)),'double');
                end
            end
        end
    end
end
fclose(fid);  
delete(h1);

% Delete the Phi_Batch files (the fluence files) to save space
delete(fullfile(pathnamePhis,'Phi_*'))

% Save the avNodalVol.txt file
avNodalVolume = median(GMnodalVol);
fid = fopen(fullfile(pathnamePhis(1:end-6),'avNodalVol.txt'),'w');
fprintf(fid,'%.15d',avNodalVolume);
close(fid);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol2gm,GMnodalVol] = vol2gmDirectVol(headVolumeMesh,GMNodes)

%This function performs the mapping from a tetrahedral volume mesh to the
%associated GM mesh assigning to each GM node the value of the closest 
% volume node.  Output is in the form of a sparse transformation matrix with
%dimensions NxM where M is the number of tetrahedral mesh nodes and N is
%the number of GM surface nodes;

headMeshNodalVol = nodevolume(headVolumeMesh.node(:,1:3),headVolumeMesh.elem(:,1:4));

GMnodalVol = zeros(size(GMNodes,1),1);
count = 1;
h = waitbar(0,'Calculating vol2gm transform...');
for n = 1:size(GMNodes,1)
    waitbar(n/size(GMNodes,1));
    p = GMNodes(n,:);
    
    dist = sqrt(sum((p-headVolumeMesh.node(:,1:3)).^2,2));
    [~,idx] = min(dist);
     
    i(count) = n;
    j(count) = idx;
    s(count) = 1;
    
    GMnodalVol(n,1) = headMeshNodalVol(idx);
    
    count = count+1;
end
delete(h)

vol2gm = sparse(i,j,s,size(GMNodes,1),size(headVolumeMesh.node,1));

end



