pathnameLoad = '/Users/sabri/Documents/code/ArrayDesigner/ArrayDesigner/HeadModels/MNI152_LowRes/';

% Load GM surface mesh to get n of voxels
load(fullfile(pathnameLoad,'GMSurfaceMesh.mat'))
nVoxels = size(GMSurfaceMesh.node,1);

% Read idx file - each row contains the number of voxels with significant
% values for each PMDF of the upper triangular matrix
fid = fopen(fullfile(pathnameLoad,'PMDFs','PMDFs.idx'),'r');
nPMDFs = fscanf(fid,'%d');
fclose(fid);

% Recover size of square matrix with source-det pairs (upper triangular
% matrix)
sizeMatrix = ceil(sqrt(length(nPMDFs)*2));

% for all p > q -> only upper triangular matrix
allPMDFs = cell(sizeMatrix,sizeMatrix);
fid = fopen(fullfile(pathnameLoad,'PMDFs','PMDFs.data'),'r');
k = 1;
for iR = 1:sizeMatrix
    for iC = 1:sizeMatrix
        if iC > iR 
            if nPMDFs(k) > 0
                allPMDFs{iR,iC} = zeros(nVoxels,1);
                for iG = 1:nPMDFs(k)
                    pos = fread(fid,1,'int32');
                    allPMDFs{iR,iC}(pos+1) = fread(fid,1,'double');
                end
            end
            k = k+1;
        end
    end
end
fclose(fid);

%% Example on how to read specific parts of the .data file
fid = fopen(fullfile(pathnameLoad,'PMDFs','PMDFs.data'),'r');
nBytes = nPMDFs(1)*4 + nPMDFs(1)*8; % int32 -> 4 bytes, double -> 8 bytes
fseek(fid, nBytes, 'bof'); % skip these bytes from origin
pos = fread(fid,1,'int32');
test = fread(fid,1,'double');
fclose(fid);