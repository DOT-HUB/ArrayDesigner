% This function reads the PMDFs of the feasible channels selected by AD from the
% PMDFs.data file and place them in the allPMDFs matrix
% 
% INPUTS: srcs                   sources selected by AD
%         dets                   detectors selected by AD
%         nVoxels                number of GM nodes
%         minRho                 minimun distance between optodes
%         maxRho                 maximum distance between optodes
%         scalpPos               3D positions of all possible scalp positions
%         pathnameHeadModel      pathname to the head model folder where the PMDFs folder is stored 
%
% OUTPUT: allPMDFs               squared cell matrix (N of scalp position) where each row is the source index and each column the detector index. In
%                                each cell of a viable AD channel, its PMDF is stored


function allPMDFs = readPMDFs(srcs,dets,nVoxels,minRho,maxRho,scalpPos,pathnameHeadModel)

% Read idx file - each row contains the number of voxels with significant
% values for each PMDF of the upper triangular matrix
fid = fopen(fullfile(pathnameHeadModel,'PMDFs','PMDFs.idx'),'r');
nPMDFs = fscanf(fid,'%d');
fclose(fid);

% Recover size of square matrix with source-det pairs (upper triangular
% matrix)
sizeMatrix = ceil(sqrt(length(nPMDFs)*2));
nD = length(dets);

% Get index values of the upper triangular matrix  
A = 1:length(nPMDFs);
B = tril(ones(sizeMatrix),-1);
B(B==1)= A;
B = B';

% Create all feasible source-detector couples 
k = 1;
for iR = 1:length(srcs)
    dist = sqrt(sum((scalpPos(dets,:) - repmat(scalpPos(srcs(iR),:),nD,1)).^2,2));
    for iC = 1:length(dist)
        if dist(iC) >= minRho && dist(iC) < maxRho
            if srcs(iR) > dets(iC)
                channels(k,:) = [dets(iC),srcs(iR)];
            else
                channels(k,:) = [srcs(iR),dets(iC)];
            end
            k = k+1;
        end
    end
end
% Order channels from smaller to bigger in terms of source index and then
% detector index
channelsOrd = sortrows(channels);

<<<<<<< HEAD
% Get only PMDFs of selected channels
allPMDFs = cell(sizeMatrix,sizeMatrix);
fid = fopen(fullfile(pathnameHeadModel,'PMDFs','PMDFs.data'),'r');
posVectorPre = 1;
for iCh = 1:size(channelsOrd,1)
    allPMDFs{channelsOrd(iCh,1),channelsOrd(iCh,2)} = zeros(nVoxels,1);
    
    posVector = B(channelsOrd(iCh,1),channelsOrd(iCh,2));
    nBytes = sum(nPMDFs(posVectorPre:posVector-1))*4 + sum(nPMDFs(posVectorPre:posVector-1))*8; % int32 -> 4 bytes, double -> 8 bytes
    fseek(fid, nBytes, 'cof'); % skip these bytes from current position
    
    for iG = 1:nPMDFs(posVector)
        pos = fread(fid,1,'int32');
        allPMDFs{channelsOrd(iCh,1),channelsOrd(iCh,2)}(pos+1) = fread(fid,1,'double');
    end
    
    posVectorPre = posVector+1;
end
fclose(fid);    
=======
%% Example on how to read specific parts of the .data file
fid = fopen(fullfile(pathnameLoad,'PMDFs','PMDFs.data'),'r');
nBytes = nPMDFs(1)*4 + nPMDFs(1)*8; % int32 -> 4 bytes, double -> 8 bytes
fseek(fid, nBytes, 'bof'); % skip these bytes from origin
pos = fread(fid,1,'int32');
test = fread(fid,1,'double');
fclose(fid);



>>>>>>> 3c92d8cb2b6ca88cb7e10b6eec99cf95b41c4c03
