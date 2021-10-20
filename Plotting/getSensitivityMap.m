function A = getSensitivityMap(array,posScalp,nS,minRho,maxRho,pathnameHeadModel,pathnameWeights)

%Load PMDF weights
PMDFweighting = load(pathnameWeights);

nD = length(array)-nS;
srcs = array(1:nS);
dets = array(nS+1:end);

%Load PMDFs of the channels of the AD array
pathnameMSHS = getMSHSpath(pathnameHeadModel);
load(pathnameMSHS,'gmSurfaceMesh','-mat');
%load(fullfile(pathnameHeadModel,'gmSurfaceMesh'),'-mat'); Possibly revert depending on handling of downsampling.
ALLPMDFs = readPMDFs(srcs,dets,size(gmSurfaceMesh.node,1),minRho,maxRho,posScalp,pathnameHeadModel);

if ~exist('PMDFweighting','var')
    PMDFweighting = ones(size(ALLPMDFs));
end

%Add channels to diagram
for i = 1:nS
    dist = sqrt(sum((posScalp(dets,:) - repmat(posScalp(srcs(i),:),nD,1)).^2,2));
    for j = 1:length(dist)
        if dist(j) >= minRho && dist(j) < maxRho
            tmpSrc = srcs(i);
            tmpDet = dets(j);
            if tmpSrc>tmpDet %This is because only top right of ALLPMDFs is populated
                tmpPMDF = ALLPMDFs{dets(j),srcs(i)};   
            else
                tmpPMDF = ALLPMDFs{srcs(i),dets(j)};
            end
            if ~exist('A','var')
                A = zeros(size(tmpPMDF));
            end
            %tmpPMDF(tmpPMDF < 1e-9) = 0;
            tmpPMDF = tmpPMDF.*PMDFweighting(srcs(i),dets(j));
            
            try 
                A = A+tmpPMDF;
            catch
                stop = 1;
            end
        end
    end
end

