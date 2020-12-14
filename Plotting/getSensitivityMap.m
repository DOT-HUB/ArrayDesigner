function A = getSensitivityMap(array,posScalp,nS,minRho,maxRho,pathnameHeadModel,pathnameWeights)

%Load PMDF weights
PMDFweights = load(pathnameWeights);

%Load PMDFs
PMDFweights = importdata([pathnameHeadModel '/PMDFs/PMDFs.data');

nD = length(array)-nS;
srcs = array(1:nS);
dets = array(nS+1:end);
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

