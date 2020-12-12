% Save PMDFs into .txt file
% pi pj vk value dove value è quello che leggo nel voxel vk se metto due optodi in pi e pj
% Ovviamente, mi devi scrivere solo le quadruple per cui value >= 1e-9.

pathnameLoad = 'C:\Users\brigadoi\Dropbox\Array_Designer\Fluence_Preprocessing\DS\';
pathnameSave = 'C:\Users\brigadoi\Dropbox\Array_Designer\Fluence_Preprocessing\DS\TestTxt\';

if ~exist(pathnameSave,'dir')
    mkdir(pathnameSave)
end

load(fullfile(pathnameLoad,'ALLPMDFs_DF_2_5_ToastCorrect.mat'));

fid = fopen(fullfile(pathnameSave,'PMDF_indexFile.idx'),'w');
for iR = 1:size(ALLPMDFs,2)
    fprintf('%d \n',iR)
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

fid = fopen(fullfile(pathnameSave,'PMDF_dataFile.data'),'w');
for iR = 1:size(ALLPMDFs,2)
    fprintf('%d \n',iR)
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

tmp = zeros(nTot+nTot1,2);
fid = fopen(fullfile(pathnameSave,'PMDF_dataFile.bin'));
for i = 1:nTot1
    tmp(i,1) = fread(fid,1,'int32');
    tmp(i,2) = fread(fid,1,'*double');
end
for k = 1:nTot
    tmp(i+k,1) = fread(fid,1,'int32');
    tmp(i+k,2) = fread(fid,1,'*double');
end
fclose(fid);

%%
% Save .txt also for GM and Scalp
load('Atlas_Model\10-2p5_Model\refpts_10-2p5.mat');
posScalp = refpts_10_2p5;
load('Atlas_Model\10-2p5_Model\GMDownsampledMesh_DF_2_5.mat');
posGM = GMSurfaceMesh.node(:,1:3);

fid = fopen(fullfile(pathnameSave,'scalpPos.txt'),'w');
for iV = 1:size(posScalp,1)
    fprintf(fid,'%.6f %.6f %.6f \n',posScalp(iV,:));
end
fclose(fid);          

fid = fopen(fullfile(pathnameSave,'GMPos.txt'),'w');
for iV = 1:size(posGM,1)
    fprintf(fid,'%.6f %.6f %.6f \n',posGM(iV,:));
end
fclose(fid);  
