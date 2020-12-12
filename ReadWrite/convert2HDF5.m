% Convert and save data to HDF5 format for Domenico

%load('Atlas_Model\10-2p5_Model\refpts_10-2p5.mat');
%posScalp = refpts_10_2p5;
%load('Atlas_Model\10-2p5_Model\GMDownsampledMesh_DF_2_5.mat');
%posGM = GMSurfaceMesh.node(:,1:3);
%nposGM = length(posGM);

%h5create('posScalp.h5','/posScalp',size(posScalp))
%h5write('posScalp.h5','/posScalp',posScalp);

%h5create('posGM_DF_2_5.h5','/posGM',size(posGM))
%h5write('posGM_DF_2_5.h5','/posGM',posGM);

%load('Fluence_Processing\ALLPMDFs_10-2p5.mat');
load('ALLPMDFs_DF_2_5_ToastCorrect_FORCESYMM.mat');
 
filename = 'ALLPMDFs_DF_2_5_ToastCorrect_FORCESYMM.h5';
for iC = 1:size(ALLPMDFs)
    for iR = 1:size(ALLPMDFs)
        if ~isempty(ALLPMDFs{iR,iC})
            
            label = ['/PMDF_row_' num2str(iR) '_col_' num2str(iC)];
            
            h5create(['Fluence_Preprocessing/DS/' filename],label,size(ALLPMDFs{iR,iC}))
            h5write(['Fluence_Preprocessing/DS/' filename],label,full(ALLPMDFs{iR,iC}));
            
        end
    end
end