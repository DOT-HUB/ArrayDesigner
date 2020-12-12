function [viableChannelDists] = plotResults_GRASP(ROIpathname,nposGM,minRho,maxRho,GMSurfaceMesh,ScalpSurfaceMesh,posScalp,ALLPMDFs,PMDFweighting, coverageThresh,sources,detectors,time,signal,signalPerc,coveragePerc)

nS = length(sources);
nD = length(detectors);

%Load ROIs
fid = fopen(ROIpathname);
tmp = textscan(fid,'%s',1,'delimiter','\n');
fclose(fid);

tmp = tmp{1};
ROIs = str2num(tmp{1});

ROIMasks = false(1,nposGM);
ROIMasks(ROIs) = true; %Non-downsampled

load('CMAPgreyMag.mat');
load('CMAPgreyGreen.mat');

figure('units','normalized','outerposition',[0 0 1 1],'color','w')

subplot(1,2,1);
displayIntensityOnMesh_RJC(GMSurfaceMesh,double(ROIMasks));caxis([0 1]);
colormap(greyGreen);colorbar off;
view([0 90]);

array = [sources detectors];

if ~isempty(array)
    
    A = GetSensitivityMap(array,posScalp,nS,minRho,maxRho,ALLPMDFs,PMDFweighting);
    %A(A>=coverageThresh) = 1;
    %A(A<coverageThresh) = 0;
    subplot(1,2,2);
    [viableChannelDists] = Plot_Array(array,posScalp,ScalpSurfaceMesh,GMSurfaceMesh,A,nS,minRho,maxRho,0);
    legend off;
    tit = ['[GRASP: ' num2str(length(viableChannelDists)) ' chans, ' num2str(time,'%0.1f')  ' s, ' num2str(signal,'%0.2f') ' (' num2str(signalPerc,'%0.1f') ' %) ,' num2str(coveragePerc,'%0.1f') ' %]'];
    title(tit);
else
    subplot(1,2,2);
    imagesc(-5*ones(10,10));caxis([-15 0]);axis equal;axis off;
    tit = ['[FAIL]'];
    title(tit);
end



