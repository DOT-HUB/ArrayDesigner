function plotArrayFileResults(AD)

if isstring(AD)
    AD = load(ADfilename,'-mat');
end

nS = AD.inputs.nS;
nD = AD.inputs.nD;
minRhoOpt = AD.inputs.minRhoOpt;
maxRho = AD.inputs.maxRho;

% %############################################################
%Head Model info
%Load GM and scalp pos
load([pathnameHeadModel '/GMSurfaceMesh.mat']);
load([pathnameHeadModel '/refpts_10-2p5.mat']);
load([pathnameHeadModel '/refpts_10-2p5.mat']);

GMSurfaceMesh.node = GMSurfaceMesh.node(:,1:3);
GMSurfaceMesh.face = GMSurfaceMesh.face(:,1:3);
posScalp = refpts_10_2p5;
posGM = GMSurfaceMesh.node(:,1:3);
nposGM = length(posGM);

if ~exist('ALLPMDFs','var')
    load('/Users/RCooper/Dropbox/Array_Designer/Fluence_Preprocessing/DS/ALLPMDFs_DF_2_5_ToastCorrect.mat')
end


% Plotting
PMDFweighting = ones(size(ALLPMDFs));

load('CMAPgreyMag.mat');
load('CMAPgreyGreen.mat');

array = [sources detectors];

A = GetSensitivityMap(array,posScalp,nS,minRhoOpt,maxRho,ALLPMDFs,PMDFweighting);
A(A>=coverageThresh) = 1;
A(A<coverageThresh) = 0;
[viableChannelDists] = Plot_Array(array,posScalp,ScalpSurfaceMesh,GMSurfaceMesh,A,nS,minRhoOpt,maxRho,0);
legend off;
%tit = ['[GRASP: ' num2str(length(viableChannelDists)) ' chans, ' num2str(time,'%0.1f')  ' s, ' num2str(signal,'%0.2f') ' (' num2str(signalPerc,'%0.1f') ' %) ,' num2str(coveragePerc,'%0.1f') ' %]'];
%title(tit);

subplot(6,7,39);
hist(viableChannelDists*headRatio,30);
xlabel('SD dist (mm)');
ylabel('# chans');

