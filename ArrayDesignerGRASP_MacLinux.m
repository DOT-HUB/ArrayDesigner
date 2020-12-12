function [sources,detectors,signal,signalFrac,coverageThresh,coverageFrac,runtime] = ArrayDesignerGRASP_MacLinux(ADinputs,pathnameOutput)

%Book-keeping
optLocation = '/Users/RCooper/Dropbox/Repositories/ArrayDesigner/nirsCPPMacLinux/build/nirsmain';
pathnameScalpPos = [ADinputs.pathnameHeadModel '/scalpPos.txt'];
pathnameGMPos = [ADinputs.pathnameHeadModel '/GMPos.txt'];
pathnameAvNodalVol = [ADinputs.pathnameHeadModel '/avNodalVol.txt'];
pathnamePMDFidx = [ADinputs.pathnameHeadModel '/PMDFs/PMDF.idx'];
pathnamePMDF = [ADinputs.pathnameHeadModel '/PMDFs/PMDF.data'];
pathnameResults = [pathnameOutput '/results.txt'];

%Calculate PMDF weights (this is system specific and shouldn't change much,
%so we save a file that is reloaded if already exists.
[pathnameWeights] = getPMDFWeights(ADinputs.maxGoodRho,ADinputs.maxRho,ADinputs.pathnameHeadModel);

%% Determine coverage threshold
%Calculated on the basis of the sensitivity needed toyield a percThresh intensity 
%change when a 10% mua change occurs in a cubic cm
avNodalVol = load(pathnameAvNodalVol);  %Average voronoi volume of node mapped to GM node - this is a function of the mesh used
activationVol = 10^3;                   %Activation volume
deltaMua = 0.001;                       %Mua change
percThresh = 1;                         %Threshold percentage signal change
coverageThresh = log((100+percThresh)/100)/((activationVol/avNodalVol)*deltaMua);

%% Run GRASP
disp('Running GRASP algorithm...');

tic
[~,b] = system([optLocation ' ' pathnameScalpPos ' ' pathnameGMPos ' ' pathnamePMDFidx ' ' pathnamePMDF ' ' pathnameWeights ' ' ADinputs.pathnameROI ' ' num2str(ADinputs.nS) ' ' num2str(ADinputs.nD) ' ' num2str(ADinputs.minRhoOpt) ' ' num2str(coverageThresh) ' ' num2str(ADinputs.coverageWeight) ' ' pathnameResults]);
runtime = toc;
disp(b);

% Load results are in /tmp/name of ROI
all = importdata(pathnameResults);
sources = str2num(all.textdata{1}(2:end-1));
detectors = str2num(all.textdata{2}(2:end-1));
signal = all.data(1,:);
signalFrac = all.data(2,:);
coverage = all.data(3,:);
coverageFrac = all.data(4,:);

