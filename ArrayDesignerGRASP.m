function [AD] = ArrayDesignerGRASP(AD)

%Book-keeping
[ADpath,~] = fileparts(mfilename('fullpath'));
pathnameScalpPos = fullfile(AD.inputs.pathnameHeadModel,'scalpPos.txt');
pathnameGMPos = fullfile(AD.inputs.pathnameHeadModel,'GMPos.txt');
pathnameAvNodalVol = fullfile(AD.inputs.pathnameHeadModel,'avNodalVol.txt');
pathnamePMDFidx = fullfile(AD.inputs.pathnameHeadModel,'PMDFs','PMDFs.idx');
pathnamePMDF = fullfile(AD.inputs.pathnameHeadModel,'PMDFs','PMDFs.data');
pathnameResults = fullfile(AD.inputs.pathnameOutput,'results.txt');

%Calculate PMDF weights (this is system specific and shouldn't change much,
%so we save a file that is reloaded if already exists).
[pathnameWeights] = getPMDFWeights(AD.inputs.maxGoodRho,AD.inputs.maxRho,AD.inputs.pathnameHeadModel);
AD.inputs.pathnameWeights = pathnameWeights;

%% Determine coverage threshold
%Calculated on the basis of the sensitivity needed toyield a percThresh intensity 
%change when a 10% mua change occurs in a cubic cm
avNodalVol = load(pathnameAvNodalVol);  %Average voronoi volume of node mapped to GM node - this is a function of the mesh used
activationVol = 10^3;                   %Activation volume
deltaMua = 0.001;                       %Mua change
percThresh = 1;                         %Threshold percentage signal change
coverageThresh = log((100+percThresh)/100)/((activationVol/avNodalVol)*deltaMua);
AD.inputs.coverageThresh = coverageThresh;

%% Create status box
% h = figure('Units','Normalized');
% h.MenuBar = 'none';
% h.ToolBar = 'none';
% h.Name = 'Array Designer Status:';
txtBox = uitextarea;
hFig = txtBox.Parent;
txtBox.Position = [5 5 hFig.Position(3)-10 hFig.Position(4)-10];
drawnow

%% Run GRASP
statString = {'Running Array Designer GRASP algorithm...'};
txtBox.Value = statString;
drawnow

tic
if isunix
    optLocation = fullfile(ADpath,'nirsCPPMacLinux','build','nirsmain');
    [~,b] = system([optLocation ' ' pathnameScalpPos ' ' pathnameGMPos ' ' pathnamePMDFidx ' ' pathnamePMDF ' ' pathnameWeights ' ' AD.inputs.pathnameROI ' ' num2str(AD.inputs.nS) ' ' num2str(AD.inputs.nD) ' ' num2str(AD.inputs.minRhoOpt) ' ' num2str(coverageThresh) ' ' num2str(AD.inputs.coverageWeight) ' ' pathnameResults]);
else
    %windows call
end
runtime = toc;
statString{end+1} = b;
txtBox.Value = statString;
drawnow

% Load results are in /tmp/name of ROI
all = importdata(pathnameResults);
sources = str2num(all.textdata{1}(2:end-1));
detectors = str2num(all.textdata{2}(2:end-1));
signal = all.data(1,:);
signalFrac = all.data(2,:);
coverage = all.data(3,:);
coverageFrac = all.data(4,:);

AD.results.source = sources;
AD.results.detector = detectors;
AD.results.signal = signal;
AD.results.signalPerc= signalFrac*100;
AD.results.coveragePerc = coverageFrac*100;
AD.results.runtime = runtime;

% Load pos scalp
scalpPos = importdata(fullfile(AD.inputs.pathnameHeadModel,'scalpPos.txt'));
% Calculate Sensitivity
AD.results.sensitivityMap = getSensitivityMap([AD.results.source AD.results.detector],scalpPos,AD.inputs.nS,AD.inputs.minRho,AD.inputs.maxRho,AD.inputs.pathnameHeadModel,AD.inputs.pathnameWeights);


