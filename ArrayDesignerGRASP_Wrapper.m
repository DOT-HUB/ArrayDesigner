% Wrapper for example run of Array Designer

%############################################################
%######################### INPUTS ###########################
% Define paths to Example variables
% Set pathname for scalp pos, GM pos, PMDFs and weights

%% Temporary path handling
[ADpath,~] = fileparts(mfilename('fullpath'));

% Place outputs in /Results
pathnameOutput = [ADpath '/Results'];

%% Inputs
% Use default atlas in this example
pathnameHeadModel = [ADpath '/HeadModels/MNI152_LowRes'];
% Use example ROI from paper
pathnameROI = [ADpath '/HeadModels/MNI152_LowRes/ROIs/ExampleROI_v5_DS_1.txt'];

%% Set example user inputs
AD.inputs.nS = 8;             %number of source optodes
AD.inputs.nD = 8;             %number of detector optodes
AD.inputs.minRhoOpt = 10;     %the minimum physical distance between fibres (i.e. fibre diameter)
AD.inputs.maxGoodRho = 30;    %the longest distance your fNIRS system reliably gets good quality measurements
AD.inputs.maxRho = 60;        %the longest distance your fNIRS system ever gets good quality measurements
AD.inputs.coverageWeight = 1; %a weighting factor for the importance of ROI coverage. 1 = coverage is equally as important as sensitivity.

AD.inputs.pathnameHeadModel = pathnameHeadModel;
AD.inputs.pathnameROI = pathnameROI;
AD.inputs.pathnameOutput = pathnameOutput;

%% Call the GRASP algorithm
if isunix %Mac (or Linux?)
    
    %Run GRASP
    [AD] = ArrayDesignerGRASP_MacLinux(AD,pathnameOutput);
    
    %Save results
    outname = [pathnameOutput '/' datestr(now,'YYYYmmddHHMMss') '.AD'];
    save(outname,'-struct','AD');
    
    %Plot full result
    f1 = plotFullResult(AD);
    
end



