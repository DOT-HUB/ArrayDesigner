function [hAxes] = plotChanHist(AD,varargin)

% Manage Inputs ############################### 
varInputs = inputParser;
addParameter(varInputs,'hAxes','',@ishandle);
parse(varInputs,varargin{:});
varInputs = varInputs.Results;
if isempty(varInputs.hAxes)
    hAxes = gca;
else
    hAxes = varInputs.hAxes;
end

%Add channel histogram
hh1 = histogram(hAxes,AD.results.channelDists,'FaceColor',[0 0 1],'binwidth',2.5);
hAxes.YLimMode = 'auto';
xlim(hAxes,[AD.inputs.minRho AD.inputs.maxRho]);
set(hAxes, 'XLimSpec', 'Tight','YLimSpec','Tight');