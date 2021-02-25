function [hAxes] = plotADOutputText(AD,varargin)

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

output{1} = ['Coverage = ' num2str(AD.results.coveragePerc,'%0.1f') '%'];
output{2} = ['Rel. sensitivity = ' num2str(AD.results.signalPerc,'%0.1f') '%'];
output{3} = ['Sensitivity = ' num2str(AD.results.signal,'%0.1f') ' mm'];
output{4} = ['nChannels = ' num2str(length(viableChannelDists))];
output{5} = ['Time = ' num2str(AD.results.runtime,'%0.1f') ' s'];
