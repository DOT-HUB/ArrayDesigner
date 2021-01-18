function f1 = plotFullResult(AD)

% Manage Inputs ############################### 
if isstr(AD) %path to .AD file
    AD = load(AD,'-mat');
end 

%delete me!
AD.inputs.minRho = 15;

% Plotting ###############################
f1 = figure('Units','Normalized','Position',[0 0 1 1],'Color','w');

subplot(1,2,1);
hAxis1 = plotROI(AD.inputs.pathnameHeadModel,AD.inputs.pathnameROI);
t1 = title('ROI','FontSize',20);
p1 = get(hAxis1,'Position');

% Load pos scalp
scalpPos = importdata(fullfile(AD.inputs.pathnameHeadModel,'scalpPos.txt'));

A = getSensitivityMap([AD.results.source AD.results.detector],scalpPos,AD.inputs.nS,AD.inputs.minRho,AD.inputs.maxRho,AD.inputs.pathnameHeadModel,AD.inputs.pathnameWeights);

subplot(1,2,2);
[hAxis2, viableChannelDists] = plotArray(AD,A);
t2 = title('Array Designer Solution','FontSize',20);
p2 = get(hAxis2,'Position');

%Set ROI axis to match array axis
ylim(hAxis1,hAxis2.YLim);
xlim(hAxis1,hAxis2.XLim);

%Add solution information
hAxis3 = axes('Position',[p1(1)+p1(3),0.25,p2(1) - (p1(1)+p1(3)),0.2],'visible','off');
x = 0.01;
y = 0.01;
voff = 0.1;
fSize = 11;

tit = ['Coverage = ' num2str(AD.results.coveragePerc,'%0.1f') '%'];
tmpH = text(1,1,tit);
set(tmpH,'Position',[x y 0],'FontSize',fSize);

tit = ['Rel. sensitivity = ' num2str(AD.results.signalPerc,'%0.1f') '%'];
tmpH = text(1,1,tit);
set(tmpH,'Position',[x y+voff 0],'FontSize',fSize);

tit = ['Sensitivity = ' num2str(AD.results.signal,'%0.1f') ' mm'];
tmpH = text(1,1,tit);
set(tmpH,'Position',[x y+2*voff 0],'FontSize',fSize);

tit = ['Time = ' num2str(AD.results.runtime,'%0.1f') ' s'];
tmpH = text(1,1,tit);
set(tmpH,'Position',[x y+3*voff 0],'FontSize',fSize);

tit = ['Channels = ' num2str(length(viableChannelDists))];
tmpH = text(1,1,tit);
set(tmpH,'Position',[x y+4*voff 0],'FontSize',fSize);

%Add channel histogram
hAxis4 = axes;
width = 0.1;
height = 0.1;
set(hAxis4,'position',[0.5-width/2 0.1 width height]);
axes(hAxis4);
hh1 = histogram(viableChannelDists,6,'FaceColor',[0 0 1],'binwidth',5);
axis tight;
xlabel('SD distance (mm)');
ylabel('# channels');
set(hAxis4,'YAxisLocation','right','XTick',[15 30 45]);
xlim([AD.inputs.minRho AD.inputs.maxRho]);

%Add colorbar
cb = colorbar(hAxis2);
width = 0.015;
height = 0.1;
cb.Position = [0.915, 0.1, width, height];
cb.Ticks = [1 2 3];
cb.TickLabels = {'Uncovered','Additional','Covered'};
cb.FontSize = fSize;

