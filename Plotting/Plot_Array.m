function [viableChannelDists] = Plot_Array(array,posScalp,ScalpSurfaceMesh,GMSurfaceMesh,A,nS,minRho,maxRho,scalpFlag)

%Plotting function.  Assumes Array is organised as a list of points in
%lattice 'pos', ordered such that sources = Array(1:nS) and detectors =
%Array(nS+1:end);

nD = length(array)-nS;
%figure;
%Add lattice points and src/det
%Don't plot array points
%ind = true(length(posScalp),1);
%ind(array) = false;
%scatter3(posScalp(ind,1),posScalp(ind,2),posScalp(ind,3),3,[0 0 0]);
%scatter(posScalp(roiMask,1),posScalp(roiMask,2),'g.');
%Plot scalp mesh
%binGreyGreenMap = zeros(64,3);
%binGreyGreenMap(1:end/2,:) = repmat([0.95 0.95 0.95],32,1);
%binGreyGreenMap(end/2+1:end,:) = repmat([0 1 0],32,1);

greyGreen = [0.5    0.5    0.5
    0.4940    0.5099    0.4940
    0.4860    0.5178    0.4860
    0.4781    0.5257    0.4781
    0.4701    0.5336    0.4701
    0.4621    0.5415    0.4621
    0.4542    0.5494    0.4542
    0.4462    0.5573    0.4462
    0.4382    0.5652    0.4382
    0.4303    0.5731    0.4303
    0.4223    0.5810    0.4223
    0.4143    0.5889    0.4143
    0.4063    0.5968    0.4063
    0.3984    0.6047    0.3984
    0.3904    0.6126    0.3904
    0.3824    0.6205    0.3824
    0.3745    0.6284    0.3745
    0.3665    0.6364    0.3665
    0.3585    0.6443    0.3585
    0.3506    0.6522    0.3506
    0.3426    0.6601    0.3426
    0.3346    0.6680    0.3346
    0.3267    0.6759    0.3267
    0.3187    0.6838    0.3187
    0.3107    0.6917    0.3107
    0.3028    0.6996    0.3028
    0.2948    0.7075    0.2948
    0.2868    0.7154    0.2868
    0.2789    0.7233    0.2789
    0.2709    0.7312    0.2709
    0.2629    0.7391    0.2629
    0.2550    0.7470    0.2550
    0.2470    0.7549    0.2470
    0.2390    0.7628    0.2390
    0.2311    0.7707    0.2311
    0.2231    0.7786    0.2231
    0.2151    0.7866    0.2151
    0.2072    0.7945    0.2072
    0.1992    0.8024    0.1992
    0.1912    0.8103    0.1912
    0.1833    0.8182    0.1833
    0.1753    0.8261    0.1753
    0.1673    0.8340    0.1673
    0.1594    0.8419    0.1594
    0.1514    0.8498    0.1514
    0.1434    0.8577    0.1434
    0.1354    0.8656    0.1354
    0.1275    0.8735    0.1275
    0.1195    0.8814    0.1195
    0.1115    0.8893    0.1115
    0.1036    0.8972    0.1036
    0.0956    0.9051    0.0956
    0.0876    0.9130    0.0876
    0.0797    0.9209    0.0797
    0.0717    0.9289    0.0717
    0.0637    0.9368    0.0637
    0.0558    0.9447    0.0558
    0.0478    0.9526    0.0478
    0.0398    0.9605    0.0398
    0.0319    0.9684    0.0319
    0.0239    0.9763    0.0239
    0.0159    0.9842    0.0159
    0.0080    0.9921    0.0080
         0    1.0000         0];

displayIntensityOnMesh_RJC(GMSurfaceMesh,A);caxis([0 1]);
hold on;colormap(greyGreen);colorbar off;

if ~exist('scalpFlag','var');
    scalpFlag = 1;
end
if scalpFlag == 1;
    plotmesh(ScalpSurfaceMesh.node,ScalpSurfaceMesh.face,'FaceColor','none','FaceAlpha',0,'EdgeColor','k','EdgeAlpha',0.3);
end
nchan = 0;
%Add channels to diagram
for i = 1:nS;
    dist = sqrt(sum((posScalp(array(nS+1:end),:) - repmat(posScalp(array(i),:),nD,1)).^2,2));
    for j = 1:length(dist);
        if dist(j) >= minRho && dist(j) < maxRho
            Hl = line([posScalp(array(i),1) posScalp(array(nS+j),1)],[ posScalp(array(i),2) posScalp(array(nS+j),2)],[posScalp(array(i),3) posScalp(array(nS+j),3)],'color','m','LineWidth',2);
            nchan = nchan+1;
            viableChannelDists(nchan) = dist(j);
        end
        hold on;
    end
end

Hs = scatter3(posScalp(array(1:nS),1),posScalp(array(1:nS),2),posScalp(array(1:nS),3),70,'MarkerEdgeColor','k','MarkerFaceColor','r');
Hd = scatter3(posScalp(array(nS+1:end),1),posScalp(array(nS+1:end),2),posScalp(array(nS+1:end),3),70,'MarkerEdgeColor','k','MarkerFaceColor','b');
view([0 90]);
%Add ROI centre
%scatter(roiCent(1),roiCent(2),100,'go','filled');
%axis equal tight
%Hleg = legend([Hs Hd Hl],{'Sources','Detectors','Channels'},'FontSize',14);
%set(Hleg,'Position',[0.8605    0.1564    0.0797    0.0763]);
hold off;

