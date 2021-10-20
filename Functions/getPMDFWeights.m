function [pathnameWeights] = getPMDFWeights(maxGoodRho,maxRho,pathnameHeadModel)

%##########################################################################
%Create weighting values for PMDFs on the basis of maxGoodRho, maxRho.
%We want channels that are too long to be penalized, rather than
%forcibly excluded.  This is more suitably for the linear integer solver.

%pathnameHeadModel will contain scalpPos.txt and GMpos.txt and the PMDFs in
%some organised format.

% DEBUG
%pathnameHeadModel = 'HeadModels/MNI152_LowRes';
%maxGoodRho = 30;
%maxRho = 45;

%%
[~,HeadModelName,~] = fileparts(pathnameHeadModel);
pathnameWeights = fullfile(pathnameHeadModel,'PMDFweightings',['PMDFweighting_' num2str(maxGoodRho) '-' num2str(maxRho) 'mm_' HeadModelName '.txt']);

if ~exist(fullfile(pathnameHeadModel,'PMDFweightings'),'dir')
    mkdir(fullfile(pathnameHeadModel,'PMDFweightings'));
end

if ~isfile(pathnameWeights)
    disp('Producing PDMF weights file');
    % Load PMDF Maxima
    % PMDFMaxima = load([pathnameHeadModel '/PMDFMaxima.txt']);
    %% TEMPORARY hack until PMDF filetypes are sorted.
    %load('/Users/RCooper/Dropbox/Projects/Array_Designer/Fluence_Preprocessing/DS/ALLPMDFs_DF_2_5_ToastCorrect.mat','PMDFnormFactor');
    load(fullfile(pathnameHeadModel,'PMDFs','PMDFNormFactor.mat'))
    PMDFMaxima = PMDFnormFactor; clear PMDFnormFactor
    %%
    
    % Load Scalp Pos
    fid = fopen(fullfile(pathnameHeadModel,'Utils','scalpSolutionSpace_10_2p5.txt'),'r');
    tmp = textscan(fid,'%f %f %f');
    fclose(fid);
    scalpPos = [tmp{1} tmp{2} tmp{3}];
    nScalpPos = size(scalpPos,1);
    
    %Create distMAT
    distMat = zeros(nScalpPos);
    for ii = 1:nScalpPos
        for jj = 1:nScalpPos
            dist = sqrt(sum((scalpPos(ii,:) - scalpPos(jj,:)).^2,2));
            distMat(ii,jj) = dist;
        end
    end
    
    %First find average pmdfMAX value at maxGoodDist +/1 5%;
    pmdfNormMean = mean(PMDFMaxima(distMat < maxGoodRho*1.05 & distMat > maxGoodRho*0.95));
    
    %Now, mask and weight PMDFs appropriately
    PMDF_MASK = ones(nScalpPos);
    PMDF_MASK(1:nScalpPos+1:end) = 0; %Diagonal values zero;
    
    wb = waitbar(0,'Computing raw weightings...');
    for ii = 1:nScalpPos
        waitbar(ii/nScalpPos,wb);
        for jj = 1:nScalpPos
            
            dist = distMat(ii,jj);
            
            %First -  - scale PMDF if it is larger than optimum separation due to
            %system SNR.  Use ration of pmdfMAX, which drops off with the same
            %rate as optical intensity (I think).
            if dist > maxGoodRho
                rawweight = PMDFMaxima(ii,jj)/pmdfNormMean;
                PMDF_MASK(ii,jj) = rawweight;
            end
            
            %Second set PMDF to empty if it falls outside distance limits
            if  dist > maxRho
                PMDF_MASK(ii,jj) = 0;
            end
            
        end
    end
    delete(wb)
    
    %Fitting
    m1 = distMat(:) > maxGoodRho & distMat(:) < maxRho;
    m2 = PMDF_MASK(:) > 0;
    m = m1 & m2;
    distMat_toplot = distMat(m);
    PMDFMASK_toplot = log10(PMDF_MASK(m));
    p = polyfit(distMat_toplot,PMDFMASK_toplot,1);
    
    %f1 = figure;
    %subplot(1,2,1);
    %scatter(distMat_toplot,PMDFMASK_toplot);xlim([0 60]);
    %hold on;
    %plot(maxGoodRho:maxRho,[maxGoodRho:maxRho].*p(1) + p(2),'LineWidth',2);
    %xlabel('SD Distance (mm)','FontSize',14);
    %ylabel('Log10(weighting)','FontSize',14);
    %legend('Raw relative PMDF norms','Linear Fit');
    %set(gca,'FontSize',14);
    %fprintf(['Gradient = ' num2str(p(1)) '\n']);
    
    %Create weighting from fit;
    %Now, mask and weight PMDFs appropriately
    %PMDF_MASKfit = ones(nScalpPos);
    %PMDF_MASKfit(1:nScalpPos+1:end) = 0; %Diagonal values zero;
    wb = waitbar(0,'Computing fitted weightings...');
    for ii = 1:nScalpPos
        waitbar(ii/nScalpPos,wb);
        for jj = 1:nScalpPos
            
            dist = distMat(ii,jj);
            
            %First -  - scale PMDF if it is larger than optimum separation due to
            %system SNR.  Use ration of pmdfMAX, which drops off with the same
            %rate as optical intensity (I think).
            if dist > maxGoodRho
                weight = 10^(p(1)*(dist - maxGoodRho));
                PMDF_MASKfit(ii,jj) = weight;
            end
            
            %Second set PMDF to empty if it falls outside distance limits
            if  dist > maxRho
                PMDF_MASKfit(ii,jj) = 0;
            end
            
            % Add by SB to remove points on SD=0 in the figure;
% %             if dist == 0
% %                 PMDF_MASKfit(ii,jj) = NaN;
% %                 PMDF_MASK(ii,jj) = NaN;
% %             end
        end
    end
    delete(wb)
    
    %Compare functions;
    %figure(f1);
    %subplot(1,2,2);
    %scatter(distMat(:),PMDF_MASK(:));xlim([0 60]);hold on;
    %scatter(distMat(:),PMDF_MASKfit(:));xlim([0 60]);
    %ylabel('SNR weighting factor','FontSize',14);
    %xlabel('Source-detector distance [mm]','FontSize',14);
    %set(gca,'FontSize',14);
    %set(gcf,'PaperPositionMode','auto','Position',[560 483 733 465])
    %drawnow;
    
    %% Save out
    disp(['Writing PMDF weighting file to ' pathnameWeights '...']);
    writematrix(PMDF_MASKfit,pathnameWeights,'delimiter','\t');
    
else
    disp(['PDMF weights file found: ' pathnameWeights ', loading...']);
end

