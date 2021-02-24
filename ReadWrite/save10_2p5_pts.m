% This function computes the 10-2.5 locations on the head surface
% 
% INPUTS: 
%
% mshs                    :  A file containing a structure of:
%
%                           % headVolumeMesh     :   The multi-layer volume mesh structure. 
%                                                    Contains fields: node, face, elem, labels
%
%                           % gmSurfaceMesh      :   The gm surface mesh structure. 
%                                                    Contains fields: node, face.
%
%                           % scalpSurfaceMesh   :   The scalp surface mesh structure.
%                                                    Contains fields: node, face.
%
%                           % vol2gm             :   The sparse matrix mapping from head volume mesh
%                                                    space to GM surface mesh space
%
%                           % landmarks          :   A matrix containing the
%                                                    landmarks coordinate
%
%                           % tenFive            :   ten-five locations for
%                                                    the mesh (.positions
%                                                    (nx3) and .labels {n})
%
%                           % logData            :   optional
%
%                           % fileName           :   The path of the saved mshs file
%
%
% OUTPUTS:
% 
% mshs                      :  A file containing a structure of:
%
%                           % headVolumeMesh     :   The multi-layer volume mesh structure. 
%                                                    Contains fields: node, face, elem, labels
%
%                           % gmSurfaceMesh      :   The gm surface mesh structure. 
%                                                    Contains fields: node, face.
%
%                           % scalpSurfaceMesh   :   The scalp surface mesh structure.
%                                                    Contains fields: node, face.
%
%                           % vol2gm             :   The sparse matrix mapping from head volume mesh
%                                                    space to GM surface mesh space
%
%                           % landmarks          :   A matrix containing the
%                                                    landmarks coordinate
%
%                           % tenFive            :   ten-five locations for
%                                                    the mesh (.positions
%                                                    (nx3) and .labels {n})
%
%                           % ADSolutionSpace    :   ten-two.five locations for
%                                                    the mesh (.positions
%                                                    (nx3) and .labels {n})
%
%                           % logData            :   optional
%
%                           % fileName           :   The path of the saved mshs file
%

function mshs = save10_2p5_pts(mshs)

% Upsample scalpSurfaceMesh (which should already be at high resolution)
if ~isfield(mshs,'scalpSurfaceMeshUpsampled')
    mshs.scalpSurfaceMeshUpsampled.node = meshUpsample(mshs.scalpSurfaceMesh.node,mshs.scalpSurfaceMesh.face,1);
end

% Reorder refpts_10_5 to make things easier
refpts_10_5_reorder(1:9,:) = mshs.tenFive.positions(2:10,:);
refpts_10_5_reorder(10,:) = mshs.tenFive.positions(1,:); %Cz
refpts_10_5_reorder(11:length(mshs.tenFive.positions),:) = mshs.tenFive.positions(11:end,:);

refpts_10_5_reorder_label(1:9,:) = mshs.tenFive.labels(2:10,:);
refpts_10_5_reorder_label(10,:) = mshs.tenFive.labels(1,:); %Cz
refpts_10_5_reorder_label(11:length(mshs.tenFive.labels),:) = mshs.tenFive.labels(11:end,:);

refpts_10_5 = refpts_10_5_reorder;
refpts_10_5_label = refpts_10_5_reorder_label;

% Create 10-2.5 coordinates %%%%%%%%%%%%%%%%%
refpts_10_2p5 = refpts_10_5;
refpts_10_2p5_label = refpts_10_5_label;
refpts_10_2p5_tmp = refpts_10_5;
refpts_10_2p5_label_tmp = refpts_10_5_label;

Midline_Points = {'NFpz'    'Fpz'    'AFpz'    'AFz'    'AFFz'    'Fz'    'FFCz'    'FCz'    'FCCz'    'Cz'    'CCPz'    'CPz'    'CPPz'    'Pz'    'PPOz'    'POz'    'POOz'    'Oz'};
%Coronal points to bisect (keep left and right separate because of ordering
%around midline
Left_Coronal_points = {
    
'AF7h'    'AFF7h'    'F7h'    'FFT7h'    'FT7h'    'FTT7h'  'T7h'  'TTP7h'  'TP7h'    'TPP7h'    'P7h'    'PPO7h'    'PO7h'
'AF5'     'AFF5'     'F5'     'FFC5'     'FC5'     'FCC5'   'C5'   'CCP5'   'CP5'     'CPP5'     'P5'     'PPO5'     'PO5'
'AF5h'    'AFF5h'    'F5h'    'FFC5h'    'FC5h'    'FCC5h'  'C5h'  'CCP5h'  'CP5h'    'CPP5h'    'P5h'    'PPO5h'    'PO5h'
'AF3'     'AFF3'     'F3'     'FFC3'     'FC3'     'FCC3'   'C3'   'CCP3'   'CP3'     'CPP3'     'P3'     'PPO3'     'PO3'
'AF3h'    'AFF3h'    'F3h'    'FFC3h'    'FC3h'    'FCC3h'  'C3h'  'CCP3h'  'CP3h'    'CPP3h'    'P3h'    'PPO3h'    'PO3h'
'AF1'     'AFF1'     'F1'     'FFC1'     'FC1'     'FCC1'   'C1'   'CCP1'   'CP1'     'CPP1'     'P1'     'PPO1'     'PO1'
'AF1h'    'AFF1h'    'F1h'    'FFC1h'    'FC1h'    'FCC1h'  'C1h'  'CCP1h'  'CP1h'    'CPP1h'    'P1h'    'PPO1h'    'PO1h'
'AFz'     'AFFz'     'Fz'     'FFCz'     'FCz'     'FCCz'   'Cz'   'CCPz'   'CPz'     'CPPz'     'Pz'     'PPOz'     'POz'};

Right_Coronal_points = {
    'AFz'     'AFFz'     'Fz'     'FFCz'     'FCz'     'FCCz'   'Cz'   'CCPz'   'CPz'     'CPPz'     'Pz'     'PPOz'     'POz'
    'AF2h'    'AFF2h'    'F2h'    'FFC2h'    'FC2h'    'FCC2h'  'C2h'  'CCP2h'  'CP2h'    'CPP2h'    'P2h'    'PPO2h'    'PO2h'
    'AF2'     'AFF2'     'F2'     'FFC2'     'FC2'     'FCC2'   'C2'   'CCP2'   'CP2'     'CPP2'     'P2'     'PPO2'     'PO2'
    'AF4h'    'AFF4h'    'F4h'    'FFC4h'    'FC4h'    'FCC4h'  'C4h'  'CCP4h'  'CP4h'    'CPP4h'    'P4h'    'PPO4h'    'PO4h'
    'AF4'     'AFF4'     'F4'     'FFC4'     'FC4'     'FCC4'   'C4'   'CCP4'   'CP4'     'CPP4'     'P4'     'PPO4'     'PO4'
    'AF6h'    'AFF6h'    'F6h'    'FFC6h'    'FC6h'    'FCC6h'  'C6h'  'CCP6h'  'CP6h'    'CPP6h'    'P6h'    'PPO6h'    'PO6h'
    'AF6'     'AFF6'     'F6'     'FFC6'     'FC6'     'FCC6'   'C6'   'CCP6'   'CP6'     'CPP6'     'P6'     'PPO6'     'PO6'
    'AF8h'    'AFF8h'    'F8h'    'FFT8h'    'FT8h'    'FTT8h'  'T8h'  'TTP8h'  'TP8h'    'TPP8h'    'P8h'    'PPO8h'    'PO8h'};

% These loops find the midway points between consecutive coronal lines
% running from the T7 line up to midline on both left and right hemispsheres

for j = 1:size(Left_Coronal_points,2)
    count = 1;
    for i = 1:size(Left_Coronal_points,1)-1
        
        ind1 = find(strcmpi(Left_Coronal_points(i,j),refpts_10_2p5_label));
        ind2 = find(strcmpi(Left_Coronal_points(i+1,j),refpts_10_2p5_label));
        
        tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
        qpoints(count,:) = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
        
        %Place new point in sequence in refpts_10_2p5
        refpts_10_2p5(ind1+1,:) = qpoints(count,:);
        refpts_10_2p5(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:)),:) = refpts_10_2p5_tmp(ind1+1:end,:);
        
        refpts_10_2p5_label(ind1+1) = strcat(refpts_10_2p5_label(ind1,1), 'q');
        refpts_10_2p5_label(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:))) = refpts_10_2p5_label_tmp(ind1+1:end);
        
        refpts_10_2p5_tmp =  refpts_10_2p5;
        refpts_10_2p5_label_tmp =  refpts_10_2p5_label;
        
        count = count + 1;
    end
    
    Left_Coronal_qpoints(j,:,:) = qpoints;
end

for j = 1:size(Right_Coronal_points,2)
    count = 1;
    for i = 1:size(Right_Coronal_points,1)-1
        
        ind1 = find(strcmpi(Right_Coronal_points(i,j),refpts_10_2p5_label));
        ind2 = find(strcmpi(Right_Coronal_points(i+1,j),refpts_10_2p5_label));
        
        tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
        qpoints(count,:) = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
        
        %Place new point in sequence in refpts_10_2p5
        refpts_10_2p5(ind2,:) = qpoints(count,:);
        refpts_10_2p5(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:)),:) = refpts_10_2p5_tmp(ind2:end,:);
        
        refpts_10_2p5_label(ind2) = strcat(refpts_10_2p5_label(ind2,1), 'q');
        refpts_10_2p5_label(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:))) = refpts_10_2p5_label_tmp(ind2:end);
        
        refpts_10_2p5_tmp =  refpts_10_2p5;
        refpts_10_2p5_label_tmp =  refpts_10_2p5_label;
        
        count = count + 1;
    end
    
    Right_Coronal_qpoints(j,:,:) = qpoints;
end

n = size(Left_Coronal_points,1)+size(Left_Coronal_qpoints,2)-1; %number of total points now along coronal line;

%Now find points half way between coronal lines in front-back direction
for j = 1:size(Left_Coronal_points,2)-1
    qpoints = [];
    ind_first_line1 = find(strcmpi(Left_Coronal_points(1,j),refpts_10_2p5_label));
    ind_first_line2 = find(strcmpi(Left_Coronal_points(1,j+1),refpts_10_2p5_label));
    
    for i = 1:n
        tmp = refpts_10_2p5(ind_first_line1+(i-1),:) + (refpts_10_2p5(ind_first_line2+(i-1),:)-refpts_10_2p5(ind_first_line1+(i-1),:))./2;
        qpoints(i,:) = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
        
        %Get label for new points
        tmp = Get_10_2p5_label(refpts_10_2p5_label{ind_first_line1+(i-1)},refpts_10_2p5_label{ind_first_line2+(i-1)});
        if strcmp('ERR',tmp)
            pause;
        end
        qlabels{i} = tmp;
    end
    
    %Place new point in sequence in refpts_10_2p5
    %This will need to insert whole new coronal line
    
    refpts_10_2p5(end+1:end+n,:) = qpoints; %Just place at end of sequence
    
    refpts_10_2p5_label(end+1:end+n) = qlabels;
    
    refpts_10_2p5_tmp =  refpts_10_2p5;
    refpts_10_2p5_label_tmp =  refpts_10_2p5_label;
end

n = size(Right_Coronal_points,1)+size(Right_Coronal_qpoints,2)-1; %number of total points now along coronal line;

%Now find points half way between coronal lines in front-back direction
for j = 1:size(Right_Coronal_points,2)-1
    qpoints = [];
    ind_first_line1 = find(strcmpi(strcat(Right_Coronal_points(2,j),'q'),refpts_10_2p5_label));
    ind_first_line2 = find(strcmpi(strcat(Right_Coronal_points(2,j+1),'q'),refpts_10_2p5_label));
    
    for i = 1:n
        tmp = refpts_10_2p5(ind_first_line1+(i-1),:) + (refpts_10_2p5(ind_first_line2+(i-1),:)-refpts_10_2p5(ind_first_line1+(i-1),:))./2;
        qpoints(i,:) = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
        
        %Get label for new points
        qlabels{i} = Get_10_2p5_label(refpts_10_2p5_label{ind_first_line1+(i-1)},refpts_10_2p5_label{ind_first_line2+(i-1)});
    end
    
    %Place new point in sequence in refpts_10_2p5
    %This will need to insert whole new coronal line
    
    refpts_10_2p5(end+1:end+n,:) = qpoints; %Just place at end of sequence
    
    refpts_10_2p5_label(end+1:end+n) = qlabels;
    
    refpts_10_2p5_tmp =  refpts_10_2p5;
    refpts_10_2p5_label_tmp =  refpts_10_2p5_label;
end

%Now add additional midline points
for i = 1:length(Midline_Points)-1
    ind1 = find(strcmpi(Midline_Points(i),refpts_10_2p5_label));
    ind2 = find(strcmpi(Midline_Points(i+1),refpts_10_2p5_label));
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    %Place new point in sequence in refpts_10_2p5
    refpts_10_2p5(ind1+1,:) = qpoint;
    refpts_10_2p5(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind2:end,:)),:) = refpts_10_2p5_tmp(ind1+1:end,:);
    
    refpts_10_2p5_label{ind1+1} = Get_10_2p5_label(refpts_10_2p5_label{ind1},refpts_10_2p5_label{ind2});
    refpts_10_2p5_label(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind2:end,:))) = refpts_10_2p5_label_tmp(ind2:end);
    
    refpts_10_2p5_tmp =  refpts_10_2p5;
    refpts_10_2p5_label_tmp =  refpts_10_2p5_label;
end

% Add points to lower crown ring in quarters
Lower_Ring_Left_points = {'NFpz'
    'NFp1h'
    'NFp1'
    'AFp9h'
    'AF9h'
    'AFF9h'
    'F9h'
    'FFT9h'
    'FT9h'
    'FTT9h'
    'T9h'
    'TTP9h'
    'TP9h'
    'TPP9h'
    'P9h'
    'PPO9h'
    'PO9h'
    'POO9h'
    'OI1'
    'OI1h'
    'OIz'};


Lower_Ring_Right_points = {'NFpz'
    'NFp2h'
    'NFp2'
    'AFp10h'
    'AF10h'
    'AFF10h'
    'F10h'
    'FFT10h'
    'FT10h'
    'FTT10h'
    'T10h'
    'TTP10h'
    'TP10h'
    'TPP10h'
    'P10h'
    'PPO10h'
    'PO10h'
    'POO10h'
    'OI2'
    'OI2h'
    'OIz'};

Upper_Ring_Left_points = {'Fpz'
    'Fp1h'
    'Fp1'
    'AFp7'
    'AF7'
    'AFF7'
    'F7'
    'FFT7'
    'FT7'
    'FTT7'
    'T7'
    'TTP7'
    'TP7'
    'TPP7'
    'P7'
    'PPO7'
    'PO7'
    'POO7'
    'O1'
    'O1h'
    'Oz'};


Upper_Ring_Right_points = {'Fpz'
    'Fp2h'
    'Fp2'
    'AFp8'
    'AF8'
    'AFF8'
    'F8'
    'FFT8'
    'FT8'
    'FTT8'
    'T8'
    'TTP8'
    'TP8'
    'TPP8'
    'P8'
    'PPO8'
    'PO8'
    'POO8'
    'O2'
    'O2h'
    'Oz'};

cross_points = [1 11];
% interpolate along lines
for i = 1:length(Lower_Ring_Left_points)-1
    
    ind1 = find(strcmpi(Lower_Ring_Left_points(i),refpts_10_2p5_label));
    ind2 = find(strcmpi(Lower_Ring_Left_points(i+1),refpts_10_2p5_label));
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    
    %Place new point in sequence in refpts_10_2p5
    if ismember(i,cross_points) %If first point (ind1) is a cross point, then point should be index ind2-1
        refpts_10_2p5(ind2,:) = qpoint;
        refpts_10_2p5(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:)),:) = refpts_10_2p5_tmp(ind2:end,:);
        
        if i == 1
            tmp2 = refpts_10_2p5_label{ind2};
            refpts_10_2p5_label{ind2} = strcat(tmp2(1:end-1),'q');
            refpts_10_2p5_label(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:))) = refpts_10_2p5_label_tmp(ind2:end);            
        else
            refpts_10_2p5_label{ind2} = strcat(refpts_10_2p5_label{ind2},'q');
            refpts_10_2p5_label(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:))) = refpts_10_2p5_label_tmp(ind2:end);
        end
    else
        refpts_10_2p5(ind1+1,:) = qpoint;
        refpts_10_2p5(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:)),:) = refpts_10_2p5_tmp(ind1+1:end,:);
        
        refpts_10_2p5_label{ind1+1,:} = Get_10_2p5_label(refpts_10_2p5_label{ind1},refpts_10_2p5_label{ind2});
        refpts_10_2p5_label(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:))) = refpts_10_2p5_label_tmp(ind1+1:end);
        
    end
    
    refpts_10_2p5_tmp =  refpts_10_2p5;
    refpts_10_2p5_label_tmp =  refpts_10_2p5_label;
end

for i = 1:length(Lower_Ring_Right_points)-1
    ind1 = find(strcmpi(Lower_Ring_Right_points(i),refpts_10_2p5_label));
    ind2 = find(strcmpi(Lower_Ring_Right_points(i+1),refpts_10_2p5_label));
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    %Place new point in sequence in refpts_10_2p5
    if ismember(i,cross_points) %If first point (ind1) is a cross point, then point should be index ind2-1
        refpts_10_2p5(ind2,:) = qpoint;
        refpts_10_2p5(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:)),:) = refpts_10_2p5_tmp(ind2:end,:);
        
        if i == 1
            tmp2 = refpts_10_2p5_label{ind2};
            refpts_10_2p5_label{ind2} = strcat(tmp2(1:end-1),'q');
            refpts_10_2p5_label(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:))) = refpts_10_2p5_label_tmp(ind2:end);            
        else
            refpts_10_2p5_label{ind2} = strcat(refpts_10_2p5_label{ind2},'q');
            refpts_10_2p5_label(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:))) = refpts_10_2p5_label_tmp(ind2:end);
        end
        
    else
        refpts_10_2p5(ind1+1,:) = qpoint;
        refpts_10_2p5(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:)),:) = refpts_10_2p5_tmp(ind1+1:end,:);
        
        refpts_10_2p5_label{ind1+1,:} = Get_10_2p5_label(refpts_10_2p5_label{ind1},refpts_10_2p5_label{ind2});
        refpts_10_2p5_label(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:))) = refpts_10_2p5_label_tmp(ind1+1:end);
        
    end
    
    refpts_10_2p5_tmp =  refpts_10_2p5;
    refpts_10_2p5_label_tmp =  refpts_10_2p5_label;
end


for i = 1:length(Upper_Ring_Left_points)-1
    
    ind1 = find(strcmpi(Upper_Ring_Left_points(i),refpts_10_2p5_label));
    ind2 = find(strcmpi(Upper_Ring_Left_points(i+1),refpts_10_2p5_label));
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    
    %Place new point in sequence in refpts_10_2p5
    if ismember(i,cross_points) %If first point (ind1) is a cross point, then point should be index ind2-1
        refpts_10_2p5(ind2,:) = qpoint;
        refpts_10_2p5(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:)),:) = refpts_10_2p5_tmp(ind2:end,:);
        
        if i == 1
            tmp2 = refpts_10_2p5_label{ind2};
            refpts_10_2p5_label{ind2} = strcat(tmp2(1:end-1),'q');
            refpts_10_2p5_label(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:))) = refpts_10_2p5_label_tmp(ind2:end);            
        else
            refpts_10_2p5_label{ind2} = strcat(refpts_10_2p5_label{ind2},'q');
            refpts_10_2p5_label(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:))) = refpts_10_2p5_label_tmp(ind2:end);
        end
        
    else
        refpts_10_2p5(ind1+1,:) = qpoint;
        refpts_10_2p5(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:)),:) = refpts_10_2p5_tmp(ind1+1:end,:);
        
        refpts_10_2p5_label{ind1+1,:} = Get_10_2p5_label(refpts_10_2p5_label{ind1},refpts_10_2p5_label{ind2});
        refpts_10_2p5_label(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:))) = refpts_10_2p5_label_tmp(ind1+1:end);
        
    end
    
    refpts_10_2p5_tmp =  refpts_10_2p5;
    refpts_10_2p5_label_tmp =  refpts_10_2p5_label;
end

for i = 1:length(Upper_Ring_Right_points)-1
    ind1 = find(strcmpi(Upper_Ring_Right_points(i),refpts_10_2p5_label));
    ind2 = find(strcmpi(Upper_Ring_Right_points(i+1),refpts_10_2p5_label));
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    %Place new point in sequence in refpts_10_2p5
    if ismember(i,cross_points) %If first point (ind1) is a cross point, then point should be index ind2-1
        refpts_10_2p5(ind2,:) = qpoint;
        refpts_10_2p5(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:)),:) = refpts_10_2p5_tmp(ind2:end,:);
        
        if i == 1
            tmp2 = refpts_10_2p5_label{ind2};
            refpts_10_2p5_label{ind2} = strcat(tmp2(1:end-1),'q');
            refpts_10_2p5_label(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:))) = refpts_10_2p5_label_tmp(ind2:end);            
        else
            refpts_10_2p5_label{ind2} = strcat(refpts_10_2p5_label{ind2},'q');
            refpts_10_2p5_label(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:))) = refpts_10_2p5_label_tmp(ind2:end);
        end
        
    else
        refpts_10_2p5(ind1+1,:) = qpoint;
        refpts_10_2p5(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:)),:) = refpts_10_2p5_tmp(ind1+1:end,:);
        
        refpts_10_2p5_label{ind1+1,:} = Get_10_2p5_label(refpts_10_2p5_label{ind1},refpts_10_2p5_label{ind2});
        refpts_10_2p5_label(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:))) = refpts_10_2p5_label_tmp(ind1+1:end);
        
    end
    
    refpts_10_2p5_tmp =  refpts_10_2p5;
    refpts_10_2p5_label_tmp =  refpts_10_2p5_label;
end

%Add points between two rings
start_end_points_lower_left = {'NFp1q' 'FTTT9h' 'T9h' 'TTP9hq' 'OI1hq' 'OIz'};
start_end_points_upper_left = {'Fp1q' 'FTTT7' 'T7' 'TTP7q' 'O1hq' 'Oz'};

start_end_points_lower_right = {'NFp2q' 'FTTT10h' 'T10h' 'TTP10hq' 'OI2hq' 'OIz'};
start_end_points_upper_right = {'Fp2q' 'FTTT8' 'T8' 'TTP8q' 'O2hq' 'Oz'};

for i = 1:length(start_end_points_lower_left)
    ind_lower_left(i) = find(strcmpi(start_end_points_lower_left(i),refpts_10_2p5_label));
    ind_lower_right(i) = find(strcmpi(start_end_points_lower_right(i),refpts_10_2p5_label));  
    ind_upper_left(i) = find(strcmpi(start_end_points_upper_left(i),refpts_10_2p5_label));
    ind_upper_right(i) = find(strcmpi(start_end_points_upper_right(i),refpts_10_2p5_label));
end

ind_lower_left = [ind_lower_left(1):ind_lower_left(2) ind_lower_left(3) ind_lower_left(4):ind_lower_left(5) ind_lower_left(6)];
ind_lower_right = [ind_lower_right(1):ind_lower_right(2) ind_lower_right(3) ind_lower_right(4):ind_lower_right(5) ind_lower_right(6)];
ind_upper_left = [ind_upper_left(1):ind_upper_left(2) ind_upper_left(3) ind_upper_left(4):ind_upper_left(5) ind_upper_left(6)];
ind_upper_right = [ind_upper_right(1):ind_upper_right(2) ind_upper_right(3) ind_upper_right(4):ind_upper_right(5) ind_upper_right(6)];

for i = 1:length(ind_lower_left)
    ind1 = ind_lower_left(i);
    ind2 = ind_upper_left(i);
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    refpts_10_2p5(end+1,:) = qpoint;
    refpts_10_2p5_label{end+1,:} = strcat(refpts_10_2p5_label{ind1},'q');
end

for i = 1:length(ind_lower_right)
    ind1 = ind_lower_right(i);
    ind2 = ind_upper_right(i);
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    refpts_10_2p5(end+1,:) = qpoint;
    refpts_10_2p5_label{end+1,:} = strcat(refpts_10_2p5_label{ind1},'q');
end

%Occipital region
occ_line_left =  {'POO7'    'POO5'    'POO3'    'POO1'    'POOz'};  
occ_line_right = {'POOz'    'POO2'    'POO4'    'POO6'    'POO8'};
frontal_line_left = {'AFp7'  'AFp5'  'AFp3'    'AFp1'    'AFpz'};
frontal_line_right = {'AFpz' 'AFp2'  'AFp4'    'AFp6'    'AFp8'};

%Update tmp
refpts_10_2p5_tmp =  refpts_10_2p5;
refpts_10_2p5_label_tmp =  refpts_10_2p5_label;

for i = 1:length(occ_line_left)-1
    ind1 = find(strcmpi(occ_line_left(i),refpts_10_2p5_label));
    ind2 = find(strcmpi(occ_line_left(i+1),refpts_10_2p5_label));
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    if i == 4
        refpts_10_2p5(ind1+1,:) = qpoint;
        refpts_10_2p5(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:)),:) = refpts_10_2p5_tmp(ind1+1:end,:);
        
        refpts_10_2p5_label{ind1+1} = strcat(refpts_10_2p5_label{ind1},'h');
        refpts_10_2p5_label(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:))) = refpts_10_2p5_label_tmp(ind1+1:end);
    else
        
        refpts_10_2p5(ind2,:) = qpoint;
        refpts_10_2p5(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:)),:) = refpts_10_2p5_tmp(ind2:end,:);
        
        refpts_10_2p5_label{ind2} = strcat(refpts_10_2p5_label{ind1},'h');
        refpts_10_2p5_label(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:))) = refpts_10_2p5_label_tmp(ind2:end);
        
    end
    
    refpts_10_2p5_tmp =  refpts_10_2p5;
    refpts_10_2p5_label_tmp =  refpts_10_2p5_label;
end

for i = 1:length(occ_line_right)-1
    ind1 = find(strcmpi(occ_line_right(i),refpts_10_2p5_label));
    ind2 = find(strcmpi(occ_line_right(i+1),refpts_10_2p5_label));
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    if i == 4
        refpts_10_2p5(ind1+1,:) = qpoint;
        refpts_10_2p5(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:)),:) = refpts_10_2p5_tmp(ind1+1:end,:);
        
        refpts_10_2p5_label{ind1+1} = strcat(refpts_10_2p5_label{ind2},'h');
        refpts_10_2p5_label(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:))) = refpts_10_2p5_label_tmp(ind1+1:end);
    else
        
        refpts_10_2p5(ind2,:) = qpoint;
        refpts_10_2p5(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:)),:) = refpts_10_2p5_tmp(ind2:end,:);
        
        refpts_10_2p5_label{ind2} = strcat(refpts_10_2p5_label{ind2},'h');
        refpts_10_2p5_label(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:))) = refpts_10_2p5_label_tmp(ind2:end);
        
    end
    
    refpts_10_2p5_tmp =  refpts_10_2p5;
    refpts_10_2p5_label_tmp =  refpts_10_2p5_label;
end

for i = 1:length(frontal_line_left)-1
    ind1 = find(strcmpi(frontal_line_left(i),refpts_10_2p5_label));
    ind2 = find(strcmpi(frontal_line_left(i+1),refpts_10_2p5_label));
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    if i == 4
        refpts_10_2p5(ind1+1,:) = qpoint;
        refpts_10_2p5(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:)),:) = refpts_10_2p5_tmp(ind1+1:end,:);
        
        refpts_10_2p5_label{ind1+1} = strcat(refpts_10_2p5_label{ind1},'h');
        refpts_10_2p5_label(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:))) = refpts_10_2p5_label_tmp(ind1+1:end);
    else
        
        refpts_10_2p5(ind2,:) = qpoint;
        refpts_10_2p5(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:)),:) = refpts_10_2p5_tmp(ind2:end,:);
        
        refpts_10_2p5_label{ind2} = strcat(refpts_10_2p5_label{ind1},'h');
        refpts_10_2p5_label(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:))) = refpts_10_2p5_label_tmp(ind2:end);
        
    end
    
    refpts_10_2p5_tmp =  refpts_10_2p5;
    refpts_10_2p5_label_tmp =  refpts_10_2p5_label;
end

for i = 1:length(frontal_line_right)-1
    ind1 = find(strcmpi(frontal_line_right(i),refpts_10_2p5_label));
    ind2 = find(strcmpi(frontal_line_right(i+1),refpts_10_2p5_label));
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    if i == 4
        refpts_10_2p5(ind1+1,:) = qpoint;
        refpts_10_2p5(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:)),:) = refpts_10_2p5_tmp(ind1+1:end,:);
        
        refpts_10_2p5_label{ind1+1} = strcat(refpts_10_2p5_label{ind2},'h');
        refpts_10_2p5_label(ind1+2:ind1+1+length(refpts_10_2p5_tmp(ind1+1:end,:))) = refpts_10_2p5_label_tmp(ind1+1:end);
    else
        
        refpts_10_2p5(ind2,:) = qpoint;
        refpts_10_2p5(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:)),:) = refpts_10_2p5_tmp(ind2:end,:);
        
        refpts_10_2p5_label{ind2} = strcat(refpts_10_2p5_label{ind2},'h');
        refpts_10_2p5_label(ind2+1:ind2+length(refpts_10_2p5_tmp(ind2:end,:))) = refpts_10_2p5_label_tmp(ind2:end);
        
    end
    
    refpts_10_2p5_tmp =  refpts_10_2p5;
    refpts_10_2p5_label_tmp =  refpts_10_2p5_label;
end

%Left and Right temporal edges
Temp_edge_upper_left = {'AF7h'
    'AFAFF7h'
    'AFF7h'
    'AFFF7h'
    'F7h'
    'FFFT7h'
    'FFT7h'
    'FFTFT7h'
    'FT7h'
    'FTFTT7h'
    'FTT7h'
    'FTTT7h'
    'T7h'
    'TTTP7h'
    'TTP7h'
    'TTPTP7h'
    'TP7h'
    'TPTPP7h'
    'TPP7h'
    'TPPP7h'
    'P7h'
    'PPPO7h'
    'PPO7h'
    'PPOPO7h'
    'PO7h'};

Temp_edge_lower_left = {'AF7'
    'AFAFF7'
    'AFF7'
    'AFFF7'
    'F7'
    'FFFT7'
    'FFT7'
    'FFTFT7'
    'FT7'
    'FTFTT7'
    'FTT7'
    'FTTT7'
    'T7'
    'TTP7q'
    'TTP7'
    'TTPTP7'
    'TP7'
    'TPTPP7'
    'TPP7'
    'TPPP7'
    'P7'
    'PPPO7'
    'PPO7'
    'PPOPO7'
    'PO7'};
    
Temp_edge_upper_right = {'AF8h'
    'AFAFF8h'
    'AFF8h'
    'AFFF8h'
    'F8h'
    'FFFT8h'
    'FFT8h'
    'FFTFT8h'
    'FT8h'
    'FTFTT8h'
    'FTT8h'
    'FTTT8h'
    'T8h'
    'TTTP8h'
    'TTP8h'
    'TTPTP8h'
    'TP8h'
    'TPTPP8h'
    'TPP8h'
    'TPPP8h'
    'P8h'
    'PPPO8h'
    'PPO8h'
    'PPOPO8h'
    'PO8h'};

Temp_edge_lower_right = {'AF8'
    'AFAFF8'
    'AFF8'
    'AFFF8'
    'F8'
    'FFFT8'
    'FFT8'
    'FFTFT8'
    'FT8'
    'FTFTT8'
    'FTT8'
    'FTTT8'
    'T8'
    'TTP8q'
    'TTP8'
    'TTPTP8'
    'TP8'
    'TPTPP8'
    'TPP8'
    'TPPP8'
    'P8'
    'PPPO8'
    'PPO8'
    'PPOPO8'
    'PO8'};

%Calculate Left and right temporal edges
for i = 1:length(Temp_edge_lower_left)
    ind1 = find(strcmpi(Temp_edge_lower_left(i),refpts_10_2p5_label));
    ind2 = find(strcmpi(Temp_edge_upper_left(i),refpts_10_2p5_label));
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    refpts_10_2p5(end+1,:) = qpoint;
    refpts_10_2p5_label{end+1} = strcat(refpts_10_2p5_label{ind1},'q');
end

for i = 1:length(Temp_edge_lower_right)
    ind1 = find(strcmpi(Temp_edge_lower_right(i),refpts_10_2p5_label));
    ind2 = find(strcmpi(Temp_edge_upper_right(i),refpts_10_2p5_label));
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    refpts_10_2p5(end+1,:) = qpoint;
    refpts_10_2p5_label{end+1} = strcat(refpts_10_2p5_label{ind1},'q');
end

Frontal_edge_left_lower = { 'AFp7h'    'AFp5'    'AFp5h'    'AFp3'    'AFp3h'    'AFp1'    'AFp1h'};
Frontal_edge_left_upper = { 'AF7h'    'AF5'    'AF5h'    'AF3'    'AF3h'    'AF1'    'AF1h'};
Frontal_edge_right_lower = {'AFp2h'    'AFp2'    'AFp4h'    'AFp4'    'AFp6h'    'AFp6'    'AFp8h'};
Frontal_edge_right_upper = {'AF2h'    'AF2'    'AF4h'    'AF4'    'AF6h'    'AF6'    'AF8h'};

for i = 1:length(Frontal_edge_left_lower)
    ind1 = find(strcmpi(Frontal_edge_left_lower(i),refpts_10_2p5_label));
    ind2 = find(strcmpi(Frontal_edge_left_upper(i),refpts_10_2p5_label));
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    refpts_10_2p5(end+1,:) = qpoint;
    refpts_10_2p5_label{end+1} = strcat(refpts_10_2p5_label{ind1},'q');
end
for i = 1:length(Frontal_edge_right_lower)
    ind1 = find(strcmpi(Frontal_edge_right_lower(i),refpts_10_2p5_label));
    ind2 = find(strcmpi(Frontal_edge_right_upper(i),refpts_10_2p5_label));
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    refpts_10_2p5(end+1,:) = qpoint;
    refpts_10_2p5_label{end+1} = strcat(refpts_10_2p5_label{ind1},'q');
end


Occ_edge_left_lower = {'POO7'    'POO7h'    'POO5'    'POO5h'    'POO3'    'POO3h'    'POO1'    'POO1h'};
Occ_edge_left_upper = { 'PO7'    'PO7h'    'PO5'    'PO5h'    'PO3'    'PO3h'    'PO1'    'PO1h'};
Occ_edge_right_lower = {'POO2h'    'POO2'    'POO4h'    'POO4'    'POO6h'    'POO6'    'POO8h'    'POO8'};
Occ_edge_right_upper = {'PO2h'    'PO2'    'PO4h'    'PO4'    'PO6h'    'PO6'    'PO8h'    'PO8'};

for i = 1:length(Occ_edge_left_lower)
    ind1 = find(strcmpi(Occ_edge_left_lower(i),refpts_10_2p5_label));
    ind2 = find(strcmpi(Occ_edge_left_upper(i),refpts_10_2p5_label));
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    refpts_10_2p5(end+1,:) = qpoint;
    refpts_10_2p5_label{end+1} = strcat(refpts_10_2p5_label{ind1},'q');
end
for i = 1:length(Occ_edge_right_lower)
    ind1 = find(strcmpi(Occ_edge_right_lower(i),refpts_10_2p5_label));
    ind2 = find(strcmpi(Occ_edge_right_upper(i),refpts_10_2p5_label));
    
    tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
    qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
    
    refpts_10_2p5(end+1,:) = qpoint;
    refpts_10_2p5_label{end+1} = strcat(refpts_10_2p5_label{ind1},'q');
end

%Fill in frontal gaps manually
ind1 = find(strcmpi('Fp1hq',refpts_10_2p5_label));
ind2 = find(strcmpi('AFp5',refpts_10_2p5_label));
tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
refpts_10_2p5(end+1,:) = qpoint;
refpts_10_2p5_label{end+1} = 'FAFp5';

ind1 = find(strcmpi('Fp1h',refpts_10_2p5_label));
ind2 = find(strcmpi('AFp3',refpts_10_2p5_label));
tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
refpts_10_2p5(end+1,:) = qpoint;
refpts_10_2p5_label{end+1} = 'FAFp3';

ind1 = find(strcmpi('Fp1q',refpts_10_2p5_label));
ind2 = find(strcmpi('AFp1',refpts_10_2p5_label));
tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
refpts_10_2p5(end+1,:) = qpoint;
refpts_10_2p5_label{end+1} = 'FAFp1';

ind1 = find(strcmpi('Fp2q',refpts_10_2p5_label));
ind2 = find(strcmpi('AFp2',refpts_10_2p5_label));
tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
refpts_10_2p5(end+1,:) = qpoint;
refpts_10_2p5_label{end+1} = 'FAFp2';

ind1 = find(strcmpi('Fp2h',refpts_10_2p5_label));
ind2 = find(strcmpi('AFp4',refpts_10_2p5_label));
tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
refpts_10_2p5(end+1,:) = qpoint;
refpts_10_2p5_label{end+1} = 'FAFp4';

ind1 = find(strcmpi('Fp2hq',refpts_10_2p5_label));
ind2 = find(strcmpi('AFp6',refpts_10_2p5_label));
tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
refpts_10_2p5(end+1,:) = qpoint;
refpts_10_2p5_label{end+1} = 'FAFp6';

%Fill in occipital gaps manually
ind1 = find(strcmpi('O1q',refpts_10_2p5_label));
ind2 = find(strcmpi('POO5',refpts_10_2p5_label));
tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
refpts_10_2p5(end+1,:) = qpoint;
refpts_10_2p5_label{end+1} = 'POOO5';

ind1 = find(strcmpi('O1h',refpts_10_2p5_label));
ind2 = find(strcmpi('POO3',refpts_10_2p5_label));
tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
refpts_10_2p5(end+1,:) = qpoint;
refpts_10_2p5_label{end+1} = 'POOO3';

ind1 = find(strcmpi('O1hq',refpts_10_2p5_label));
ind2 = find(strcmpi('POO1',refpts_10_2p5_label));
tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
refpts_10_2p5(end+1,:) = qpoint;
refpts_10_2p5_label{end+1} = 'POOO1';

ind1 = find(strcmpi('O2hq',refpts_10_2p5_label));
ind2 = find(strcmpi('POO2',refpts_10_2p5_label));
tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
refpts_10_2p5(end+1,:) = qpoint;
refpts_10_2p5_label{end+1} = 'POOO2';

ind1 = find(strcmpi('O2h',refpts_10_2p5_label));
ind2 = find(strcmpi('POO4',refpts_10_2p5_label));
tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
refpts_10_2p5(end+1,:) = qpoint;
refpts_10_2p5_label{end+1} = 'POOO4';

ind1 = find(strcmpi('O2q',refpts_10_2p5_label));
ind2 = find(strcmpi('POO6',refpts_10_2p5_label));
tmp = refpts_10_2p5(ind1,:) + (refpts_10_2p5(ind2,:)-refpts_10_2p5(ind1,:))./2;
qpoint = nearestNode(tmp,mshs.scalpSurfaceMeshUpsampled.node);
refpts_10_2p5(end+1,:) = qpoint;
refpts_10_2p5_label{end+1} = 'POOO6';

% Save results in mshs
mshs.ADSolutionSpace.positions = refpts_10_2p5;
mshs.ADSolutionSpace.labels = refpts_10_2p5_label;

%PLOT ###########################
% figure;
% plotmesh(refpts_10_2p5,'r.','MarkerSize',15);
% for ii = 1:length(refpts_10_2p5)
%     hold on;text(refpts_10_2p5(ii,1),refpts_10_2p5(ii,2),refpts_10_2p5(ii,3)+3,refpts_10_2p5_label(ii));
% end
% 
% figure;
% plotmesh(refpts_10_2p5,'r.','MarkerSize',15);
% for ii = 1:length(refpts_10_2p5);
%     hold on;text(refpts_10_2p5(ii,1),refpts_10_2p5(ii,2),refpts_10_2p5(ii,3)+3,num2str(ii));
% end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nnode,ind,offset] = nearestNode(point,nodes)

% This function outputs the node (1x3) present in the list nodes (Mx3) that is closest
% to the 3D location specified by point (1x3) and the offset in whatever
% dimensions are provided.

% Inputs (point,nodes) : the single point and the list of nodes
% Outputs [nnode, offset]: the index of the specific node in the list
% 'nodes' which is nearest to the input point, and the offset (euclidean
% error).

if size(point,2)~=3
    point = point';
    if size(point,2)~=3
        error('Dimensions of input point not 1x3');
    end
end

if size(nodes,2)~=3
    nodes = nodes';
    if size(nodes,2)~=3
        error('Dimensions of input nodes not Mx3');
    end
end

M = size(nodes,1);
dists = sqrt(sum((nodes - repmat(point,M,1)).^2,2));
[offset,ind] = min(dists);

nnode = nodes(ind,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newlabel = Get_10_2p5_label(lab1,lab2)

%Creates new label for point between 10-5 locations lab1 and lab2.
newlabel = 'Err';

LowCase1 = find(isstrprop(lab1,'lower'),1,'first'); if isempty(LowCase1); LowCase1 = 1e6; end
LowCase2 = find(isstrprop(lab2,'lower'),1,'first'); if isempty(LowCase2); LowCase2 = 1e6; end

DigCase1 = find(isstrprop(lab1,'digit'),1,'first'); if isempty(DigCase1); DigCase1 = 1e6; end
DigCase2 = find(isstrprop(lab2,'digit'),1,'first'); if isempty(DigCase2); DigCase2 = 1e6; end

ind1 = min([LowCase1 DigCase1]);
ind2 = min([LowCase2 DigCase2]);

hstr1 = lab1(1:ind1-1);
hstr2 = lab2(1:ind2-1);
fstr1 = lab1(ind1:end);
fstr2 = lab2(ind2:end);

z1 = strfind(lab1,'z');
z2 = strfind(lab2,'z');

if strcmp(hstr1,hstr2) && ~strcmp(fstr1,fstr2) %If header is the same but footer not, it means the new label is along a row, so just the first label with a 'q' added
    newlabel = strcat(lab1,'q');
elseif ~strcmp(hstr1,hstr2) && strcmp(fstr1,fstr2) %If footer is the same but header not, it means the new label is between rows, so need to combine header and keep foter
    newlabel = strcat(hstr1,hstr2,fstr1);
elseif z2 & z1
   newlabel = strcat(lab1(1:z1-1),lab2(1:z2-1),'z');
end

if strcmp(newlabel,'Err')
    teststop = 1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [newnode,newface] = meshUpsample(node,face,nodeonlyflag)

% This function upsample the scalp surface mesh by dividing each triangle
% into 4 triangles

if nodeonlyflag
    h = waitbar(0,'Upsampling Mesh...');
    for i = 1:length(face)
        
        if mod(i,50) == 0
            waitbar(i/length(face));
        end
        
        nodetmp(1,:) = node(face(i,1),:);
        nodetmp(2,:) = node(face(i,2),:);
        nodetmp(3,:) = node(face(i,3),:);
        node(end+1,:) = mean(nodetmp([1 2],:));
        node(end+1,:) = mean(nodetmp([2 3],:));
        node(end+1,:) = mean(nodetmp([3 1],:));
        node(end+1,:) = mean(nodetmp([1 2 3],:));
        
    end
    close(h)
    newnode = node;
    newface = [];
else
    
    ncount = 1;
    facecount = 1;
    newnode(1,1:3) = [1e6 1e6 1e6];
    newface = [];
    
    h = waitbar(0,'Upsampling Mesh...');
    for i = 1:length(face)
        
        if mod(i,50) == 0
            waitbar(i/length(face));
        end
        
        %Find 6 nodes in current set
        node_set_pos(1,:) = node(face(i,1),:);
        node_set_pos(2,:) = node(face(i,2),:);
        node_set_pos(3,:) = node(face(i,3),:);
        node_set_pos(4,:) = mean(node_set_pos([1 2],:));
        node_set_pos(5,:) = mean(node_set_pos([2 3],:));
        node_set_pos(6,:) = mean(node_set_pos([3 1],:));
        
        ncount_in = 0;
        newnode_in = [];
        for j = 1:6  %Check
            ind = find(ismember(newnode,node_set_pos(j,:),'rows'));
            if ind
                node_set_ind(j) = ind;
            else
                ncount_in = ncount_in + 1;
                node_set_ind(j) = ncount+ncount_in;
                newnode_in(ncount_in,:) = node_set_pos(j,:);
                
            end
        end
        
        newface(facecount,:) =   [node_set_ind(1) node_set_ind(4) node_set_ind(6)];
        newface(facecount+1,:) = [node_set_ind(2) node_set_ind(4) node_set_ind(5)];
        newface(facecount+2,:) = [node_set_ind(3) node_set_ind(5) node_set_ind(6)];
        newface(facecount+3,:) = [node_set_ind(4) node_set_ind(5) node_set_ind(6)];
        
        facecount = facecount+4;
        
        if ncount_in
            newnode(ncount+1:ncount+ncount_in,:) = newnode_in;
        end
        ncount = length(newnode);
    end
    close(h)
    newnode = newnode(2:end,:);
    newface = newface - 1;
end

end