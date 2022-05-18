% a1_SegmentQuantifyMovie_Wulan_v5.m
% Wulan 2017-2019. 
% v4, this is to add a control region inside of the nucleus and use that
% for normalization and compare results with normalization to the whole nucleus. 
% v5, Wulan, April 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%   Before Running This Script      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *** TRANSFORM TO TIF AND ALIGN STACKS: 
% Use FIJI macro to process all FRAP movie files with .czi extension: 
% 1. transform to tiff "batch_convert_czi_to_tif.ijm"
% 2. align stacks using StackReg "batch_align_stacks_StackReg.ijm"  
% Name the tif files with a "_StackReg" at the end. 

% For a single file, plugin -> 'Registration' -> 'StackReg'-'Rigid body'.  
% In this way the movie is mostly aligned. You can manually align it 
% further by set WhatToDo = 2 in this matlab script. 

% *** CACULATE THE BLEACH CENTER POSITION:  
% To extract the information of the bleaching ROI, use the python code that
% Alec has helped me to write up. It reads the czi files using library
% czifiles, and then extract the metadata and extract the metric position of
% the rectangualr ROI (the whole field of view) and the metric position of
% the circular ROI (bleaching spot). Calculate the distance and divide it
% by the pixel size, to calculate the position of the bleaching spot in
% pixel unit. 
% April 2020, open with fiji bioFormats, select "display ome-xml metadata"
% and "display ROIs", to get the bleach spot center value. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   Description      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The purpose of this script is to analyze INDIVIDUAL FRAP movie, mannually 
% set the bleach spot center and background spot center. Then quantify the
% intensity of the circular bleached region and normalize it to the intensity
% of the entire nucleus. 

% Analyze one file in the folder at one time because you will need to
% manually change the positions of the bleach center and background
% center, and manually set the intensity threshold. This is done
% INDIVIDUALLY for every single FRAP movie. Save the information in txt file
% in each experiment folder for future reference. 

% Output is a mat file including all the parameters and intensity data. It is saved
% in the folder of './FRAP_workspaces_v2/' by default.  The file name is SampleName. 
% So it is wise to include the date in the sample name.
% Use the mat file for the next steps. 

% For each folder of tiff, after analysis save the parameters in a text
% file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; 
close all;

%% Directory and shared parameters
% folder containing the tiff files to be analyzed
directory = '/Users/wdeng/OneDrive - 北京大学生物医学前沿创新中心/BigRawData2020/2017-2019FRAPdata/mat_FRAP_2020April/KI14_DE_2Hz_short260/';
input_path = [directory, '3_tiff_StackReg\'];
output_path = [directory, '4_mat_FRAPdata\'];

disp('Listing all files in the directory: '); dir([input_path,'*.tif']);

WhatToDo = 1;  % WhatToDo = 1 to determine the bleach spot location. WhatToDo = 2 to adjust nucleus movement (Do not use if already use ImageJ to align frap video before this script, ImageJ alignment is recommended). 

% Muller, Wach, McNally, BioPhys J, 2008: Can ignore Gaussian shape of photobleach spot when the recovery is slow. 

ExtraRadius = 0; % ExtraRadius = 10; % UNIT IS PIXEL  % don't use extra radius for mitotic cells % add extra pixels out of CircleRadius calculating intensity before bleaching, ZW always keep 0.
ExpFactor = 12; %exp(-endFrame/(ExpFactor*endFrame)) = 0.8187 
                %e.g. 2.2*exp(-480/(12*480))  e.g 3.3*exp(-310/(12*310)) = 3.036
                %Adjust the intensity threshold used across the movie to account for photobleaching. Assume ~18% photobleach over the 300 frames.
MakeMovie = 0;  % usually keep it 0, unless you want to adjust nucleus movement

% Define the FRAP circle:
FirstBleachFrame = 21;   % the first frame after bleaching

% Zuhui modified time series, two stage of recording
% FrameRate = 2;  %frames per second
timing.Exposure = 0.142;
timing.Interval = [0 0 0 10];
timing.Duration = [2.84 152.35 29.86 2*60+30];
timing.Loops = [20 2 210 16];
timing.bleachT = flip(-[timing.Exposure:timing.Exposure:timing.Exposure*timing.Loops(1)],2);
timing.recordT1 = [0:timing.Exposure:timing.Exposure*(timing.Loops(3)-1)];
timing.recordT2 = [timing.Exposure*(timing.Loops(3)-1)+timing.Interval(4):timing.Interval(4):timing.Exposure*timing.Loops(3)+timing.Interval(4)*timing.Loops(4)];
time = [timing.bleachT timing.recordT1 timing.recordT2];

% ZWmodified: bleach ROI is a 1um*1um square under 170 nm/pixel ~6 pixel,
% equal to SmallCircleRadius
CircleRadius = 6;  % UNIT IS PIXEL, in my FRAP practice,pixel size is 100 nm/pixel, the bleaching spot is a circle with 10pixel diameter. 
SmallCircleRadius = 3; % UNIT IS PIXEL, in my FRAP practice,pixel size is 100 nm/pixel. 

% these parameters are usually not used. 
FrameUpdate = 20; %Update every 10th frame; %Define FRAP Cirlce movement % this is used for align the video manually in this matlab code. 
U2OS_background = 0;   %For some of the U2-OS movies, the background was hidden. Need to adjust the background manually (Otherwise FRAP recovery is negative in first time-point). 

% make output folder if it does not already exist
if exist(output_path) == 7
    disp('The given output folder exists. MATLAB workspaces will be save to:'); disp(output_path);
else
    mkdir(output_path); disp('The given output folder did not exist, but was just created. MATLAB workspaces will be save to:'); disp(output_path);
end
%% Invidual tiff video 
SampleName = 'nodrug-008_StackReg';  % without .tif
SaveVal = 1; % Change to 1 if you want to save the variables into a matlab file
IntThres = 800;  
BleachCentre = [298, 86];    % green circle
BackgroundCentre = [23, 81];  % white
ControlCentre = [316, 39];  % yellow



% intensity threshold to segment the nucleus.
% the center of the bleaching spot
% select a circle outside of the nucleus as background
% select a circle inside of the nucleus as an additional control region besides of the entire nucleus. 
%% Movement of the nucleus (optional, usually do not use it) 
    %%%%% The default movement is none, but if you need to change the
    %%%%% movement of the cell, change here %%%%%%%%%%%%%%%%%%%%%%%%%
 
Movement = [0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; ...
            0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; ...
            0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; ...
            0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; ...
            0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0; 0,0]; 

%% Load the image video
% Read the image file
file = strcat(SampleName,'.tif');
[stack, nbImages] = tiffread([input_path, file]);   %nbimages contains the number of images read
                                              %stack is a structure array
% Convert to a matrix
images = zeros(stack(1,1).height, stack(1,1).width, size(stack,2));  % height and width of the FRAP images, size is the number of frames in the movie %
ImHeight = stack(1,1).height;
ImWidth = stack(1,1).width;
fprintf('\nTiff width: %d; height: %d\n', ImWidth, ImHeight);
for i=1:size(images,3)
    images(:,:,i) = stack(1,i).data;
end

% Define time
% time = 1/FrameRate.*(-(FirstBleachFrame-1):(size(images, 3)-FirstBleachFrame));

    
% Make smoothed images for segmentation of the nucleus
smooth_filter = fspecial('gaussian', [6 6], 2);   %Perform gaussian smoothing % fspecial: Create predefined 2-D filter % 'gaussian' Gaussian lowpass filter
smooth_images = zeros(stack(1,1).height, stack(1,1).width, size(stack,2));
for i=1:size(images,3)
    smooth_images(:,:,i) = imfilter(images(:,:,i), smooth_filter);
end

%% WhatToDo = 1: Determine bleach regions
if WhatToDo == 1
    
    % creat 2D coordinates for all the regions of interest
    [rr cc] = meshgrid(1:ImWidth,1:ImHeight);  % Create 2-D grid coordinates for the entire image area
    BW_circle = sqrt((rr-BleachCentre(1,1)).^2+(cc-BleachCentre(1,2)).^2)<=CircleRadius;  % Create 2-D grid coordinates for the bleach spot with certain radius (here it is 10 pixel)
    Background_circle = sqrt((rr-BackgroundCentre(1,1)).^2+(cc-BackgroundCentre(1,2)).^2)<=CircleRadius; % Create 2-D grid coordinates for the background spot with certain radius (here it is 10 pixel)
    Control_circle = sqrt((rr-ControlCentre(1,1)).^2+(cc-ControlCentre(1,2)).^2)<=CircleRadius * 1.5;  % Create 2-D grid coordinates for the nucleus control spot with certain radius (here it is 10 pixel  * 1.5)
    
    %---------------------- Figure 1: define bleach center ------------------------%   
    
    % plot out the images before and after bleaching and overlay with regions of interest, so we can examine the bleaching center 
    figure('position',[10 80 ImWidth*1 ImHeight*3]); %[x y width height]
    
    subplot(3,1,1);
    overlay1 = imoverlay(mat2gray(images(:,:,FirstBleachFrame - 1)), bwperim(BW_circle), [.3 1 .3]);  % green circle    
    overlay2 = imoverlay(overlay1, bwperim(Background_circle), [1 1 1]); % white circle 
    overlay3 = imoverlay(overlay2, bwperim(Control_circle), [1 1 0]); % yellow circle
    imagesc(overlay3);
    title('Before bleaching');
    
    subplot(3,1,2);
    overlay1 = imoverlay(mat2gray(images(:,:,FirstBleachFrame)), bwperim(BW_circle), [.3 1 .3]);  
    overlay2 = imoverlay(overlay1, bwperim(Background_circle), [1 1 1]); 
    overlay3 = imoverlay(overlay2, bwperim(Control_circle), [1 1 0]);
    imagesc(overlay3);
    title('After bleaching'); 
    
    % subtract before image from after image, help to reveal the bleach spot
    subplot(3,1,3);
    z = imsubtract(stack(FirstBleachFrame).data, stack(FirstBleachFrame-1).data);  
    zz = imadjust (z, [], [], 0.1);  % gamma = 0.1
    overlay0 = imoverlay(zz, bwperim(BW_circle), [.3 1 .3]);  
    imagesc(overlay0);  % scaled color 
    title('After - Before');
    impixelinfo;

    %---------------------- Figure 2: nucleus segmentation with threshold ------------------------%

    figure('position',[10 40 ImWidth*1 ImHeight*2]); %[x y width height]
    
    subplot(2,1,1);   
    Nucleus_binary = smooth_images(:,:,FirstBleachFrame) > IntThres*exp(-FirstBleachFrame/(ExpFactor*length(time))); % nucleus segmentation of the bleached frame
    Nucleus_binary_filled = imfill(Nucleus_binary, 'holes');  % fill in any small holes in the middle of the nucleus
    overlay1 = imoverlay(mat2gray(images(:,:,FirstBleachFrame)), bwperim(Nucleus_binary_filled), [1 0 0]);  % red outline of the nucleus
    overlay2 = imoverlay(overlay1, bwperim(Background_circle), [1 1 1]); 
    overlay3 = imoverlay(overlay2, bwperim(Control_circle), [1 1 0]);
    overlay4 = imoverlay(overlay3, bwperim(BW_circle), [.3 1 .3]); 
    imagesc(overlay4);
    title('Nucleus segmentation: bleaching frame');
    
    subplot(2,1,2);
    lastframe = size(images, 3);
    Nucleus_binary = smooth_images(:,:,lastframe) > IntThres*exp(-lastframe/(ExpFactor*length(time)));  % nucleus segmentation of the last frame
    Nucleus_binary_filled = imfill(Nucleus_binary, 'holes');     % fill in any small holes in the middle of the nucleus
    overlay1 = imoverlay(mat2gray(images(:,:,lastframe)), bwperim(Nucleus_binary_filled), [1 0 0]);
    overlay2 = imoverlay(overlay1, bwperim(Background_circle), [1 1 1]); 
    overlay3 = imoverlay(overlay2, bwperim(Control_circle), [1 1 0]);
    overlay4 = imoverlay(overlay3, bwperim(BW_circle), [.3 1 .3]); 
    imagesc(overlay4);
    title('Nucleus segmentation: last frame');
    impixelinfo;
    
end

%% WhatToDo = 2: Adjust Nucleus Movement
if WhatToDo == 2
    FRAP_circle = zeros(1, round(size(images,3)/FrameUpdate));  % default FrameUpdate = 10
    %Figure out how much the bleach circle moved
    for n=1:length(FRAP_circle)
        CumMove = Movement(1:n,:); %The cumulative movement of the circle
        CurrCentre = BleachCentre - [sum(CumMove(:,1)), sum(CumMove(:,2))];
        currBackgroundCentre = BackgroundCentre;
        
        %Update the circle
        [rr cc] = meshgrid(1:ImWidth,1:ImHeight);
        BW_circle = sqrt((rr-CurrCentre(1,1)).^2+(cc-CurrCentre(1,2)).^2)<=CircleRadius;
        % <= less or equal to
        %Find the pixels inside the circle:
        [row,col,val] = find(BW_circle); clear val
        FRAP_vals = zeros(1,length(row));
        for i=1:length(row)
            FRAP_vals(1,i) = images(row(i), col(i), n*FrameUpdate);
        end
        FRAP_circle(1,n) = mean(FRAP_vals);
        %Plot the circle fit
        figure('position',[100 100 500 250]); %[x y width height]
        subplot(1,2,1);
        overlay1 = imoverlay(mat2gray(images(:,:,n*FrameUpdate)), bwperim(BW_circle), [.3 1 .3]);
        imagesc(overlay1);
        title(['FRAP spot in frame ', num2str(n*FrameUpdate)]);
        
        %Threshold to identify the nucleus
        Nucleus_binary = smooth_images(:,:,n*FrameUpdate) > IntThres*exp(-n*FrameUpdate/(ExpFactor*length(time)));
        %Now fill in any small holes in the middle of the nucleus:
        Nucleus_binary_filled = imfill(Nucleus_binary, 'holes');
        
        subplot(1,2,2);
        overlay2 = imoverlay(mat2gray(images(:,:,n*FrameUpdate)), bwperim(Nucleus_binary_filled), [1 0 0]);
        imagesc(overlay2);
        title(['Segmented nucleus in frame ', num2str(n*FrameUpdate)]);
        
        print('-dpng',['./images/', file, '_FRAME_', num2str(n*FrameUpdate), '.png' ]);
        close;
    end
    figure; 
    hold on;
    plot(1:length(FRAP_circle), FRAP_circle);
    
    hold off;
end

%% Calculate and plot the intensity in all segmentations 
% 1. Calculate the intensity in the bleaching spot in the intial frames before bleaching %%%%%
InitSpotIntensity = zeros(1,FirstBleachFrame-1);
for n = 1:FirstBleachFrame-1
    CumMove = Movement(1:max([round(n/FrameUpdate) 1]),:); %The cumulative movement of the circle
    CurrCentre = BleachCentre - [sum(CumMove(:,1)), sum(CumMove(:,2))];
    %Update the circle
    [rr cc] = meshgrid(1:ImWidth,1:ImHeight);
    InitSpot_circle = sqrt((rr-CurrCentre(1,1)).^2+(cc-CurrCentre(1,2)).^2)<=CircleRadius+ExtraRadius;
    %Find pixels
    [row,col,val] = find(InitSpot_circle); clear val
        %[row,col] = find(___) returns the row and column subscripts of each nonzero element 
        %in array X using any of the input arguments in previous syntaxes.
        %[row,col,v] = find(___) also returns vector v, which contains the nonzero elements of X.
    InitSpot_vals = zeros(1,length(row));
    for i=1:length(row)        
        InitSpot_vals(1,i) = images(row(i), col(i), n);
    end
    InitSpotIntensity(1,n) = mean(InitSpot_vals);
end

% 2. Calculate intensities in all segmentations %%%%
FRAP_intensity = zeros(1, size(images, 3));
SmallCircle_intensity = zeros(1, size(images, 3));
Nuc_intensity = zeros(1, size(images, 3)); 
Control_intensity = zeros(1, size(images, 3)); % added v4
Black_intensity = zeros(1, size(images, 3));

for n = 1:size(images, 3)
    CumMove = Movement(1:max([round(n/FrameUpdate) 1]),:); %The cumulative movement of the circle
    CurrCentre = BleachCentre - [sum(CumMove(:,1)), sum(CumMove(:,2))];
    currControlCentre = ControlCentre - [sum(CumMove(:,1)), sum(CumMove(:,2))];
    currBackgroundCentre = BackgroundCentre;
    
    % 2.1 Update the circle
    [rr cc] = meshgrid(1:ImWidth,1:ImHeight);
    BW_circle = sqrt((rr-CurrCentre(1,1)).^2+(cc-CurrCentre(1,2)).^2)<=CircleRadius;
    Small_circle = sqrt((rr-CurrCentre(1,1)).^2+(cc-CurrCentre(1,2)).^2)<=SmallCircleRadius;
    Black_circle = sqrt((rr-BackgroundCentre(1,1)).^2+(cc-BackgroundCentre(1,2)).^2)<=CircleRadius;
    Control_circle = sqrt((rr-ControlCentre(1,1)).^2+(cc-ControlCentre(1,2)).^2)<=CircleRadius *2;     % added v4

    
    % 2.2 Find the pixels inside the FRAP circle and smaller FRAP circle:
    [row,col,val] = find(BW_circle); clear val
    FRAP_vals = zeros(1,length(row));
    for i=1:length(row)        
        FRAP_vals(1,i) = images(row(i), col(i), n);
    end
    FRAP_intensity(1,n) = mean(FRAP_vals);
    
    [row2,col2,val] = find(Small_circle); clear val
    Small_vals = zeros(1,length(row2));
    for i=1:length(row2)        
        Small_vals(1,i) = images(row2(i), col2(i), n);
    end
    SmallCircle_intensity(1,n) = mean(Small_vals);
    
    % 2.3 The nucleus 
    % Segment the nucleus 
    Nucleus_binary = smooth_images(:,:,n) > IntThres*exp(-n/(ExpFactor*length(time)));     %Threshold to identify the nucleus is adjusted over the course of photobleaching
    % Now fill in any small holes in the middle of the nucleus:
    Nucleus_binary_filled = imfill(Nucleus_binary, 'holes');
    % Find the pixels inside of the nucelus
    [nuc_row,nuc_col,val] = find(Nucleus_binary_filled); clear val
    nuc_vals = zeros(1,length(nuc_row)); 
    for i=1:length(nuc_row)        
        nuc_vals(1,i) = images(nuc_row(i), nuc_col(i), n);
    end    
    Nuc_intensity(1,n) = mean(nuc_vals);

    
    % 2.4 Find the pixels inside the Control circle:  (added v4) 
    [row3,col3,val] = find(Control_circle); clear val
    Control_vals = zeros(1,length(row3));
    for i=1:length(row3)        
        Control_vals(1,i) = images(row3(i), col3(i), n);
    end
    Control_intensity(1,n) = mean(Control_vals);
    

    % 2.5 Find the pixels inside the background circles:
    [row4,col4,val] = find(Black_circle); clear val
    Black_vals = zeros(1,length(row4));
    for i=1:length(row4)        
        Black_vals(1,i) = images(row4(i), col4(i), n);
    end

    if U2OS_background > 0
        Black_intensity(1,n) = mean(Black_vals) - U2OS_background; 
    else
        Black_intensity(1,n) = mean(Black_vals);
    end
end


% 3. PLOT all the raw intensity quantifications    %%%% 


%----------------------Figure 3, subplot 1  -----------------------------%

fig = figure('PaperOrientation','landscape','Units','inches','position',[1.7812 2.0625 12.6667 7.0625]); %[x y width height]
subplot(2,3,1);
hold on;
plot(time, FRAP_intensity, 'o', 'Color', [237/255, 28/255, 36/255], 'MarkerSize', 3);
plot(time, SmallCircle_intensity, 'o', 'Color', [0/255, 153/255, 204/255] , 'MarkerSize', 3);
plot(time, Black_intensity, 'ko', 'MarkerSize', 3);
plot(time, Nuc_intensity, 'o', 'Color', [100/255, 151/255, 75/255] , 'MarkerSize', 3);
plot(time, Control_intensity, 'o', 'Color', [100/255, 151/255, 255/255] , 'MarkerSize', 3);

% Below are line plots
%plot(time(FirstBleachFrame:end), FRAP_intensity(FirstBleachFrame:end), '-', 'Color', [237/255, 28/255, 36/255], 'LineWidth', 1);
%plot(time(FirstBleachFrame:end), SmallCircle_intensity(FirstBleachFrame:end), '-', 'Color', [0/255, 153/255, 204/255], 'LineWidth', 1);
%plot(time(FirstBleachFrame:end), Black_intensity(FirstBleachFrame:end), 'k-', 'LineWidth', 1);
%plot(time(FirstBleachFrame:end), Back1_intensity(FirstBleachFrame:end), '-', 'Color', [0/255, 161/255, 75/255], 'LineWidth', 1);
%plot(time(FirstBleachFrame:end), Back2_intensity(FirstBleachFrame:end), '-', 'Color', [0/255, 161/255, 75/255], 'LineWidth', 1);
%plot(time(FirstBleachFrame:end), Back3_intensity(FirstBleachFrame:end), '-', 'Color', [0/255, 161/255, 75/255], 'LineWidth', 1);

legend(sprintf('FRAP: r=%d',CircleRadius), sprintf('FRAP: r=%d',SmallCircleRadius), 'Black bg', 'whole Nuc', 'control', 'Location', 'NorthEast');
axis([min(time) max(time) 0 1.55*max(FRAP_intensity)]);
title('Raw FRAP quantification', 'FontSize',10, 'FontName', 'Helvetica');
ylabel('hFOXA2-Halo TMR fluorescence (AU)', 'FontSize',10, 'FontName', 'Helvetica');
xlabel('time (seconds)', 'FontSize',10, 'FontName', 'Helvetica');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalization Analysis and plot the curves

%%%%%%%%%%% PERFORM ANALYSIS AND PHOTOBLEACHING CORRECTION %%%%%%%%%%%%%%%%

%Rationale for the analysis: first of all, subtract the black background.
%Second, use the total nuclear fluorescence to correct for
%photobleaching and reduction in cell fluorescence overall. Furthermore, the cell 
%moves a bit so it is important to take into account that cell total and
%local fluorescence can change a bit. In some of the movies, cell moves
%quite a bit so this is important to do. 

%So do the following:
%1. subtract black background from all. 
%2. Divide the FRAP circle by the average nuclear background throughout. 
%3. Normalize the ratio based on the initial fluorescence before bleaching.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1. SUBTRACT BLACK BACKGROUND
bsFRAP_intensity = FRAP_intensity - Black_intensity;
bsSmallCircle_intensity = SmallCircle_intensity - Black_intensity;
bsNuc_intensity = Nuc_intensity - Black_intensity;
bsControl_intensity = Control_intensity - Black_intensity;
InitSpotIntensity = InitSpotIntensity - Black_intensity(1:FirstBleachFrame-1);

%2. ADJUST FOR PHOTOBLEACHING TO THE WHOLE NUCLEUS 
adjust_nuc_InitSpotIntensity = InitSpotIntensity./(bsNuc_intensity(1:FirstBleachFrame-1));
% Element-wise right division
% X = A./B performs right-array division by dividing each element of A by the corresponding element of B.
% yy = smooth(y,span) sets the span of the moving average to span. span must be odd.
adjust_FRAP = bsFRAP_intensity ./smooth(bsNuc_intensity, 3)';
adjust_smallFRAP = bsSmallCircle_intensity ./smooth(bsNuc_intensity, 3)';

%---------- Figure 3, subplot 2, normalized to the nucleus -------------%
subplot(2,3,2);
hold on;
plot(time, adjust_FRAP, 'o', 'Color', [237/255, 28/255, 36/255], 'MarkerSize', 3);
plot(time, adjust_smallFRAP, 'o', 'Color', [0/255, 153/255, 204/255] , 'MarkerSize', 3);

plot(time(FirstBleachFrame:end), smooth(adjust_FRAP(FirstBleachFrame:end),7)', '-', 'Color', [237/255, 28/255, 36/255], 'LineWidth', 2);
plot(time(FirstBleachFrame:end), smooth(adjust_smallFRAP(FirstBleachFrame:end),7)', '-', 'Color', [0/255, 153/255, 204/255], 'LineWidth', 2);
legend(sprintf('FRAP: r=%d',CircleRadius), sprintf('FRAP: r=%d',SmallCircleRadius),  'Location', 'SouthEast');
axis([min(time) max(time) 0 1.05*max(adjust_smallFRAP)]);
title('Nuclear normalized FRAP', 'FontSize',12, 'FontName', 'Helvetica');
ylabel('Relative hFOXA2-Halo TMR fluorescence', 'FontSize',12, 'FontName', 'Helvetica');
xlabel('Time (seconds)', 'FontSize',12, 'FontName', 'Helvetica');
hold off;

% 3. NORMALIZE (photobleach corrected by the whole nucleus) BASED ON INITIAL FLUORESCENCE 
Init_norm_FRAP = adjust_FRAP./mean(adjust_nuc_InitSpotIntensity);
Init_norm_smallFRAP = adjust_smallFRAP./mean(adjust_nuc_InitSpotIntensity);
%---------- Figure 3, subplot 3, normalized again to the initial intensity before bleaching -------------%
subplot(2,3,3);
hold on;
plot(time, Init_norm_FRAP, 'o', 'Color', [237/255, 28/255, 36/255], 'MarkerSize', 3);
plot(time, Init_norm_smallFRAP, 'o', 'Color', [0/255, 153/255, 204/255] , 'MarkerSize', 3);
plot([min(time) max(time)], [1 1], 'k--', 'LineWidth', 1);
plot(time(FirstBleachFrame:end), smooth(Init_norm_FRAP(FirstBleachFrame:end),7)', '-', 'Color', [237/255, 28/255, 36/255], 'LineWidth', 2);
plot(time(FirstBleachFrame:end), smooth(Init_norm_smallFRAP(FirstBleachFrame:end),7)', '-', 'Color', [0/255, 153/255, 204/255], 'LineWidth', 2);
legend(sprintf('FRAP: r=%d',CircleRadius), sprintf('FRAP: r=%d',SmallCircleRadius),  'Location', 'SouthEast');
axis([min(time) max(time) 0 1.15]);
title('Initial spot normalized FRAP - nucleus', 'FontSize',10, 'FontName', 'Helvetica');
ylabel('relative hFOXA2-Halo TMR fluorescence', 'FontSize',10, 'FontName', 'Helvetica');
xlabel('time (seconds)', 'FontSize',10, 'FontName', 'Helvetica');
hold off;

%4. ADJUST FOR PHOTOBLEACHING TO THE CONTROL REGION 
adjust_control_InitSpotIntensity = InitSpotIntensity./(bsControl_intensity(1:FirstBleachFrame-1));
%./ Element�wise right division
% X = A./B performs right-array division by dividing each element of A by the corresponding element of B.
% yy = smooth(y,span) sets the span of the moving average to span. span must be odd.
adjust_control_FRAP = bsFRAP_intensity ./smooth(bsControl_intensity, 3)';
adjust_control_smallFRAP = bsSmallCircle_intensity ./smooth(bsControl_intensity, 3)';

%---------- Figure 3, subplot 2&3 , normalized to the nucleus -------------%
subplot(2,3,5);
hold on;
plot(time, adjust_control_FRAP, 'o', 'Color', [237/255, 28/255, 36/255], 'MarkerSize', 3);
plot(time, adjust_control_smallFRAP, 'o', 'Color', [0/255, 153/255, 204/255] , 'MarkerSize', 3);

plot(time(FirstBleachFrame:end), smooth(adjust_control_FRAP(FirstBleachFrame:end),7)', '-', 'Color', [237/255, 28/255, 36/255], 'LineWidth', 2);
plot(time(FirstBleachFrame:end), smooth(adjust_control_smallFRAP(FirstBleachFrame:end),7)', '-', 'Color', [0/255, 153/255, 204/255], 'LineWidth', 2);
legend(sprintf('FRAP: r=%d',CircleRadius), sprintf('FRAP: r=%d',SmallCircleRadius),  'Location', 'SouthEast');
axis([min(time) max(time) 0 1.05*max(adjust_control_smallFRAP)]);
title('Control region normalized FRAP', 'FontSize',12, 'FontName', 'Helvetica');
ylabel('Relative hFOXA2-Halo TMR fluorescence', 'FontSize',12, 'FontName', 'Helvetica');
xlabel('Time (seconds)', 'FontSize',12, 'FontName', 'Helvetica');
hold off;

% 5. NORMALIZE (photobleach corrected by the CONTROL REGION) BASED ON INITIAL FLUORESCENCE 
Init_norm_control_FRAP = adjust_control_FRAP./mean(adjust_control_InitSpotIntensity);
Init_norm_control_smallFRAP = adjust_control_smallFRAP./mean(adjust_control_InitSpotIntensity);

subplot(2,3,6);
hold on;
plot(time, Init_norm_control_FRAP, 'o', 'Color', [237/255, 28/255, 36/255], 'MarkerSize', 3);
plot(time, Init_norm_control_smallFRAP, 'o', 'Color', [0/255, 153/255, 204/255] , 'MarkerSize', 3);
plot([min(time) max(time)], [1 1], 'k--', 'LineWidth', 1);
plot(time(FirstBleachFrame:end), smooth(Init_norm_control_FRAP(FirstBleachFrame:end),7)', '-', 'Color', [237/255, 28/255, 36/255], 'LineWidth', 2);
plot(time(FirstBleachFrame:end), smooth(Init_norm_control_smallFRAP(FirstBleachFrame:end),7)', '-', 'Color', [0/255, 153/255, 204/255], 'LineWidth', 2);
legend(sprintf('FRAP: r=%d',CircleRadius), sprintf('FRAP: r=%d',SmallCircleRadius),  'Location', 'SouthEast');
axis([min(time) max(time) 0 1.15]);
title('Initial spot normalized FRAP - control', 'FontSize',10, 'FontName', 'Helvetica');
ylabel('relative hFOXA2-Halo TMR fluorescence', 'FontSize',10, 'FontName', 'Helvetica');
xlabel('time (seconds)', 'FontSize',10, 'FontName', 'Helvetica');
hold off;

%suptitle(SampleName); % add title above all subplots.  % for matlab after 2018b, use sgtitle


%% Make a movie
if MakeMovie == 1
    %Adjust CumMove:
    FullCumMove = zeros(length(time), 2);
    for i=1:size(CumMove,1)
        CurrIndex = FrameUpdate + (i-1)*FrameUpdate;
        FullCumMove(CurrIndex,:) = CumMove(i,:);
    end
    prebleach = mean(adjust_smallFRAP(1:FirstBleachFrame-1));
    adjust_smallFRAP2 = [adjust_smallFRAP(1:FirstBleachFrame-1)./prebleach adjust_smallFRAP(FirstBleachFrame:end)];
           
    %Figure out how much the bleach circle moved
    for n=1:size(FullCumMove,1)
        %Update the FRAP circle centre
        CurrCentre = BleachCentre - [sum(FullCumMove(1:n,1)), sum(FullCumMove(1:n,2))];

        %Update the circle
        [rr cc] = meshgrid(1:ImWidth,1:ImHeight);
        BW_circle = sqrt((rr-CurrCentre(1,1)).^2+(cc-CurrCentre(1,2)).^2)<=CircleRadius;
        
        %Plot the circle fit
        figure('position',[100 100 650 250]); %[x y width height]
        subplot(1,2,1);
        overlay1 = imoverlay(mat2gray(images(:,:,n)), bwperim(BW_circle), [.3 1 .3]);
        imagesc(overlay1);
        title(['hFOXA2-Halo FRAP image at ', num2str(n-FirstBleachFrame), ' seconds']);
        set(gca,'XTick',[]); set(gca,'YTick',[]);
        
        subplot(1,2,2);
        hold on;
        plot(time(1:n), Init_norm_smallFRAP(1:n), 'o', 'Color', [237/255, 28/255, 36/255] , 'MarkerSize', 6);
        axis([min(time) max(time) 0 1.15]);
        title('Quanfied FRAP recovery', 'FontSize',10, 'FontName', 'Helvetica');
        ylabel('relative TMR fluorescence', 'FontSize',10, 'FontName', 'Helvetica');
        xlabel('time (seconds)', 'FontSize',10, 'FontName', 'Helvetica');
        hold off;
        print('-dpng', '-r300',['./images/', file, '_FRAP_Movie_', num2str(n), '.png' ]);
        close;
    end
end


%% Save results into a matlab file

if SaveVal == true
    print(fig,[output_path, SampleName '_FRAPcurve'],'-dpdf','-r0');
    save([output_path, SampleName, '.mat'], 'input_path', 'output_path', 'ExtraRadius', 'ExpFactor', 'FirstBleachFrame', 'timing', 'CircleRadius', 'SmallCircleRadius',...
        'SampleName','IntThres', 'BleachCentre', 'BackgroundCentre', 'ControlCentre', ...
        'time','FRAP_intensity', 'SmallCircle_intensity', 'adjust_FRAP', 'adjust_smallFRAP',...
        'bsNuc_intensity', 'Init_norm_FRAP', 'Init_norm_smallFRAP', 'Init_norm_control_FRAP', 'Init_norm_control_smallFRAP');   
end


