% a2_AnalyzeQuantifiedFRAPdata_v3_Wulan_v2.m
% Based on Sejr Hansen, Jan 2016  
% Modified by Wulan Deng, Nov. 2017 - 2020
% Latest wulan's version, April 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   Description      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The purpose of this script is to COMBINE data from multiple analyzed and
% normalized FRAP curve OF THE SAME CONDITION AND GENOTYPE,
% and to do either one reaction-dominant fitting, or
% two reaction-dominant fitting, using either linear or logarithmically 
% spaced time axis. This script will generate the figure and fitting 
% results on the figures. These figures will be saved into PDF files in the
% folder of SavePDFs_FRAP. 

%%%%%%%%%%%%%% Modified by Wulan Deng, Nov. 2017 %%%%%%%%%%%%%%%%%%%%%%%%
%--1. Added SavePDF function 
%%%%%%%%%%%%%% Modified by Wulan Deng, May. 2018%%%%%%%%%%%%%%%%%%%%%%%%%
%--1. Added a parameter of the end frame of the video to be analyzed
%--2. Added a 3-exp fit
%--3. Simplified the output data on the figure output 
%%%%%%%%%%%%%% Modified by Wulan Deng, May. 2019 %%%%%%%%%%%%%%%%%%%%%%%%
%--1. Read mat files from a folder with two specific key words and 
%--2. Generate the ToKeep list automatically. Can modify it manually.  
%--3. Enhance the figure appearance.

%   v3, logarithmically spaced timepoints for model fitting

% April 2020, add a module to plot a lot of subplots, make new figures when
% the subplots are full. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all; 

%% Path, Keywords, Parameters (CHANGE)
% input folder containning FRAP mat of many cells
input_path  = '/Users/wdeng/OneDrive - 北京大学生物医学前沿创新中心/BigRawData2020/2017-2019FRAPdata/mat_FRAP_2020April/KI14_DE_2Hz_short260/4_mat_FRAPdata/';
% output folder
FolderName = 'MergedFRAP_KI14_DE_2Hz_260'; 
output_path = ['5_fitdata/', FolderName, '/'];
% make output folder if it does not already exist
if exist(output_path) == 7
    disp('The given output folder exists. MATLAB workspaces will be save to:'); 
    disp(output_path);
else
    mkdir(output_path); disp('The given output folder did not exist, but was just created. MATLAB workspaces will be save to:'); disp(output_path);
end

% set parameters 
keyword_1 = 'KI';  
keyword_2 = 'DE';
ExpName = FolderName;   % {' '} is to insert space

%ExpName = '20190419 KI14 DE';  % saved PDF name 
select = 0;  % default 0; select = 0 means using FRAP from all cells, select = 1 means using selected FRAP data
SavePDF = 1;  % 1 means you want to save resutls in pdf. 0 means don't save. 
SingleCellFit = 1;  % 1 mean you will do fitting on the single cell plot, 0 means not. 
LastFitFrame = 260;  % Set last fit frame if you want to customerize it, otherwise set it to 0, and it will use all the collected data for fitting. 
                     %  I set it to 260 for all the frap analysis on Apri. 2020

% other parameters (that usually don't change) 
n_hist = 250;   %   n_hist = number of time points for fitting, used in generating logspacing of the time axis                   
DoModelFitting = 1;  % 1 means do model fitting. or 0 means don't do fitting
FirstFitFrame = 21;  % this is the first after-bleaching frame
LongMovie = 0;   % so far I haven't worked on make movie, so leave it as 0. 



% read all the files in the directory and with the keywords in the filename
files = dir([input_path, '*.mat']); i = 1; disp('All Samples = ');% mat_files = dir([input_path,'*_Tracked.mat']);
for iter = 1: length(files)
    if strfind(files(iter).name, keyword_1) & strfind(files(iter).name, keyword_2)  % . && doesn't work here, only work on scalars; & works on array.  strfind returns the location of the keyword 
        Workspaces (i) = cellstr(files(iter).name); % cellstr transfer a string to cell.  Alternative method:   Workspaces (i) = {files(iter).name};
        disp(Workspaces(i));
        i = i + 1;
    end
end

if select == 1 % select = 1 means that you choose which files to process, put the value to 0 if you don't want to include certain file. 
    
    % change this as needed, exclude cells that move too much. %%%%%% CHANGE %%%%
    ToKeep = [0,0,3,4,5,6,7,8,0,10,11,12,13];
    
else
    n = length(Workspaces);
    ToKeep = linspace(1, n, n)   % ToKeep(i) >= 1 means selecting this file for analysis
end



%% Individual FRAP plot and fitting 

iter = 1;  % this is to iterate FRAP data file
fig_iter = 0; % this is to iterate for ploting new figure.  added April. 2019 by Wulan


for i=1:length(Workspaces)
    if ToKeep(i) >= 1
    
    % ----- Step 1: Load the data of one cell FRAP ------
    % ----- Also compile FRAP data from all cells  ------
        %load([input_path, Workspaces{i}]);
        load([input_path, Workspaces{i}],'-regexp', ['^(?!', 'input_path', '$|', 'output_path',')\w']);  % load all except input and output path
        if LastFitFrame == 0  % if this value is not set, then use all the collected data by set it to the actual last frame
            LastFitFrame = length(Init_norm_smallFRAP); 
        end
      
        if iter == 1
            FRAP_data = Init_norm_smallFRAP(1:LastFitFrame);  % this is to load the first data files 
                                                              % Init_norm_smallFRAP is the FRAP data of the smaller dot r=6; so all the fitting is using data from the small circle
        else
            FRAP_data = vertcat(FRAP_data, Init_norm_smallFRAP(1:LastFitFrame));  % this is to add on the second data file and beyond % vertcat: concatenate arrays vertically
        end
        iter = iter + 1;

    % ----- Step 2: LOG-model fitting with INDIVIDUAL cells -----
    % ----- Added by Wulan, May 2018  

        if DoModelFitting == 1 && SingleCellFit == 1 
            %Variable for fitting:
            xTime1 = time(FirstFitFrame:LastFitFrame);  % time is in the mat file
            xTime2 = time(FirstBleachFrame:LastFitFrame);
            %Model fitting using logarithmically spaced vector (ie.timepoints)
            log_time = logspace(0, log10(time(LastFitFrame)), n_hist);  % logspace: generate logarithmically spaced vector (ie.timepoints) y = logspace(a,b,n) generates n points between decades 10^a and 10^b.

            yFRAP1single = Init_norm_smallFRAP(FirstFitFrame:LastFitFrame);  %%%
            yFRAP2single = Init_norm_smallFRAP(FirstBleachFrame:LastFitFrame);   %%%FirstBleachFrame is in the mat file
            %Now find the closest y_value
            log_yFRAPsingle = zeros(1,length(log_time));
            for k = 1:length(log_yFRAPsingle)
                temp_diff = abs(log_time(1,k)-xTime2); % abs: absolute value
                [val, idx] = min(temp_diff); %min: Smallest elements in array
                log_yFRAPsingle(1,k) = yFRAP2single(1,idx);
            end
            %Now repeat fitting using log-sampled data
            f2_log = fittype('1 - A*exp(-ka*x) - B*exp(-kb*x)');
            [TwoExp_fit_log, TwoExp_param_log] = fit(log_time', log_yFRAPsingle', f2_log, 'Lower', [0 0.1 0 0], 'Upper', [1 10 1 10], 'StartPoint', [0.5 1.1 0.5 0.025]); 
            TwoExp_CI_log = confint(TwoExp_fit_log);
            xFit2_log_single = 0:0.25:max(xTime2);
            yFit2_log_single = 1 - TwoExp_fit_log.A .* exp(-TwoExp_fit_log.ka .* xFit2_log_single) - TwoExp_fit_log.B .* exp(-TwoExp_fit_log.kb .* xFit2_log_single);

            Single_text(1) = {'2-Exp fit - log-sampled : 1-A*exp(-ka*t)-B*exp(-kb*t)'};
            Single_text(2) = {['A = ', num2str(TwoExp_fit_log.A), '     ( ', num2str(TwoExp_CI_log(1,1)), ', ', num2str(TwoExp_CI_log(2,1)), ' )']};
            Single_text(3) = {['1/ka = ', num2str(1/TwoExp_fit_log.ka), 's','     ( ', num2str(1/TwoExp_CI_log(1,3)), ', ', num2str(1/TwoExp_CI_log(2,3)), ' )']};
            Single_text(4) = {['B = ', num2str(TwoExp_fit_log.B), '     ( ', num2str(TwoExp_CI_log(1,2)), ', ', num2str(TwoExp_CI_log(2,2)), ' )']};
            Single_text(5) = {['1/kb = ', num2str(1/TwoExp_fit_log.kb), 's','     ( ', num2str(1/TwoExp_CI_log(1,4)), ', ', num2str(1/TwoExp_CI_log(2,4)), ' )']};

        end


        
    end
        % ----- Step 3: PLOT ----- 
        
        % decide on figure and subplot
             
        if mod(i, 16) == 1  % plot a new figure for every 16 subplots 
            fig_iter = fig_iter +1;
            figure('position',[100 100 1200 1200]); %[x y width height] %[left bottom width height]
        end

        subplot(4,4, (i-(fig_iter-1)*16)); % this is maximum 16 subplot, i is the index of ploted data file. 
                        % subplot(m,n,p) divides the current figure into an m-by-n grid and creates axes in the position specified by p.
        hold on;
        
        
        % plot the FRAP data
        plot(time, Init_norm_smallFRAP, 'ro', 'MarkerSize', 3);  % 'ko' means black circle. 'ro' means red circle
        c = Workspaces{i};
        title([keyword_1, keyword_2, ' ', c(strfind(c, 'cell'):end-4)], 'FontSize',10, 'FontName', 'Helvetica', 'Interpreter', 'none');
        ylabel('fluorescence (AU)', 'FontSize',10, 'FontName', 'Helvetica');
        xlabel('time (seconds)', 'FontSize',10, 'FontName', 'Helvetica');
        axis([min(time) max(time) 0 1.2]);
        
        % plot the fited line
        if DoModelFitting == 1 && SingleCellFit == 1 
            plot(xFit2_log_single, yFit2_log_single, 'k-.', 'LineWidth', 2);  % plot the log fitted line
            plot([min(time), max(time)], [1, 1], 'k--', 'LineWidth', 1);  % plot the top line of 1
            text(10,0.3,Single_text,'HorizontalAlignment','Left',  'FontSize',10, 'FontName', 'Helvetica');
        end
        hold off;
        
        
        % save figure as pdf when  (1) it is full with subplots  or (2) this is the last fig of mulitple figs. 
        if SavePDF == 1 && ( mod(i, 16) == 0 || i == length(Workspaces)) % mod: remainder after divisions.
            h=gcf;
            set(h,'Units','Inches');
            pos = get(h,'Position');
                                     %pos(1) and pos(2) are the coordinate of the legend windows
                                     %pos(3) and pos(4) are width and height of this windows
            set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
            print(h,[output_path, strcat(ExpName, '_figure1_', num2str(fig_iter)), '.pdf'], '-dpdf', '-r0');  % modify the figure name accordingly 
        end

        
end



    
if LongMovie == 1
    time = time ./60;
end 


%% Averaged FRAP plot

meanFRAP = mean(FRAP_data, 1);
stdFRAP = std(FRAP_data);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Plot the FRAP mean and std without fitting  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[100 100 900 600]); %[x y width height]
hold on;

%
time = time (1: LastFitFrame);  % added to adjust axis of time according to the manual setting of LastFitFrame.  May 2019, by Wulan

errorbar(time, meanFRAP, stdFRAP, 'o', 'Color', [237/255, 28/255, 36/255], 'MarkerSize', 2);
plot(time(FirstBleachFrame:end), smooth(meanFRAP(FirstBleachFrame:end),11), 'k-', 'LineWidth', 3);
%title(['FRAP for ', num2str(length(ToKeep)), ' cells'], 'FontSize',10, 'FontName', 'Helvetica');
title(['FRAP for ',ExpName,', ', num2str(length(ToKeep)), ' cells'], 'FontSize',14, 'FontName', 'Helvetica'); % modified by Wulan Nov. 2017
ylabel('fluorescence (AU)', 'FontSize',10, 'FontName', 'Helvetica');
xlabel('time (seconds)', 'FontSize',10, 'FontName', 'Helvetica');
axis([min(time) max(time) 0 1.2]);
legend('mean with error', 'smoothed mean', 'Location', 'SouthEast');
legend boxoff;
hold off;


if SavePDF == 1
    h=gcf;
    set(h,'Units','Inches');
    pos = get(h,'Position');
    %pos(1) and pos(2) are the coordinate of the legend windows
    %pos(3) and pos(4) are width and height of this windows
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,[output_path, strcat(ExpName, '_figure2'), '.pdf'], '-dpdf', '-r0');
end

%% FRAP 1-Exp,2-Exp, 3-Exp fit & plot figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Fitting and Plot  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if DoModelFitting == 1
    %Do a simply two-exponential fitting:
    %FRAP(t) = 1 - A*exp(-ka*t) - B*exp(-kb*t)
    
    % ----- Step 1: Calculate variable for fitting
    
    xTime1 = time(FirstFitFrame:LastFitFrame);
    yFRAP1 = meanFRAP(FirstFitFrame:LastFitFrame);
    xTime2 = time(FirstBleachFrame:LastFitFrame);
    yFRAP2 = meanFRAP(FirstBleachFrame:LastFitFrame);

    % Model fitting using logarithmically spaced vector (ie.timepoints)
    log_time = logspace(0, log10(time(LastFitFrame)), n_hist);  % logspace: generate logarithmically spaced vector (ie.timepoints) y = logspace(a,b,n) generates n points between decades 10^a and 10^b.
    %Now find the closest y_value
    log_yFRAP = zeros(1,length(log_time));
    for k = 1:length(log_yFRAP)
        temp_diff = abs(log_time(1,k)-xTime2); % abs: absolute value
        [val, idx] = min(temp_diff); %min: Smallest elements in array
        log_yFRAP(1,k) = yFRAP2(1,idx);
    end
    
    % ----- Step 2: Four types of fitting 
    % ----- OneExp, TwoExp, ThreeExp, TwoExp_log  
    % ----- Each uses function, f1, f2, f3, f2_log  respectively
 
    f1 = fittype('1 - A*exp(-ka*x) ');
    [OneExp_fit, OneExp_param] = fit(xTime1', yFRAP1', f1, 'Lower', [0.1 0.0001], 'Upper', [0.9 0.5], 'StartPoint', [0.35 0.025]); 
    OneExp_CI = confint(OneExp_fit);
    
    f2 = fittype('1 - A*exp(-ka*x) - B*exp(-kb*x)');
    [TwoExp_fit, TwoExp_param] = fit(xTime2', yFRAP2', f2, 'Lower', [0.05 0 0.05 0], 'Upper', [0.8 4 2 0.05], 'StartPoint', [0.5 1.1 0.5 0.0005]);  % 'Lower', [0.05 0.05 0.05 0], 'Upper', [0.8 4 2 0.05], 'StartPoint', [0.5 1.1 0.5 0.0005]); 
    TwoExp_CI = confint(TwoExp_fit);
    
    %%3-exp fit %%%added by Wulan May 2018, copied from PLOT_MultipleFRAPcurves.m %%
    f3 = fittype('1 - A*exp(-ka*x) - B*exp(-kb*x) - C*exp(-kc*x)');
    [ThreeExp_fit, ThreeExp_param] = fit(xTime2', yFRAP2', f3, 'Lower', [0.1 0.1 0.1 0.1 0.001 0.0001], 'Upper', [0.8 0.8 0.8 10 0.8 0.05], 'StartPoint', [0.3 0.2 0.2 0.25 0.01 0.0009]); 
    ThreeExp_CI = confint(ThreeExp_fit);  
    % CONFINT  Confidence intervals for the coefficients of a fit result object; coeffnames(f3) can show the content, as a cell array
    % B = A' computes the complex conjugate transpose of A.

    f2_log = fittype('1 - A*exp(-ka*x) - B*exp(-kb*x)');
    [TwoExp_fit_log, TwoExp_param_log] = fit(log_time', log_yFRAP', f2_log, 'Lower', [0 0 0 0], 'Upper', [0.9 1 10 100], 'StartPoint', [0.5 1.1 0.5 0.025]);  % original 'Lower', [0 0.1 0 0], 'Upper', [1 10 1 10], 'StartPoint', [0.5 1.1 0.5 0.025]); 
    TwoExp_CI_log = confint(TwoExp_fit_log);
    % coeffnames(f2_log) % it is in the order of A, B, ka, kb % Get the coefficient names and order using the coeffnames function.
    
    % ----- Each generates a set of x, y values for plotting     
    
    xFit1 = min(xTime1):0.25:max(xTime1);
    yFit1 = 1 -  OneExp_fit.A.* exp(-OneExp_fit.ka .* xFit1);
    
    xFit2 = 0:0.25:max(xTime2);
    yFit2 = 1 - TwoExp_fit.A .* exp(-TwoExp_fit.ka .* xFit2) - TwoExp_fit.B .* exp(-TwoExp_fit.kb .* xFit2);
    
    xFit3 = min(xTime2)-1:0.25:max(xTime2);
    yFit3 = 1 - ThreeExp_fit.A .* exp(-ThreeExp_fit.ka .* xFit3) - ThreeExp_fit.B .* exp(-ThreeExp_fit.kb .* xFit3) - ThreeExp_fit.C .* exp(-ThreeExp_fit.kc .* xFit3);
    
    xFit2_log = 0:0.25:max(xTime2);
    yFit2_log = 1 - TwoExp_fit_log.A .* exp(-TwoExp_fit_log.ka .* xFit2_log) - TwoExp_fit_log.B .* exp(-TwoExp_fit_log.kb .* xFit2_log);
        
    % ----- Each generates corresponding text on the figure   
    
    Fit1_text(1) = {'1-Exp fit: 1-A*exp(-ka*t) '};
    Fit1_text(2) = {['A = ', num2str(OneExp_fit.A), '     ( ', num2str(OneExp_CI(1,1)), ', ', num2str(OneExp_CI(2,1)), ' )']};
    Fit1_text(3) = {['1/ka = ', num2str(1/OneExp_fit.ka), 's','     ( ', num2str(1/OneExp_CI(1,2)), ', ', num2str(1/OneExp_CI(2,2)), ' )']};

    Fit2_text(1) = {'2-Exp fit: 1-A*exp(-ka*t)-B*exp(-kb*t)'};
    Fit2_text(2) = {['A = ', num2str(TwoExp_fit.A), '     ( ', num2str(TwoExp_CI(1,1)), ', ', num2str(TwoExp_CI(2,1)), ' )']};
    Fit2_text(3) = {['1/ka = ', num2str(1/TwoExp_fit.ka), 's','     ( ', num2str(1/TwoExp_CI(1,3)), ', ', num2str(1/TwoExp_CI(2,3)), ' )']};
    Fit2_text(4) = {['B = ', num2str(TwoExp_fit.B), '     ( ', num2str(TwoExp_CI(1,2)), ', ', num2str(TwoExp_CI(2,2)), ' )']};
    Fit2_text(5) = {['1/kb = ', num2str(1/TwoExp_fit.kb), 's','     ( ', num2str(1/TwoExp_CI(1,4)), ', ', num2str(1/TwoExp_CI(2,4)), ' )']};
     
    Fit3_text(1) = {'3-Exp fit: 1 - A*exp(-ka*x) - B*exp(-kb*x) - C*exp(-kc*x)'};
    Fit3_text(2) = {['A = ', num2str(ThreeExp_fit.A), '     ( ', num2str(ThreeExp_CI(1,1)), ', ', num2str(ThreeExp_CI(2,1)), ' )']};
    Fit3_text(3) = {['1/ka = ', num2str(1/ThreeExp_fit.ka), 's','     ( ', num2str(1/ThreeExp_CI(1,4)), ', ', num2str(1/ThreeExp_CI(2,4)), ' )']};
    Fit3_text(4) = {['B = ', num2str(ThreeExp_fit.B), '     ( ', num2str(ThreeExp_CI(1,2)), ', ', num2str(ThreeExp_CI(2,2)), ' )']};
    Fit3_text(5) = {['1/kb = ', num2str(1/ThreeExp_fit.kb), 's','     ( ', num2str(1/ThreeExp_CI(1,5)), ', ', num2str(1/ThreeExp_CI(2,5)), ' )']};
    Fit3_text(6) = {['C = ', num2str(ThreeExp_fit.B), '     ( ', num2str(ThreeExp_CI(1,3)), ', ', num2str(ThreeExp_CI(2,3)), ' )']};
    Fit3_text(7) = {['1/kc = ', num2str(1/ThreeExp_fit.kb), 's','     ( ', num2str(1/ThreeExp_CI(1,6)), ', ', num2str(1/ThreeExp_CI(2,6)), ' )']};
    
    % this is for the Fit2_log
    Fit2_text(7) = {'2-Exp fit - log-sampled : 1-A*exp(-ka*t)-B*exp(-kb*t)'};
    Fit2_text(8) = {['A = ', num2str(TwoExp_fit_log.A), '     ( ', num2str(TwoExp_CI_log(1,1)), ', ', num2str(TwoExp_CI_log(2,1)), ' )']};
    Fit2_text(9) = {['1/ka = ', num2str(1/TwoExp_fit_log.ka), 's','     ( ', num2str(1/TwoExp_CI_log(1,3)), ', ', num2str(1/TwoExp_CI_log(2,3)), ' )']};
    Fit2_text(10) = {['B = ', num2str(TwoExp_fit_log.B), '     ( ', num2str(TwoExp_CI_log(1,2)), ', ', num2str(TwoExp_CI_log(2,2)), ' )']};
    Fit2_text(11) = {['1/kb = ', num2str(1/TwoExp_fit_log.kb), 's','     ( ', num2str(1/TwoExp_CI_log(1,4)), ', ', num2str(1/TwoExp_CI_log(2,4)), ' )']};

    
    % ----- Step 3: Plotting 
    % ----- OneExp, TwoExp, ThreeExp, TwoExp_log  
    % ----- Each uses function, f1, f2, f3, f2_log  respectively    
    
    % plot 1-Exp, 2-Exp 
    figure('position',[1200 700 900 600]); %[x y width height]
    hold on;
    
    plot(time, meanFRAP, 'o', 'Color', [237/255, 28/255, 36/255], 'MarkerSize', 6);
    plot(xFit1, yFit1, 'k--', 'LineWidth', 3);
    plot(xFit2, yFit2, 'k-', 'LineWidth', 2);
    plot(xFit2_log, yFit2_log, 'k-.', 'LineWidth', 2);
    plot([min(time), max(time)], [1, 1], 'k--', 'LineWidth', 1);
    if LongMovie == 1
        text(5,0.2,Fit1_text,'HorizontalAlignment','Left', 'FontSize',12, 'FontName', 'Helvetica'); %text(x,y,txt) x, y is based on the axis value
        text(20,0.6,Fit2_text,'HorizontalAlignment','Left', 'FontSize',12, 'FontName', 'Helvetica');
    else
        text_x = time(round(length(time)/3));  % added to automatically adjust the position of text. May 2019, by Wulan
        text(text_x,0.6,Fit1_text,'HorizontalAlignment','Left', 'FontSize',14, 'FontName', 'Helvetica');%text(x,y,txt) x, y is based on the axis value
        text(text_x,0.3,Fit2_text,'HorizontalAlignment','Left', 'FontSize',14, 'FontName', 'Helvetica');
    end
    title(strcat({'Exponential fit: TMR FRAP for '},ExpName), 'FontSize',14, 'FontName', 'Helvetica'); % modified by Wulan Nov. 2017
    ylabel('FRAP recovery', 'FontSize',14, 'FontName', 'Helvetica');
    xlabel('Time (seconds)', 'FontSize',14, 'FontName', 'Helvetica');
    axis([min(time) max(time) 0 1.08]);
    legend({'FRAP Data', '1-Exp model fit', '2-Exp model fit', 'log sample 2Exp fit'}, 'Location', 'SouthEast', 'FontSize',12, 'FontName', 'Helvetica'); % modified to bigger font by Wulan Nov. 2017
    legend boxoff;
    hold off;

    % Save as PDF
    if SavePDF == 1
        h=gcf;
        set(h,'Units','Inches');
        pos = get(h,'Position');
        %pos(1) and pos(2) are the coordinate of the legend windows
        %pos(3) and pos(4) are width and height of this windows
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(h,[output_path, strcat(ExpName, '_figure3'), '.pdf'], '-dpdf', '-r0');
    end
    

    % --  PLOT THE 3-EXP FITS (Optional) %%%
    % this section is optional, only if you want to check how it would be
    % looking like if it is fit with three reaction-domainant model, which
    % in most of the cases is not interpretable. So only do it out of
    % curiosity, don't over interpretate. 
    
    figure('position',[1200 100 900 600]); %[x y width height]
    hold on;
    plot(time, meanFRAP, 'o', 'Color', [237/255, 28/255, 36/255], 'MarkerSize', 6);
    plot(xFit3, yFit3, 'k--', 'LineWidth', 3);
    plot([min(time), max(time)], [1, 1], 'k--', 'LineWidth', 1);
    text(text_x,0.6,Fit3_text,'HorizontalAlignment','Left', 'FontSize',14, 'FontName', 'Helvetica');  
    % text_x = time(round(length(time)/3));  % added to automatically adjust the position of text. May 2019, by Wulan
    title(strcat({'3-Exponential fit: TMR FRAP for '},ExpName), 'FontSize',14, 'FontName', 'Helvetica'); % modified by Wulan Nov. 2017
    ylabel('FRAP recovery', 'FontSize',14, 'FontName', 'Helvetica');
    xlabel('Time (seconds)', 'FontSize',14, 'FontName', 'Helvetica');
    axis([min(time) max(time) 0 1.08]);
    hold off;

    % Save as PDF 
    if SavePDF == 1
        h=gcf;
        set(h,'Units','Inches');
        pos = get(h,'Position');
        %pos(1) and pos(2) are the coordinate of the legend windows
        %pos(3) and pos(4) are width and height of this windows
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(h,[output_path, strcat(ExpName, '_figure4'), '.pdf'], '-dpdf', '-r0');
    end
    
    
end

%% Save variables  
    % ----- OneExp, TwoExp, ThreeExp, TwoExp_log      
    % add April 2020
    % save time, meanFRAP, stdFRAP
    save([output_path, ExpName, '.mat'],...  
        'time', 'meanFRAP', 'stdFRAP', 'ExpFactor', 'FrameRate','CircleRadius','ExpName', 'FirstFitFrame', 'LastFitFrame', ...     % merged raw data of the same kind of cell
        'OneExp_CI', 'OneExp_fit', 'OneExp_param', 'TwoExp_CI', 'TwoExp_CI_log', 'TwoExp_fit', 'TwoExp_fit_log', 'TwoExp_param', 'TwoExp_param_log',...
        'xFit1', 'yFit1', 'xFit2', 'yFit2', 'xFit2_log', 'yFit2_log',...  % fitting curves x and y of each kind of fit 
            'Fit1_text', 'Fit2_text', 'Fit3_text');   % description of fitting results

%

