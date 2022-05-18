% a3_PLOT_MultipleFRAPcurves_shadedErrorBar_Wulan_modelfitcurve.m
% Wulan Deng, April 7th, 2020

% Purpose: to compare frap curves from multiple samples/conditions, and visualize them on one figure. 
% The input of this script is output mat file from
% "a2_AnalyzeQuantifiedFRAPdata_v3_Wulan_v2.mat". It contains time,
% meanFRAP, stdFRAP, and fit curve of exp1, exp2, logExp2. The output of this script is figure. 

% this script call the function shadedErrorBar (add-on) 

clear; clc; close all; 


%% Path, Keywords, Parameters (CHANGE)
% output_path
output_path = '/Users/wdeng/Downloads/';

% to load input samples into a structure, you have dir first, and then load
% input_path = ('/Users/wdeng/OneDrive - 北京大学生物医学前沿创新中心/BigRawData2020/2017-2019FRAPdata/mat_FRAP_2020April/'); 
% % each sample has its own folder
% s1 = dir([input_path , 'KI14_DE_2Hz_long330/', 'combine_analysis/', '*.mat']);  % structure 
% s2 = dir([input_path , 'dDBD-ES/', 'combine_analysis/', '*.mat']) ; 
% % or you can dir directly
s1 = dir('/Users/wdeng/Downloads/MergedFRAP_KI14_DE_2Hz_330/MergedFRAP_KI14_DE_2Hz_330.mat');
s1 = dir('/Users/wdeng/Downloads/MergedFRAP_dCTD_DE_2Hz_260/MergedFRAP_dCTD_DE_2Hz_260.mat');
s2 = dir('/Users/wdeng/Downloads/MergedFRAP_dCdN_DE_2Hz_200/MergedFRAP_dCdN_DE_2Hz_200.mat');
output_name = 'Frap_KI14DE_vs_dCTD_fit';

% load mat file and set legend names
s1_stru = load([s1.folder, '/', s1.name]); s1_stru.legend = 'FOXA2-Halo DE';  
s2_stru = load([s2.folder, '/', s2.name]); s2_stru.legend = 'dCTD-Halo DE';  %'exo \DeltaDBD-Halo'
%s3_stru = load([s3.folder, '/', s3.name]); s3_stru.legend =  ;
%s4_stru = load([s4.folder, '/', s4.name]); s4_stru.legend =  ;
%s5_stru = load([s5.folder, '/', s5.name]); s5_stru.legend =  ;

color_plate = {[0.8510 0.3255 0.3098], [0.3608 0.7216 0.3608], [0.8510 0.3255 0.3098]};  % light red and green
%color_plate = {[0.8510 0.3255 0.3098], [0.2588 0.5451 0.7922]};  % light red and blue

%% plot with matlab default errorbar function. 
figure;
% plot the merged raw data
errorbar(s1_stru.time, s1_stru.meanFRAP, s1_stru.stdFRAP, '.', 'MarkerSize', 15, 'color',color_plate{1});
hold on;
errorbar(s2_stru.time, s2_stru.meanFRAP, s2_stru.stdFRAP, '.', 'MarkerSize', 15, 'color',color_plate{2});
hold on;
% plot fitted curves 
plot(s1_stru.xFit1, s1_stru.yFit1, 'k--', 'LineWidth', 2);
plot(s1_stru.xFit2, s1_stru.yFit2, 'k-', 'LineWidth', 2);
plot(s1_stru.xFit2_log, s1_stru.yFit2_log, 'k-.', 'LineWidth', 2);
plot([min(s1_stru.time), max(s1_stru.time)], [1, 1], 'k--', 'LineWidth', 1);
% overlay annotation text

text_x = s1_stru.time(round(length(s1_stru.time)/3));  % added to automatically adjust the position of text. May 2019, by Wulan

plottext = {'2-Exp fit: 1-A*exp(-ka*t)-B*exp(-kb*t)', 'Ts = 1/kb = 122.6 s', '95% CI: [116.2; 129.8]'};  % this is from a2 script.  \tau
text(text_x, 0.7, plottext, 'FontSize', 20); 
hold off;

ax = gca; ax.FontSize = 20; ax.LineWidth = 2; ax.Box = 'on'; ax.XLim(1) =  min(s1_stru.time); ax.XLim(2) = 160;
title('\rm FRAP dynamics (2 Hz)', 'FontSize',20, 'FontName', 'Helvetica'); % modified by Wulan Nov. 2017
ylabel('Relative Fluorescence Intensity', 'FontSize',20, 'FontName', 'Helvetica');
xlabel('Recovery time (sec)', 'FontSize',20, 'FontName', 'Helvetica');
[~, objh] = legend({s1_stru.legend, s2_stru.legend, 'model fit'}, 'Location', 'South',  'FontSize',20, 'FontName', 'Helvetica'); % modified to bigger font by Wulan Nov. 2017
legend boxoff; objhl = findobj(objh, 'type', 'line'); set(objhl, 'Markersize', 30)  % note that even if you plot(x,y,'.') it's a "line" plot 

% save figure
h = gcf; h.OuterPosition = [ 100 100 700 500]; % [left bottom width height]
set(h, 'Units', 'Inches');
pos = get(h, 'Position')'; % pos(3) and pos(4) are the width and heigth of this window
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(h, [output_path, output_name, '_1.pdf'], '-dpdf', '-r0');

hold off;
%% plot with the shadedErrorBar function. It is a function written by others

figure;
s = shadedErrorBar(s1_stru.time, s1_stru.meanFRAP, s1_stru.stdFRAP, 'lineProps',{'-o', 'MarkerFaceColor', color_plate{1}, 'color', color_plate{1}, 'linewidth', 1}, 'transparent', 1); set(s.edge, 'color', 'none');
hold on;
s = shadedErrorBar(s2_stru.time, s2_stru.meanFRAP, s2_stru.stdFRAP, 'lineProps',{'-o', 'MarkerFaceColor', color_plate{2}, 'color', color_plate{2}, 'linewidth', 1},'transparent', 1); set(s.edge, 'color', 'none');% {'b-', 'linewidth', 3, 'markerfacecolor', 'blue'});  % 'b-o'
plot([min(s1_stru.time), max(s1_stru.time)], [1, 1], 'k--', 'LineWidth', 1);
hold on; 
plot(s1_stru.xFit2, s1_stru.yFit2, 'k-', 'LineWidth', 2);
text(text_x, 0.7, plottext, 'FontSize', 20); 
        %text(text_x,0.8,s1_stru.Fit1_text,'HorizontalAlignment','Left', 'FontSize',10, 'FontName', 'Helvetica');%text(x,y,txt) x, y is based on the axis value
        %text(text_x,0.5,s1_stru.Fit2_text,'HorizontalAlignment','Left', 'FontSize',10, 'FontName', 'Helvetica'); %1:5 is about Fit2, 6:11 is aobut Fit2-log sampled

hold off; 
% figure cosmetics
ax = gca; ax.FontSize = 16; ax.LineWidth = 2; ax.Box = 'on'; ax.XLim(1) =  min(s1_stru.time); ax.XLim(2) = 160;
title('\rm FRAP dynamics (2 Hz)', 'FontSize',20, 'FontName', 'Helvetica'); % modified by Wulan Nov. 2017
ylabel('Relative Fluorescence Intensity', 'FontSize',20, 'FontName', 'Helvetica');
xlabel('Time (sec)', 'FontSize',20, 'FontName', 'Helvetica');
legend({s1_stru.legend, s2_stru.legend, 'model fit'}, 'Location', 'South', 'FontSize',20, 'FontName', 'Helvetica'); % modified to bigger font by Wulan Nov. 2017
legend boxoff;

% save figure
h = gcf;h.OuterPosition = [ 100 100 700 500];
set(h, 'Units', 'Inches');  % [left bottom width height]
pos = get(h, 'Position')'; % pos(3) and pos(4) are the width and heigth of this window
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(h, [output_path, output_name, '_2.pdf'], '-dpdf', '-r0');

hold off;



