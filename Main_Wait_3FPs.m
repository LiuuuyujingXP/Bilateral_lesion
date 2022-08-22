% _________________________________________________________________________
% File:               Main_Wait_3FPs.m
% Created on:         Oct 6, 2021
% Created by:         Yu Chen
% Last revised on:    Oct 6, 2021
% Last revised by:    Yu Chen
% _________________________________________________________________________
% Required Functions:
%   med_to_tec_new
%   med_to_protocol
%   med_to_tec_fp
%   med_DataExtract
%   med_LearningPlot_Individual
%   med_3FPsPlot_Individual
% _________________________________________________________________________
% Required Packages:
% 'gramm' by Pierre Morel
% _________________________________________________________________________
% Protocol:
% 1. Add the folder path of this program to MATLAB's search path
% 2. Put all MED files to be analyzed in a specific folder
% 3. Run this code (Main_Wait_3FPs.m)
% 4. Select the folder containing MED files in dialog box
% 5. Extracted data are saved in selected folder as .mat file, and backup 
%   codes & result figures are saved in corresponding subfolders
% _________________________________________________________________________
% Caution:
% 1. Suggest put these codes into a specific folder and don't mixed with
%   data files.
% 2. There are probably some outputs in terminal.
%   - Orphan Press: The presses couldn't be classified into
%       correct/premature/late/inter-trial(dark), and they are finally
%       incorporated into inter-trial type. The index is for all presses.
%   - Mismatch FP or RW: The FP & RW data in Wait 1&2 are inferred by
%       specific rule. If the inferred values are not the same as the
%       recorded values that are calculated by events, the terminal would 
%       show the index of the mismatch value (index of all presses).
%% Initiate
clear;clc;
rng('default'); % Set the random seed for reproducibility of the results
%% Select The Path Containing Data & Saving Results
tarPath = uigetdir2('C:\Users\liuyujing\Me\Work\AnalysisCode','Select the directory containing data file');%选择工作目录

savePath = uigetdir('C:\Users\liuyujing\Me\Work\AnalysisCode','Select the directory to save codes and results');
resultSavePath = fullfile(savePath,'Result');
if ~exist(resultSavePath,'dir')
    mkdir(resultSavePath);
end
addpath(pwd);
%% Extract & Save Data
%addpath(pwd);

btAll2d = {};
global lesionBoundary;
lesionBoundary = NaN(length(tarPath),1);
sessNearBound = 5;
sessNearBound_pre = 5;
sessNearBound_post = 10;
plotRange = cell(length(tarPath),1);
plotRange_SBS = cell(length(tarPath),1);
trialNearBound_LP = 50;
trialNearBound = 500;

for ipath=1:length(tarPath)
    sbjPath = tarPath{ipath};
    dataPath = fullfile(sbjPath,'Pre_Post_Lesion');
    cd(dataPath);
    
    FileNames = arrayfun(@(x)x.name, dir('*Subject*.txt'), 'UniformOutput', false);
    clear bAll;
    btAll = cell(1,length(FileNames));
    for i=1:length(FileNames)
        %bAll(i)=track_training_progress_advanced(FileNames{i}); % classic method     
        %[bAll(i),btAll{i}] = med_DataExtract(FileNames{i},[0,1]); %plot
        [bAll(i),btAll{i}] = med_DataExtract(FileNames{i},[0,0]); %don't plot
    end
    savename = 'bmixedAll_' + upper(btAll{1}.Subject(1));
    save(savename, 'bAll','btAll')

    % Calculate lesionBoundary & plotRange
    Dateincludeyear = [];
    for j=1:length(bAll)
        b = bAll(j);
        Dateincludeyear = [Dateincludeyear; b.Metadata.Date];
    end
    Dateincludeyear = cellstr(Dateincludeyear);
    t = datetime(Dateincludeyear,'InputFormat','yyyyMMdd');
    Maxtimeduration = max(caldays(caldiff(t)));
    sessBeforeles=find(caldays(caldiff(t))==Maxtimeduration,1,'last');
    lesionBoundary(ipath,1) = str2double(Dateincludeyear{sessBeforeles}); %session before lession
    plotRange(ipath,1) = {sessBeforeles-sessNearBound+1:sessBeforeles+sessNearBound};
    plotRange_SBS(ipath,1) = {sessBeforeles-sessNearBound_pre+1:sessBeforeles+sessNearBound_post};
    
    btAll2d(end+1,1:length(btAll)) = btAll;
end

save(fullfile(resultSavePath,'bmixedAllsbj.mat'),'btAll2d','plotRange','plotRange_SBS','trialNearBound','lesionBoundary');
close all;
%% Analysis & Plot
cd(resultSavePath);
load('bmixedAllsbj.mat');
% plotbAll(bAll)
med_Heatmap_Group(btAll2d,plotRange,trialNearBound);
med_LearningPlot_Individual(btAll2d,plotRange_SBS);
med_plotLee2Vigor(btAll2d,plotRange,trialNearBound);
%% Code Backup
codeSavePath = fullfile(savePath,'Code');
if ~exist(codeSavePath,'dir')
    mkdir(codeSavePath);
end

curCodePath = mfilename('fullpath');
[codeFolder,~] = fileparts(curCodePath);

cd(codeFolder);
exeCodes = dir('*.m');
for i=1:size(exeCodes,1)
    copyfile(exeCodes(i).name,codeSavePath)
end
cd(savePath);
%% Function
function [pathname] = uigetdir2(start_path, dialog_title)
% Pick multiple directories and/or files

import javax.swing.JFileChooser;

%if nargin == 0 || start_path == '' || start_path == 0 % Allow a null argument.
    %start_path = pwd;
%end

jchooser = javaObjectEDT('javax.swing.JFileChooser', start_path);

jchooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
if nargin > 1
    jchooser.setDialogTitle(dialog_title);
end

jchooser.setMultiSelectionEnabled(true);

status = jchooser.showOpenDialog([]);

if status == JFileChooser.APPROVE_OPTION
    jFile = jchooser.getSelectedFiles();
	pathname{size(jFile, 1)}=[];
    for i=1:size(jFile, 1)
		pathname{i} = char(jFile(i).getAbsolutePath);
	end
	
elseif status == JFileChooser.CANCEL_OPTION
    pathname = [];
else
    error('Error occured while picking file.');
end
end