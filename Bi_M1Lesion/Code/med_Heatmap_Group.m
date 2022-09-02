function med_Heatmap_Group(btAll2d,plotRange,trialNearBound)
%%
global lesionBoundary edges_RT edges_HT edges_RelT smo_win;
altMethod = {'mean','median','geomean'};
cenMethod = altMethod{2};
grandCen = 'mean';
grandErr = 'sem';

edges_RT = 0:0.025:0.6; % Reaction Time (only correct)
edges_RelT = 0:0.05:3; % Realease Time (correct + late trials)
edges_HT = 0:0.05:3.5; % Hold Time (all trials)

smo_win = 8; % smoothdata('gaussian'), 
prdName = ["Pre","Post"];
fplist = [0.5 1.0 1.5];
%% Data processing
btAll2d_use = {};
for i=1:size(plotRange,1)
    if find(cellfun(@(x)~isempty(x),btAll2d(i,:)),1,'last') < plotRange{i,1}(end)
        plotRange{i,1} = plotRange{i,1}(1):find(cellfun(@(x)~isempty(x),btAll2d(i,:)),1,'last');
    end
    btAll2d_use(i,1:length(btAll2d(i,plotRange{i,1}))) = btAll2d(i,plotRange{i,1});
end

EES = struct;
[EES.SBS,TBT] = packData(btAll2d_use,cenMethod);
EES.Pre = estTBT_3FPs(TBT(TBT.Task=="ThreeFPsMixedBpod" & TBT.Group=="Pre",:),cenMethod);
EES.Post = estTBT_3FPs(TBT(TBT.Task=="ThreeFPsMixedBpod" & TBT.Group=="Post",:),cenMethod);
%EES.Pre.Group = strcat(EES.Pre.Group,'-Pre');
%EES.Post.Group = strcat(EES.Post.Group,'-Post');
EES.Pre = addvars(EES.Pre,repelem("Pre",length(EES.Pre.Subject))','Before','Task','NewVariableNames','Group');
EES.Post = addvars(EES.Post,repelem("Post",length(EES.Post.Subject))','Before','Task','NewVariableNames','Group');
EES.TBT = [EES.Pre;EES.Post];

EES_val = removevars(EES.TBT,{'Subject','Task'});
GAS = grpstats(EES_val,'Group',{grandCen,grandErr});

% lession ± ntrials
EES_trial = struct;
[EES_trial.SBS,TBT_trial] = packData(btAll2d,cenMethod);
%TBT_trial_use_pre = TBT_trial(find(TBT_trial.Task=="ThreeFPsMixedBpod" & TBT_trial.Type~=0 & (TBT_trial.FP==0.5 | TBT_trial.FP==1.0 | TBT_trial.FP==1.5) & TBT_trial.Group=="Pre",trialNearBound,'last'),:);
%Session_pre = unique(TBT_trial_use_pre.Session,'stable'); 
%TBT_trial_sessionBoundary_pre = TBT_trial(find(TBT_trial.Task=="ThreeFPsMixedBpod" & TBT_trial.Type~=0 & (TBT_trial.FP==0.5 | TBT_trial.FP==1.0 | TBT_trial.FP==1.5) & TBT_trial.Session==Session_pre(1)),:);
%nSession_pre = length(Session_pre)-1+sum(TBT_trial_use_pre.Session==Session_pre(1))/size(TBT_trial_sessionBoundary_pre,1);
%TBT_trial_use_post = TBT_trial(find(TBT_trial.Task=="ThreeFPsMixedBpod" & TBT_trial.Type~=0 & (TBT_trial.FP==0.5 | TBT_trial.FP==1.0 | TBT_trial.FP==1.5) & TBT_trial.Group=="Post",trialNearBound,'first'),:);
%Session_post = unique(TBT_trial_use_post.Session,'stable'); 
%TBT_trial_sessionBoundary_post = TBT_trial(find(TBT_trial.Task=="ThreeFPsMixedBpod" & TBT_trial.Type~=0 & (TBT_trial.FP==0.5 | TBT_trial.FP==1.0 | TBT_trial.FP==1.5) & TBT_trial.Session==Session_post(end)),:);
%nSession_post = length(Session_post)-1+sum(TBT_trial_use_post.Session==Session_post(end))/size(TBT_trial_sessionBoundary_post,1);
[TBT_trial_use_pre,nSession_pre,TBT_trial_use_post,nSession_post] = calculateTrials(TBT_trial,trialNearBound);
EES_trial.Pre = estTBT_3FPs(TBT_trial_use_pre,cenMethod);
EES_trial.Post = estTBT_3FPs(TBT_trial_use_post,cenMethod);
%EES.Pre.Group = strcat(EES.Pre.Group,'-Pre');
%EES.Post.Group = strcat(EES.Post.Group,'-Post');
EES_trial.Pre = addvars(EES_trial.Pre,repelem("Pre",length(EES_trial.Pre.Subject))','Before','Task','NewVariableNames','Group');
EES_trial.Post = addvars(EES_trial.Post,repelem("Post",length(EES_trial.Post.Subject))','Before','Task','NewVariableNames','Group');
EES_trial.TBT = [EES_trial.Pre;EES_trial.Post];

EES_trial_val = removevars(EES_trial.TBT,{'Subject','Task'});
GAS_trial = grpstats(EES_trial_val,'Group',{grandCen,grandErr});

xedges = struct;
xedges.edges_RT = edges_RT;
xedges.edges_RelT = edges_RelT;
xedges.edges_HT = edges_HT;
xedges.RT = movmean(edges_RT,2,'Endpoints','discard');
xedges.RelT = movmean(edges_RelT,2,'Endpoints','discard');
xedges.HT = movmean(edges_HT,2,'Endpoints','discard');

save('VarsToPlot.mat','EES','GAS','EES_trial','GAS_trial');
%% Plot
load('VarsToPlot.mat');

cTab20 = [0.0901960784313726,0.466666666666667,0.701960784313725;0.682352941176471,0.780392156862745,0.901960784313726;0.960784313725490,0.498039215686275,0.137254901960784;0.988235294117647,0.729411764705882,0.470588235294118;0.152941176470588,0.631372549019608,0.278431372549020;0.611764705882353,0.811764705882353,0.533333333333333;0.843137254901961,0.149019607843137,0.172549019607843;0.964705882352941,0.588235294117647,0.592156862745098;0.564705882352941,0.403921568627451,0.674509803921569;0.768627450980392,0.690196078431373,0.827450980392157;0.549019607843137,0.337254901960784,0.290196078431373;0.768627450980392,0.607843137254902,0.576470588235294;0.847058823529412,0.474509803921569,0.698039215686275;0.956862745098039,0.709803921568628,0.807843137254902;0.501960784313726,0.501960784313726,0.501960784313726;0.780392156862745,0.780392156862745,0.776470588235294;0.737254901960784,0.745098039215686,0.196078431372549;0.854901960784314,0.862745098039216,0.549019607843137;0.113725490196078,0.737254901960784,0.803921568627451;0.627450980392157,0.843137254901961,0.890196078431373];
cRed = cTab20(7,:);
cRed2 = cTab20(8,:);
cGreen = cTab20(5,:);
cGreen2 = cTab20(6,:);
cBlue = cTab20(1,:);
cBlue2 = cTab20(2,:);
cGray = cTab20(15,:);
cGray2 = cTab20(16,:);
cOrange = cTab20(3,:);
cOrange2 = cTab20(4,:);

colormap_diff = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});

SBS = EES.SBS;
%SBSs = grpstats(removevars(SBS,{'Subject','Date','Task'}),{'Group','Session'},{cenMethod,errMethod});
%Period = repelem("NaN",size(SBS,1))';
%Period(SBS.Date<=lesionBoundary) = "Pre";
%Period(SBS.Date>lesionBoundary) = "Post";
%SBSp = addvars(SBS,Period,'After','Subject');                SBSp==SBS
SBSp = grpstats(removevars(SBS,{'Session','Date','Task'}),{'Subject','Period'},{grandCen,grandErr});

%hf = figure(44); clf(hf,'reset');
%set(hf, 'name', 'Lesion effect', 'units', 'centimeters', 'position', [1 1 17 11.5],...
%    'PaperPositionMode', 'auto','renderer','painter'); % 生科论文要求版面裁掉边距还剩，宽度15.8cm,高度24.2cm
hf = figure(815); clf(hf,'reset');
%set(hf, 'name', 'Distribution heatmap', 'unit', 'centimeters', 'position',[1 1 12.5 13], 'paperpositionmode', 'auto' ,'renderer','painter');
set(hf, 'name', 'Distribution heatmap', 'unit', 'centimeters', 'position',[1 1 14.5 13], 'paperpositionmode', 'auto' ,'renderer','painter');

size1 = [4 5];
% ys = [1 5.2 8.8 12.5 16.3,20.2]; % yStart
ys = [7.5 1]; % yStart
%xs1 = [1.5 7]; % xStart
%xs1 = [1.5 7 13];
%xs1 = [1.5 1.5 8];
xs1 = [1.5 2 9];



% RT include late/Release time
%ha1d = axes('unit', 'centimeters', 'position', [xs1(1) ys(1) size1], 'nextplot', 'add',...
%   'xlim',[xedges.edges_RelT(1),1],'xtick',[xedges.edges_RelT(1):0.25:1],...
%    'xticklabel',cellstr(string((xedges.edges_RelT(1):0.25:1)*1000)),'fontsize',7,'fontname','Arial',...
%    'ylim', [0.5 size(EES_trial.TBT, 1)+0.5],'ytick', [0.5:size(EES_trial.TBT, 1)+0.5]+0.5, 'yticklabel', {'Pre', 'Post'}, 'ticklength', [.02 .025]);

% do a little upsamnpling x10 
%New_EES_trial_TBT_RelTdist = zeros(size(EES_trial.TBT.RelTdist));
%New_EES_trial_TBT_RelTdist(1:2:end-1,:) = EES_trial.Pre.RelTdist;
%New_EES_trial_TBT_RelTdist(2:2:end,:) = EES_trial.Post.RelTdist;
%NewCenters_RelT = linspace(xedges.edges_RelT(1), xedges.edges_RelT(end-1), 10*length(xedges.RelT));
%NewDistribution_RelT = interp1(xedges.edges_RelT(1:end-1), New_EES_trial_TBT_RelTdist', NewCenters_RelT, 'v5cubic');
%NewDistribution_RelT = NewDistribution_RelT';
%himageRelT = imagesc(NewCenters_RelT, [1:size(NewDistribution_RelT, 1)], NewDistribution_RelT, [0 max(NewDistribution_RelT(:))]);
%colormap(ha1d, fake_parula);
%colorbar('Units','centimeters','Position',[xs1(1)+size1(1)+0.04 ys(1) 0.25 size1(2)], 'FontSize', 6,'FontName','Arial');

%for i = 2:size(EES_trial.Pre, 1)
%    line([0 3], 0.5+2*(i-1)*[1 1], 'color', 'r', 'linewidth', 1, 'linestyle', '-')
%     text(2300, 2*i-0.5, Animalnames{i}, 'fontsize', 10, 'color', 'w');
%end
%for i = 1:size(EES_trial.Pre, 1)
%    text(0.7, 2*i-0.5, EES_trial.Pre.Subject{i}, 'fontsize', 7,'fontname','Arial', 'color', 'w');
%end
%xlabel('RT include late (ms)','fontname','Arial', 'fontsize', 7);
%title('Total','fontname','Arial', 'fontsize', 7);

% RT include late/Release time    FP = 0.5
%ha2d = axes('unit', 'centimeters', 'position', [xs1(2) ys(1) size1], 'nextplot', 'add',...
%   'xlim',[xedges.edges_RelT(1),1],'xtick',[xedges.edges_RelT(1):0.25:1],...
%    'xticklabel',cellstr(string((xedges.edges_RelT(1):0.25:1)*1000)),'fontsize',7,'fontname','Arial',...
%    'ylim', [0.5 size(EES_trial.TBT, 1)+0.5],'ytick', [0.5:size(EES_trial.TBT, 1)+0.5]+0.5, 'yticklabel', {'Pre', 'Post'}, 'ticklength', [.02 .025]);

% do a little upsamnpling x10 
%New_EES_trial_TBT_RelTdist_S = zeros(size(EES_trial.TBT.RelTdist_S));
%New_EES_trial_TBT_RelTdist_S(1:2:end-1,:) = EES_trial.Pre.RelTdist_S;
%New_EES_trial_TBT_RelTdist_S(2:2:end,:) = EES_trial.Post.RelTdist_S;
%NewCenters_RelT_S = linspace(xedges.edges_RelT(1), xedges.edges_RelT(end-1), 10*length(xedges.RelT));
%NewDistribution_RelT_S = interp1(xedges.edges_RelT(1:end-1), New_EES_trial_TBT_RelTdist_S', NewCenters_RelT_S, 'v5cubic');
%NewDistribution_RelT_S = NewDistribution_RelT_S';
%himageRelT_S = imagesc(NewCenters_RelT_S, [1:size(NewDistribution_RelT_S, 1)], NewDistribution_RelT_S, [0 max(NewDistribution_RelT_S(:))]);
%colormap(ha2d, fake_parula);

%for i = 2:size(EES_trial.Pre, 1)
%    line([0 3], 0.5+2*(i-1)*[1 1], 'color', 'r', 'linewidth', 1, 'linestyle', '-')
%     text(2300, 2*i-0.5, Animalnames{i}, 'fontsize', 10, 'color', 'w');
%end
%for i = 1:size(EES_trial.Pre, 1)
%    text(0.7, 2*i-0.5, EES_trial.Pre.Subject{i}, 'fontsize', 7,'fontname','Arial', 'color', 'w');
%end
%xlabel('RT include late (ms)','fontname','Arial', 'fontsize', 7);
%title('FP 0.5 s','fontname','Arial', 'fontsize', 7);

% RT include late/Release time    Total   sorted：median RelT
ha2d = axes('unit', 'centimeters', 'position', [xs1(2) ys(1) size1], 'nextplot', 'add',...
   'xlim',[xedges.edges_RelT(1),1],'xtick',[xedges.edges_RelT(1):0.25:1],...
    'xticklabel',cellstr(string((xedges.edges_RelT(1):0.25:1)*1000)),'fontsize',7,'fontname','Arial',...
    'ylim', [0.5 size(EES_trial.TBT, 1)+0.5],'ytick', [prctile(1:size(EES_trial.TBT, 1),25) prctile(1:size(EES_trial.TBT, 1),75)], 'yticklabel', {'Post', 'Pre'}, 'ticklength', [.02 .025]);

% do a little upsamnpling x10 
[Sorted_EES_trial_Pre,idx_sorted] = sortrows(EES_trial.Pre,'RelT','descend');
Sorted_EES_trial_Post = EES_trial.Post(idx_sorted,:);
Sorted_EES_trial_TBT = [Sorted_EES_trial_Post;Sorted_EES_trial_Pre];
New_EES_trial_TBT_RelTdist = Sorted_EES_trial_TBT.RelTdist;
NewCenters_RelT = linspace(xedges.edges_RelT(1), xedges.edges_RelT(end-1), 10*length(xedges.RelT));
NewDistribution_RelT = interp1(xedges.edges_RelT(1:end-1), New_EES_trial_TBT_RelTdist', NewCenters_RelT, 'v5cubic');
NewDistribution_RelT = NewDistribution_RelT';
himageRelT = imagesc(NewCenters_RelT, [1:size(NewDistribution_RelT, 1)], NewDistribution_RelT, [0 max(NewDistribution_RelT(:))]);
colormap(ha2d, fake_parula);
cb2 = colorbar('Units','centimeters','Position',[xs1(2)-0.25-0.6 ys(1) 0.25 size1(2)], 'FontSize', 6,'FontName','Arial');
cb2.Label.String = 'Probability (%)';
cb2.Label.FontSize = 6;
cb2.Label.FontName = 'Arial';

%for i = 2:size(EES_trial.Pre, 1)
%    line([0 3], 0.5+2*(i-1)*[1 1], 'color', 'r', 'linewidth', 1, 'linestyle', '-')
%     text(2300, 2*i-0.5, Animalnames{i}, 'fontsize', 10, 'color', 'w');
%end
line([0 3], prctile(1:size(EES_trial.TBT, 1),50)*[1 1], 'color', 'r', 'linewidth', 1, 'linestyle', '-')
for i = 1:size(EES_trial.Pre, 1)*2
    text(1.01, i, Sorted_EES_trial_TBT.Subject{i}, 'fontsize', 7,'fontname','Arial', 'color', 'k');
end
xlabel('RT include late (ms)','fontname','Arial', 'fontsize', 7);
title('Total','fontname','Arial', 'fontsize', 7);

% RT include late/Release time(pre-post)
ha5d = axes('unit', 'centimeters', 'position', [xs1(3) ys(1) size1], 'nextplot', 'add',...
   'xlim',[xedges.edges_RelT(1),1],'xtick',[xedges.edges_RelT(1):0.25:1],...
    'xticklabel',cellstr(string((xedges.edges_RelT(1):0.25:1)*1000)),'fontsize',7,'fontname','Arial',...
    'ylim', [0.5 size(EES_trial.Pre, 1)+0.5],'ytick', [], 'yticklabel', {}, 'ticklength', [.02 .025]);

% do a little upsamnpling x10 
New_EES_trial_TBT_RelTdist_PreminusPost = EES_trial.Pre.RelTdist(idx_sorted,:) - EES_trial.Post.RelTdist(idx_sorted,:);
NewCenters_RelT = linspace(xedges.edges_RelT(1), xedges.edges_RelT(end-1), 10*length(xedges.RelT));
NewDistribution_RelT_PreminusPost = interp1(xedges.edges_RelT(1:end-1), New_EES_trial_TBT_RelTdist_PreminusPost', NewCenters_RelT, 'v5cubic');
NewDistribution_RelT_PreminusPost = NewDistribution_RelT_PreminusPost';
himageRelT_PreminusPost = imagesc(NewCenters_RelT, [1:size(NewDistribution_RelT_PreminusPost, 1)], NewDistribution_RelT_PreminusPost, [-max(abs(NewDistribution_RelT_PreminusPost(:))) max(abs(NewDistribution_RelT_PreminusPost(:)))]);
colormap(ha5d, colormap_diff);
cb5 = colorbar('Units','centimeters','Position',[xs1(3)-0.25-0.1 ys(1) 0.25 size1(2)],'FontSize', 6,'FontName','Arial');
cb5.Label.String = 'Probability (%)';
cb5.Label.FontSize = 6;
cb5.Label.FontName = 'Arial';
%cb5.Label.Units = 'centimeters';
%cb5.Label.Position = [cb5.Position(1)+cb5.Position(3)/2, cb5.Position(2)+cb5.Position(4)];
%cb5.Label.Rotation = 0;
%cb5.Label.HorizontalAlignment = 'center';

%for i = 2:size(EES_trial.Pre, 1)
%    line([0 3], 0.5+2*(i-1)*[1 1], 'color', 'r', 'linewidth', 1, 'linestyle', '-')
%     text(2300, 2*i-0.5, Animalnames{i}, 'fontsize', 10, 'color', 'w');
%end
Sorted_EES_trial_Pre_Subject = EES_trial.Pre.Subject(idx_sorted,:);
for i = 1:size(EES_trial.Pre, 1)
    text(1.01, i, Sorted_EES_trial_Pre_Subject{i}, 'fontsize', 7,'fontname','Arial', 'color', 'k');
end
xlabel('RT include late (ms)','fontname','Arial', 'fontsize', 7);
title('Total(Pre - Post)','fontname','Arial', 'fontsize', 7);

% RT include late/Release time    FP = 1.0
%ha3d = axes('unit', 'centimeters', 'position', [xs1(1) ys(2) size1], 'nextplot', 'add',...
%   'xlim',[xedges.edges_RelT(1),1],'xtick',[xedges.edges_RelT(1):0.25:1],...
%    'xticklabel',cellstr(string((xedges.edges_RelT(1):0.25:1)*1000)),'fontsize',7,'fontname','Arial',...
%    'ylim', [0.5 size(EES_trial.TBT, 1)+0.5],'ytick', [0.5:size(EES_trial.TBT, 1)+0.5]+0.5, 'yticklabel', {'Pre', 'Post'}, 'ticklength', [.02 .025]);

% do a little upsamnpling x10 
%New_EES_trial_TBT_RelTdist_M = zeros(size(EES_trial.TBT.RelTdist_M));
%New_EES_trial_TBT_RelTdist_M(1:2:end-1,:) = EES_trial.Pre.RelTdist_M;
%New_EES_trial_TBT_RelTdist_M(2:2:end,:) = EES_trial.Post.RelTdist_M;
%NewCenters_RelT_M = linspace(xedges.edges_RelT(1), xedges.edges_RelT(end-1), 10*length(xedges.RelT));
%NewDistribution_RelT_M = interp1(xedges.edges_RelT(1:end-1), New_EES_trial_TBT_RelTdist_M', NewCenters_RelT_M, 'v5cubic');
%NewDistribution_RelT_M = NewDistribution_RelT_M';
%himageRelT_M = imagesc(NewCenters_RelT_M, [1:size(NewDistribution_RelT_M, 1)], NewDistribution_RelT_M, [0 max(NewDistribution_RelT_M(:))]);
%colormap(ha3d, fake_parula);

%for i = 2:size(EES_trial.Pre, 1)
%    line([0 3], 0.5+2*(i-1)*[1 1], 'color', 'r', 'linewidth', 1, 'linestyle', '-')
%     text(2300, 2*i-0.5, Animalnames{i}, 'fontsize', 10, 'color', 'w');
%end
%for i = 1:size(EES_trial.Pre, 1)
%    text(0.7, 2*i-0.5, EES_trial.Pre.Subject{i}, 'fontsize', 7,'fontname','Arial', 'color', 'w');
%end
%xlabel('RT include late (ms)','fontname','Arial', 'fontsize', 7);
%title('FP 1.0 s','fontname','Arial', 'fontsize', 7);

% Press duration    FP = 1.0
%ha3d = axes('unit', 'centimeters', 'position', [xs1(1) ys(2) size1], 'nextplot', 'add',...
%   'xlim',[xedges.edges_HT(1),3],'xtick',[xedges.edges_HT(1):1:3],...
%    'xticklabel',cellstr(string((xedges.edges_HT(1):1:3)*1000)),'fontsize',7,'fontname','Arial',...
%    'ylim', [0.5 size(EES_trial.TBT, 1)+0.5],'ytick', [0.5:size(EES_trial.TBT, 1)+0.5]+0.5, 'yticklabel', {'Pre', 'Post'}, 'ticklength', [.02 .025]);

% do a little upsamnpling x10 
%New_EES_trial_TBT_HTdist_M = zeros(size(EES_trial.TBT.HTdist_M));
%New_EES_trial_TBT_HTdist_M(1:2:end-1,:) = EES_trial.Pre.HTdist_M;
%New_EES_trial_TBT_HTdist_M(2:2:end,:) = EES_trial.Post.HTdist_M;
%NewCenters_HT_M = linspace(xedges.edges_HT(1), xedges.edges_HT(end-1), 10*length(xedges.HT));
%NewDistribution_HT_M = interp1(xedges.edges_HT(1:end-1), New_EES_trial_TBT_HTdist_M', NewCenters_HT_M, 'v5cubic');
%NewDistribution_HT_M = NewDistribution_HT_M';
%himageHT_M = imagesc(NewCenters_HT_M, [1:size(NewDistribution_HT_M, 1)], NewDistribution_HT_M, [0 max(NewDistribution_HT_M(:))]);
%colormap(ha3d, fake_parula);
%colorbar('Units','centimeters','Position',[xs1(1)+size1(1)+0.04 ys(2) 0.25 size1(2)], 'FontSize', 6,'FontName','Arial');

%for i = 2:size(EES_trial.Pre, 1)
%    line([0 3], 0.5+2*(i-1)*[1 1], 'color', 'r', 'linewidth', 1, 'linestyle', '-')
%     text(2300, 2*i-0.5, Animalnames{i}, 'fontsize', 10, 'color', 'w');
%end
%for i = 1:size(EES_trial.Pre, 1)
%    text(0.7*3, 2*i-0.5, EES_trial.Pre.Subject{i}, 'fontsize', 7,'fontname','Arial', 'color', 'w');
%    line(fplist(2)*[1 1], 0.5+[0 2]+2*(i-1), 'color', 'w', 'linewidth', 2, 'linestyle', ':');
%end
%xlabel('Press duration (ms)','fontname','Arial', 'fontsize', 7);
%title('FP 1.0 s','fontname','Arial', 'fontsize', 7);

% RT include late/Release time    FP = 1.5
%ha4d = axes('unit', 'centimeters', 'position', [xs1(2) ys(2) size1], 'nextplot', 'add',...
%   'xlim',[xedges.edges_RelT(1),1],'xtick',[xedges.edges_RelT(1):0.25:1],...
%    'xticklabel',cellstr(string((xedges.edges_RelT(1):0.25:1)*1000)),'fontsize',7,'fontname','Arial',...
%    'ylim', [0.5 size(EES_trial.TBT, 1)+0.5],'ytick', [0.5:size(EES_trial.TBT, 1)+0.5]+0.5, 'yticklabel', {'Pre', 'Post'}, 'ticklength', [.02 .025]);

% do a little upsamnpling x10 
%New_EES_trial_TBT_RelTdist_L = zeros(size(EES_trial.TBT.RelTdist_L));
%New_EES_trial_TBT_RelTdist_L(1:2:end-1,:) = EES_trial.Pre.RelTdist_L;
%New_EES_trial_TBT_RelTdist_L(2:2:end,:) = EES_trial.Post.RelTdist_L;
%NewCenters_RelT_L = linspace(xedges.edges_RelT(1), xedges.edges_RelT(end-1), 10*length(xedges.RelT));
%NewDistribution_RelT_L = interp1(xedges.edges_RelT(1:end-1), New_EES_trial_TBT_RelTdist_L', NewCenters_RelT_L, 'v5cubic');
%NewDistribution_RelT_L = NewDistribution_RelT_L';
%himageRelT_L = imagesc(NewCenters_RelT_L, [1:size(NewDistribution_RelT_L, 1)], NewDistribution_RelT_L, [0 max(NewDistribution_RelT_L(:))]);
%colormap(ha4d, fake_parula);

%for i = 2:size(EES_trial.Pre, 1)
%    line([0 3], 0.5+2*(i-1)*[1 1], 'color', 'r', 'linewidth', 1, 'linestyle', '-')
%     text(2300, 2*i-0.5, Animalnames{i}, 'fontsize', 10, 'color', 'w');
%end
%for i = 1:size(EES_trial.Pre, 1)
%    text(0.7, 2*i-0.5, EES_trial.Pre.Subject{i}, 'fontsize', 7,'fontname','Arial', 'color', 'w');
%end
%xlabel('RT include late (ms)','fontname','Arial', 'fontsize', 7);
%title('FP 1.5 s','fontname','Arial', 'fontsize', 7);

% Press duration    FP = 1.0   sorted：median HT
ha4d = axes('unit', 'centimeters', 'position', [xs1(2) ys(2) size1], 'nextplot', 'add',...
   'xlim',[xedges.edges_HT(1),3],'xtick',[xedges.edges_HT(1):1:3],...
    'xticklabel',cellstr(string((xedges.edges_HT(1):1:3)*1000)),'fontsize',7,'fontname','Arial',...
    'ylim', [0.5 size(EES_trial.TBT, 1)+0.5],'ytick', [prctile(1:size(EES_trial.TBT, 1),25) prctile(1:size(EES_trial.TBT, 1),75)], 'yticklabel', {'Post', 'Pre'}, 'ticklength', [.02 .025]);

% do a little upsamnpling x10 
[HT_M_Sorted_EES_trial_Pre,idx_sorted_HT_M] = sortrows(EES_trial.Pre,'HT_M','descend');
HT_M_Sorted_EES_trial_Post = EES_trial.Post(idx_sorted_HT_M,:);
HT_M_Sorted_EES_trial_TBT = [HT_M_Sorted_EES_trial_Post;HT_M_Sorted_EES_trial_Pre];
New_EES_trial_TBT_HTdist_M = HT_M_Sorted_EES_trial_TBT.HTdist_M;
NewCenters_HT = linspace(xedges.edges_HT(1), xedges.edges_HT(end-1), 10*length(xedges.HT));
NewDistribution_HT_M = interp1(xedges.edges_HT(1:end-1), New_EES_trial_TBT_HTdist_M', NewCenters_HT, 'v5cubic');
NewDistribution_HT_M = NewDistribution_HT_M';
himageHT_M = imagesc(NewCenters_HT, [1:size(NewDistribution_HT_M, 1)], NewDistribution_HT_M, [0 max(NewDistribution_HT_M(:))]);
colormap(ha4d, fake_parula);
cb4 = colorbar('Units','centimeters','Position',[xs1(2)-0.25-0.6 ys(2) 0.25 size1(2)], 'FontSize', 6,'FontName','Arial');
cb4.Label.String = 'Probability (%)';
cb4.Label.FontSize = 6;
cb4.Label.FontName = 'Arial';

%for i = 2:size(EES_trial.Pre, 1)
%    line([0 3], 0.5+2*(i-1)*[1 1], 'color', 'r', 'linewidth', 1, 'linestyle', '-')
%     text(2300, 2*i-0.5, Animalnames{i}, 'fontsize', 10, 'color', 'w');
%end
line([0 3], prctile(1:size(EES_trial.TBT, 1),50)*[1 1], 'color', 'r', 'linewidth', 1, 'linestyle', '-')
for i = 1:size(EES_trial.Pre, 1)*2
    text(1.01*3, i, HT_M_Sorted_EES_trial_TBT.Subject{i}, 'fontsize', 7,'fontname','Arial', 'color', 'k');
    line(fplist(2)*[1 1], 0.5+[0 1]+1*(i-1), 'color', 'w', 'linewidth', 2, 'linestyle', ':');
end
xlabel('Press duration (ms)','fontname','Arial', 'fontsize', 7);
title('FP 1.0 s','fontname','Arial', 'fontsize', 7);

% Press duration(pre-post)    FP = 1.0
ha6d = axes('unit', 'centimeters', 'position', [xs1(3) ys(2) size1], 'nextplot', 'add',...
   'xlim',[xedges.edges_HT(1),3],'xtick',[xedges.edges_HT(1):1:3],...
    'xticklabel',cellstr(string((xedges.edges_HT(1):1:3)*1000)),'fontsize',7,'fontname','Arial',...
    'ylim', [0.5 size(EES_trial.Pre, 1)+0.5],'ytick', [], 'yticklabel', {}, 'ticklength', [.02 .025]);

% do a little upsamnpling x10 
New_EES_trial_TBT_HTdist_M_PreminusPost = EES_trial.Pre.HTdist_M(idx_sorted_HT_M,:) - EES_trial.Post.HTdist_M(idx_sorted_HT_M,:);
NewCenters_HT_M = linspace(xedges.edges_HT(1), xedges.edges_HT(end-1), 10*length(xedges.HT));
NewDistribution_HT_M_PreminusPost = interp1(xedges.edges_HT(1:end-1), New_EES_trial_TBT_HTdist_M_PreminusPost', NewCenters_HT_M, 'v5cubic');
NewDistribution_HT_M_PreminusPost = NewDistribution_HT_M_PreminusPost';
himageHT_M_PreminusPost = imagesc(NewCenters_HT_M, [1:size(NewDistribution_HT_M_PreminusPost, 1)], NewDistribution_HT_M_PreminusPost, [-max(abs(NewDistribution_HT_M_PreminusPost(:))) max(abs(NewDistribution_HT_M_PreminusPost(:)))]);
colormap(ha6d, colormap_diff);
cb6 = colorbar('Units','centimeters','Position',[xs1(3)-0.25-0.1 ys(2) 0.25 size1(2)],'FontSize', 6,'FontName','Arial');
cb6.Label.String = 'Probability (%)';
cb6.Label.FontSize = 6;
cb6.Label.FontName = 'Arial';

%for i = 2:size(EES_trial.Pre, 1)
%    line([0 3], 0.5+2*(i-1)*[1 1], 'color', 'r', 'linewidth', 1, 'linestyle', '-')
%     text(2300, 2*i-0.5, Animalnames{i}, 'fontsize', 10, 'color', 'w');
%end
Sorted_EES_trial_Pre_Subject = EES_trial.Pre.Subject(idx_sorted_HT_M,:);
for i = 1:size(EES_trial.Pre, 1)
    text(1.01*3, i, Sorted_EES_trial_Pre_Subject{i}, 'fontsize', 7,'fontname','Arial', 'color', 'k');
    line(fplist(2)*[1 1], 0.5+[0 1]+1*(i-1), 'color', 'k', 'linewidth', 2, 'linestyle', ':');
end
xlabel('Press duration (ms)','fontname','Arial', 'fontsize', 7);
title('FP 1.0 s(Pre-Post)','fontname','Arial', 'fontsize', 7);
%% Save
savename = fullfile(pwd,'MiceLesion_RelT_Heatmap');
saveas(hf,savename,'fig');
print(hf,'-dpng',savename);
print(hf,'-depsc2',savename);
end
%% Functions
function [SBS,TBT] = packData(btAll2d,cenMethod)
global lesionBoundary;
SBS = table; % session by session data
TBT = table; % trial by trial data
for i=1:size(btAll2d,1)
    for j=1:size(btAll2d,2)
        T = btAll2d{i,j};
        if isempty(T)
            continue
        end
        if T.Date(1) <= lesionBoundary(i,1)
            Period = "Pre";
        else
            Period = "Post";
        end
        tempT = estSBS(T,j,cenMethod);
        tempTTTT = addvars(tempT,Period,'After','Subject');
        SBS = [SBS;tempTTTT];
        
        nrow = size(T,1);
        if nrow>1
            %tempT = addvars(T,repelem(grp(i),nrow)','Before','Date','NewVariableNames','Group');
            tempTT = addvars(T,repelem(j,nrow)','After','Date','NewVariableNames','Session');
            tempTTT = addvars(tempTT,repelem(Period,nrow)','Before','Date','NewVariableNames','Group');
            TBT = [TBT;tempTTT];
        end
    end
end
end

function outT = estSBS(data,session,cenMethod)

outT = table;
if isempty(data)
    return;
end

sbj = data.Subject(1);
task = data.Task(1);

t = struct;
t.Subject = sbj;
t.Date = data.Date(1);
t.Session = session;
t.Task = task;
tdata = data(data.Type~=0,:); % trial data

% do not eliminate the warm up trials?
t.nTrial = length(tdata.iTrial);
t.Dark   = sum(data.Type==0)./size(data,1);
t.Cor  = sum(tdata.Type==1)./t.nTrial;
t.Pre  = sum(tdata.Type==-1)./t.nTrial;
t.Late = sum(tdata.Type==-2)./t.nTrial;

t.Cor_S = sum(tdata.Type==1 & abs(tdata.FP-0.5)<1e-4)./sum(abs(tdata.FP-0.5)<1e-4);
t.Cor_M = sum(tdata.Type==1 & abs(tdata.FP-1.0)<1e-4)./sum(abs(tdata.FP-1.0)<1e-4);
t.Cor_L = sum(tdata.Type==1 & abs(tdata.FP-1.5)<1e-4)./sum(abs(tdata.FP-1.5)<1e-4);
t.Pre_S = sum(tdata.Type==-1 & abs(tdata.FP-0.5)<1e-4)./sum(abs(tdata.FP-0.5)<1e-4);
t.Pre_M = sum(tdata.Type==-1 & abs(tdata.FP-1.0)<1e-4)./sum(abs(tdata.FP-1.0)<1e-4);
t.Pre_L = sum(tdata.Type==-1 & abs(tdata.FP-1.5)<1e-4)./sum(abs(tdata.FP-1.5)<1e-4);
t.Late_S = sum(tdata.Type==-2 & abs(tdata.FP-0.5)<1e-4)./sum(abs(tdata.FP-0.5)<1e-4);
t.Late_M = sum(tdata.Type==-2 & abs(tdata.FP-1.0)<1e-4)./sum(abs(tdata.FP-1.0)<1e-4);
t.Late_L = sum(tdata.Type==-2 & abs(tdata.FP-1.5)<1e-4)./sum(abs(tdata.FP-1.5)<1e-4);

t.maxFP = max(tdata.FP);
t.t2mFP = find(tdata.FP==t.maxFP,1,'first');
t.minRW = min(tdata.RW);
t.t2mRW = find(tdata.RW==t.minRW,1,'first');

%     if string(taskfilter)=="Wait2"
%         t2cri = find(abs(tdata.FP-1.5)<1e-4 & abs(tdata.RW-0.6)<1e-4,1,'first');
%     else
%         t2cri = find(abs(tdata.FP-1.5)<1e-4,1,'first');
%     end
%     if isempty(t2cri)
%         t2cInv = 0;
%     elseif t2cri==1
%         t2cInv = NaN;
%     else
%         t2cInv = 1./t2cri;
%     end
%     t.t2cInv = t2cInv;

% PressDur(Cor, Pre & Late)? RT(Cor)?
switch cenMethod
    case 'mean'
        %t.HT = mean(tdata.PressDur);
        %t.HT = mean(rmoutliers(tdata.PressDur,'mean'),'omitnan');
        %t.RT = mean(tdata(tdata.Type==1,:).RT);
        %RT = rmoutliers(tdata(tdata.Type==1,:).RT,'mean');
        %t.RT = mean(RT(RT>=0.1),'omitnan');
        [HT_rmvd, HT_idxrmv] = rmoutliers_custome(tdata.PressDur);
        t.HT = mean(HT_rmvd,'omitnan');
        [RT_rmvd, RT_idxrmv] = rmoutliers_custome(tdata(tdata.Type==1,:).RT);
        t.RT = mean(RT_rmvd(RT_rmvd>=0.1),'omitnan');
    case 'median'
        %t.HT = median(tdata.PressDur);
        %t.RT = median(tdata(tdata.Type==1,:).RT);
        %t.HT = median(rmoutliers(tdata.PressDur,'median'),'omitnan');
        %RT = rmoutliers(tdata(tdata.Type==1,:).RT,'median');
        %t.RT = median(RT(RT>=0.1),'omitnan');
        [HT_rmvd, HT_idxrmv] = rmoutliers_custome(tdata.PressDur);
        t.HT = median(HT_rmvd,'omitnan');
        [RT_rmvd, RT_idxrmv] = rmoutliers_custome(tdata(tdata.Type==1,:).RT);
        t.RT = median(RT_rmvd(RT_rmvd>=0.1),'omitnan');
    case 'geomean'
        %t.HT = geomean(tdata.PressDur);
        %t.RT = geomean(tdata(tdata.Type==1 & tdata.RT>0,:).RT);
        %t.HT = geomean(rmoutliers(tdata.PressDur,'quartiles'),'omitnan');
        %RT = rmoutliers(tdata(tdata.Type==1 & tdata.RT>0,:).RT,'quartiles');
        %t.RT = geomean(RT(RT>=0.1),'omitnan');
        [HT_rmvd, HT_idxrmv] = rmoutliers_custome(tdata.PressDur);
        t.HT = geomean(HT_rmvd,'omitnan');
        [RT_rmvd, RT_idxrmv] = rmoutliers_custome(tdata(tdata.Type==1 & tdata.RT>0,:).RT);
        t.RT = geomean(RT_rmvd(RT_rmvd>=0.1),'omitnan');
end
outT = [outT;struct2table(t)];

end

function outT = estTBT_3FPs(TBT,cenMethod)
global edges_RT edges_HT edges_RelT smo_win;
fplist = [0.5,1.0,1.5];
nboot = 1000;

outT = table;
sbjlist = unique(TBT.Subject,'stable');

for i=1:length(sbjlist)
    data = TBT(TBT.Subject==sbjlist(i),:);

    t = struct;
    t.Subject = sbjlist(i);
    t.Task = data.Task(1);
    tdata = data(data.Type~=0,:); % trial data

    t.nSession = length(unique(tdata.Session));
    t.nTrial = size(tdata,1);
    t.Dark = sum(data.Type==0)./(size(data,1));

    % do not eliminate the warm up trials?
    idxFPS = abs(tdata.FP-fplist(1))<1E-4; % small
    idxFPM = abs(tdata.FP-fplist(2))<1E-4; % medium
    idxFPL = abs(tdata.FP-fplist(3))<1E-4; % large
    idxCor = tdata.Type==1;
    idxPre = tdata.Type==-1;
    idxLate = tdata.Type==-2;

    t.Cor = sum(idxCor)./t.nTrial;
    t.Pre = sum(idxPre)./t.nTrial;
    t.Late = sum(idxLate)./t.nTrial;
%         t.Perf = [...
%             sum(idxCor)./t.nTrial,...
%             sum(idxPre)./t.nTrial,...
%             sum(idxLate)./t.nTrial];
    t.Cor_S = sum( idxFPS & idxCor )./sum(idxFPS);
    t.Pre_S = sum( idxFPS & idxPre )./sum(idxFPS);
    t.Late_S = sum( idxFPS & idxLate )./sum(idxFPS);
%         t.Perf_FPS = [...
%             sum( idxFPS & idxCor )./sum(idxFPS),...
%             sum( idxFPS & idxPre )./sum(idxFPS),...
%             sum( idxFPS & idxLate )./sum(idxFPS)];
    t.Cor_M = sum( idxFPM & idxCor )./sum(idxFPM);
    t.Pre_M = sum( idxFPM & idxPre )./sum(idxFPM);
    t.Late_M = sum( idxFPM & idxLate )./sum(idxFPM);
%         t.Perf_FPM = [...
%             sum( idxFPM & idxCor )./sum(idxFPM),...
%             sum( idxFPM & idxPre )./sum(idxFPM),...
%             sum( idxFPM & idxLate )./sum(idxFPM)];
    t.Cor_L = sum( idxFPL & idxCor )./sum(idxFPL);
    t.Pre_L = sum( idxFPL & idxPre )./sum(idxFPL);
    t.Late_L = sum( idxFPL & idxLate )./sum(idxFPL);
%         t.Perf_FPL = [...
%             sum( idxFPL & idxCor )./sum(idxFPL),...
%             sum( idxFPL & idxPre )./sum(idxFPL),...
%             sum( idxFPL & idxLate )./sum(idxFPL)];

%    t.RT = mean(tdata.RT(idxCor),'omitnan');
%   t.RT_sem = std(tdata.RT(idxCor),'omitnan')/sqrt(sum(~isnan(tdata.RT(idxCor))));
%    t.RT_S = mean(tdata.RT(idxCor&idxFPS));
%    t.RT_M = mean(tdata.RT(idxCor&idxFPM));
%    t.RT_L = mean(tdata.RT(idxCor&idxFPL));
    %RT = rmoutliers(tdata.RT(idxCor),cenMethod);
    %RT_S = rmoutliers(tdata.RT(idxCor&idxFPS),cenMethod);
    %RT_M = rmoutliers(tdata.RT(idxCor&idxFPM),cenMethod);
    %RT_L = rmoutliers(tdata.RT(idxCor&idxFPL),cenMethod);
    [RT, RT_idxrmv] = rmoutliers_custome(tdata.RT(idxCor));
    [RT_S, RT_idxrmv_S] = rmoutliers_custome(tdata.RT(idxCor&idxFPS));
    [RT_M, RT_idxrmv_M] = rmoutliers_custome(tdata.RT(idxCor&idxFPM));
    [RT_L, RT_idxrmv_L] = rmoutliers_custome(tdata.RT(idxCor&idxFPL));
    RT = RT(RT>=0.1);
    RT_S = RT_S(RT_S>=0.1);
    RT_M = RT_M(RT_M>=0.1);
    RT_L = RT_L(RT_L>=0.1);
    t.RT = eval(cenMethod+"(RT,'omitnan')");
    RT_CI = eval("bootci(nboot,{@"+cenMethod+",RT},'alpha',0.05)'");
    t.RT_CI = [t.RT-RT_CI(1), RT_CI(2)-t.RT];
    RT_SE = eval("std(bootstrp(nboot,@"+cenMethod+",RT))");
    t.RT_SE = [RT_SE,RT_SE];
    t.RT_S = eval(cenMethod+"(RT_S,'omitnan')");
    t.RT_M = eval(cenMethod+"(RT_M,'omitnan')");
    t.RT_L = eval(cenMethod+"(RT_L,'omitnan')");

    %RelT = rmoutliers(tdata.RT(idxCor|idxLate),cenMethod);
    %RelT_S = rmoutliers(tdata.RT((idxCor|idxLate)&idxFPS),cenMethod);
    %RelT_M = rmoutliers(tdata.RT((idxCor|idxLate)&idxFPM),cenMethod);
    %RelT_L = rmoutliers(tdata.RT((idxCor|idxLate)&idxFPL),cenMethod);
    [RelT, RelT_idxrmv] = rmoutliers_custome(tdata.RT(idxCor|idxLate));
    [RelT_S, RelT_idxrmv_S] = rmoutliers_custome(tdata.RT((idxCor|idxLate)&idxFPS));
    [RelT_M, RelT_idxrmv_M] = rmoutliers_custome(tdata.RT((idxCor|idxLate)&idxFPM));
    [RelT_L, RelT_idxrmv_L] = rmoutliers_custome(tdata.RT((idxCor|idxLate)&idxFPL));
    RelT = RelT(RelT>=0.1);
    RelT_S = RelT_S(RelT_S>=0.1);
    RelT_M = RelT_M(RelT_M>=0.1);
    RelT_L = RelT_L(RelT_L>=0.1);
    t.RelT = eval(cenMethod+"(RelT,'omitnan')");
    RelT_CI = eval("bootci(nboot,{@"+cenMethod+",RelT},'alpha',0.05)'");
    t.RelT_CI = [t.RelT-RelT_CI(1), RelT_CI(2)-t.RelT];
    RelT_SE = eval("std(bootstrp(nboot,@"+cenMethod+",RelT))");
    t.RelT_SE = [RelT_SE,RelT_SE];
    t.RelT_S = eval(cenMethod+"(RelT_S,'omitnan')");
    t.RelT_M = eval(cenMethod+"(RelT_M,'omitnan')");
    t.RelT_L = eval(cenMethod+"(RelT_L,'omitnan')");

    [HT, HT_idxrmv] = rmoutliers_custome(tdata.PressDur);
    [HT_S, HT_idxrmv_S] = rmoutliers_custome(tdata.PressDur(idxFPS));
    [HT_M, HT_idxrmv_M] = rmoutliers_custome(tdata.PressDur(idxFPM));
    [HT_L, HT_idxrmv_L] = rmoutliers_custome(tdata.PressDur(idxFPL));
    t.HT = eval(cenMethod+"(HT,'omitnan')");
    HT_CI = eval("bootci(nboot,{@"+cenMethod+",HT},'alpha',0.05)'");
    t.HT_CI = [t.HT-HT_CI(1), HT_CI(2)-t.HT];
    HT_SE = eval("std(bootstrp(nboot,@"+cenMethod+",HT))");
    t.HT_SE = [HT_SE,HT_SE];
    t.HT_S = eval(cenMethod+"(HT_S,'omitnan')");
    t.HT_M = eval(cenMethod+"(HT_M,'omitnan')");
    t.HT_L = eval(cenMethod+"(HT_L,'omitnan')");
%         t.RT_3FPs = [mean(tdata.RT(idxCor&idxFPS)),mean(tdata.RT(idxCor&idxFPM)),mean(tdata.RT(idxCor&idxFPL))];
%    edges_RT = 0:0.05:0.6;
%    t.RTdist_S = smooth(histcounts(tdata.RT(idxCor&idxFPS),edges_RT,'Normalization','probability'),3)';
%    t.RTdist_M = smooth(histcounts(tdata.RT(idxCor&idxFPM),edges_RT,'Normalization','probability'),3)';
%    t.RTdist_L = smooth(histcounts(tdata.RT(idxCor&idxFPL),edges_RT,'Normalization','probability'),3)';

%    edges_HT = 0:0.05:2.5;
%    t.HTdist = histcounts(tdata.PressDur,edges_HT,'Normalization','probability');
%    t.HTdist_S = smooth(histcounts(tdata.PressDur(idxFPS),edges_HT,'Normalization','probability'),3)';
%    t.HTdist_M = smooth(histcounts(tdata.PressDur(idxFPM),edges_HT,'Normalization','probability'),3)';
%    t.HTdist_L = smooth(histcounts(tdata.PressDur(idxFPL),edges_HT,'Normalization','probability'),3)';
    t.RTdist = smoothdata(histcounts(tdata.RT(idxCor),...
        edges_RT,'Normalization','pdf'),2,'gaussian',smo_win);
    t.RTdist_S = smoothdata(histcounts(tdata.RT(idxCor&idxFPS),...
        edges_RT,'Normalization','pdf'),2,'gaussian',smo_win);
    t.RTdist_M = smoothdata(histcounts(tdata.RT(idxCor&idxFPM),...
        edges_RT,'Normalization','pdf'),2,'gaussian',smo_win);
    t.RTdist_L = smoothdata(histcounts(tdata.RT(idxCor&idxFPL),...
        edges_RT,'Normalization','pdf'),2,'gaussian',smo_win);

    t.HTdist = smoothdata(histcounts(tdata.PressDur,...
        edges_HT,'Normalization','probability'),2,'gaussian',smo_win);
    t.HTdist_S = smoothdata(histcounts(tdata.PressDur(idxFPS),...
        edges_HT,'Normalization','probability'),2,'gaussian',smo_win);
    t.HTdist_M = smoothdata(histcounts(tdata.PressDur(idxFPM),...
        edges_HT,'Normalization','probability'),2,'gaussian',smo_win);
    t.HTdist_L = smoothdata(histcounts(tdata.PressDur(idxFPL),...
        edges_HT,'Normalization','probability'),2,'gaussian',smo_win);

    t.RelTdist = smoothdata(histcounts(tdata.RT(idxCor|idxLate),...
        edges_RelT,'Normalization','probability'),2,'gaussian',smo_win);
    t.RelTdist_S = smoothdata(histcounts(tdata.RT((idxCor|idxLate)&idxFPS),...
        edges_RelT,'Normalization','probability'),2,'gaussian',smo_win);
    t.RelTdist_M = smoothdata(histcounts(tdata.RT((idxCor|idxLate)&idxFPM),...
        edges_RelT,'Normalization','probability'),2,'gaussian',smo_win);
    t.RelTdist_L = smoothdata(histcounts(tdata.RT((idxCor|idxLate)&idxFPL),...
        edges_RelT,'Normalization','probability'),2,'gaussian',smo_win);

    t.RelTdist_CDF = smoothdata(histcounts(tdata.RT(idxCor|idxLate),...
        edges_RelT,'Normalization','cdf'),2,'gaussian',smo_win);
    t.RelTdist_S_CDF = smoothdata(histcounts(tdata.RT((idxCor|idxLate)&idxFPS),...
        edges_RelT,'Normalization','cdf'),2,'gaussian',smo_win);
    t.RelTdist_M_CDF = smoothdata(histcounts(tdata.RT((idxCor|idxLate)&idxFPM),...
        edges_RelT,'Normalization','cdf'),2,'gaussian',smo_win);
    t.RelTdist_L_CDF = smoothdata(histcounts(tdata.RT((idxCor|idxLate)&idxFPL),...
        edges_RelT,'Normalization','cdf'),2,'gaussian',smo_win);

    outT = [outT;struct2table(t)];
end

end

function [TBT_trial_use_pre,nSession_pre,TBT_trial_use_post,nSession_post] = calculateTrials(TBT_trial,trialNearBound)
fplist = [0.5,1.0,1.5];
sbjlist = unique(TBT_trial.Subject,'stable');

TBT_trial_use_pre = [];
nSession_pre = NaN(length(sbjlist),1);
TBT_trial_use_post = [];
nSession_post = NaN(length(sbjlist),1);
for i=1:length(sbjlist)
    data = TBT_trial(TBT_trial.Subject==sbjlist(i),:);
    tempT = data(find(data.Task=="ThreeFPsMixedBpod" & data.Type~=0 & ismember(data.FP,fplist) & data.Group=="Pre",trialNearBound,'last'),:);
    TBT_trial_use_pre = [TBT_trial_use_pre;tempT];
    Session_pre = unique(tempT.Session,'stable'); 
    TBT_trial_sessionBoundary_pre = data(find(data.Task=="ThreeFPsMixedBpod" & data.Type~=0 & ismember(data.FP,fplist) & data.Session==Session_pre(1)),:);
    nSession_pre(i,1) = length(Session_pre)-1+sum(tempT.Session==Session_pre(1))/size(TBT_trial_sessionBoundary_pre,1);

    tempTT = data(find(data.Task=="ThreeFPsMixedBpod" & data.Type~=0 & ismember(data.FP,fplist) & data.Group=="Post",trialNearBound,'first'),:);
    TBT_trial_use_post = [TBT_trial_use_post;tempTT];
    Session_post = unique(tempTT.Session,'stable'); 
    TBT_trial_sessionBoundary_post = data(find(data.Task=="ThreeFPsMixedBpod" & data.Type~=0 & ismember(data.FP,fplist) & data.Session==Session_post(end)),:);
    nSession_post(i,1) = length(Session_post)-1+sum(tempTT.Session==Session_post(end))/size(TBT_trial_sessionBoundary_post,1);
end
end

function symbol = pValue2symbol(p)
    if p>0.05
%         symbol = 'n.s';
        symbol = sprintf('%0.2f',p);
    elseif p>0.01 && p<=0.05
        symbol = '*';
    elseif p>0.001 && p<=0.01
        symbol = '**';
    else
        symbol = '***';
    end
end

function [dataout,indrmv] = rmoutliers_custome(datain)
[data2575] = prctile(datain, [25, 75]);
interq = data2575(2) - data2575(1);
c=5;
indrmv = find(datain>data2575(2)+interq*c | datain<data2575(1)-interq*c);
dataout = datain;
dataout(indrmv) = [];
end
