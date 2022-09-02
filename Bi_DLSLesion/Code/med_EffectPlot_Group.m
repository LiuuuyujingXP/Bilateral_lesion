function med_EffectPlot_Group(btAll2d,plotRange,trialNearBound)
%%
global lesionBoundary edges_RT edges_HT edges_RelT smo_win;
altMethod = {'mean','median','geomean'};
cenMethod = altMethod{2};
grandCen = 'mean';
grandErr = 'sem';

edges_RT = 0:0.025:0.6; % Reaction Time (only correct)
edges_RelT = 0:0.05:1; % Realease Time (correct + late trials)
edges_HT = 0:0.05:2.5; % Hold Time (all trials)

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

SBS = EES.SBS;
%SBSs = grpstats(removevars(SBS,{'Subject','Date','Task'}),{'Group','Session'},{cenMethod,errMethod});
%Period = repelem("NaN",size(SBS,1))';
%Period(SBS.Date<=lesionBoundary) = "Pre";
%Period(SBS.Date>lesionBoundary) = "Post";
%SBSp = addvars(SBS,Period,'After','Subject');                SBSp==SBS
SBSp = grpstats(removevars(SBS,{'Session','Date','Task'}),{'Subject','Period'},{grandCen,grandErr});

hf = figure(44); clf(hf,'reset');
%set(hf, 'name', 'Lesion effect', 'units', 'centimeters', 'position', [1 1 17 11.5],...
%    'PaperPositionMode', 'auto','renderer','painter'); % 生科论文要求版面裁掉边距还剩，宽度15.8cm,高度24.2cm
set(hf, 'name', 'Lesion effect', 'units', 'centimeters', 'position', [1 1 17 10.9],...
    'PaperPositionMode', 'auto','renderer','painter'); % 生科论文要求版面裁掉边距还剩，宽度15.8cm,高度24.2cm

size1 = [3.5,3.5*0.7];
size2 = [4*0.618,4*0.618];
size3 = [2.4,4*0.618];
size4 = [3.5 3.5*0.7];  
size5 = [3.5*0.7 3.5*0.7]; % compare performance
% ys = [1 5.2 8.8 12.5 16.3,20.2]; % yStart
%ys = [1 8.8 5.2 12.5 16.3,20.2]; % yStart
ys = [4.6 11.8 8.2 1 16.3,20.2]; % yStart
xs1 = [1.5 6 11 15.5]; % xStart
xs2 = [1.5 6.3 11.2 16];
xs3 = [1 5.4 10.5 13.2];



%le1 = legend({'Correct','Premature','Late'},'Fontsize',6,'units','centimeters','Position',[16.13,5,1,1]);% [4.7,2.8,1,1]
%le1.ItemTokenSize = [12,22];
%le1.Position = le1.Position + [0.025 0.045 0 0];
%legend('boxoff');



% subplot3: x:cor/pre/late * Pre/Post * Session/Trial, y:%, line&thickness: S/M/L
ha3 = axes;
set(ha3, 'units', 'centimeters', 'position', [xs2(3) ys(1) size5], 'nextplot', 'add','tickDir', 'out',...
    'xtick',[],'xticklabel',{},'xticklabelRotation',-45,'fontsize',7,'fontname','Arial',...
    'ytick',0:0.5:1,'yticklabel',{'0','50','100'},'xlim', [3.75 7.5],'ylim', [0 1],'ticklength', [0.02 0.025]);
% 'xtick',[1.25,3.25,5.75,7.75,10.25,12.25]
%plot([0.5 1.25],[GAS.mean_Cor_S(GAS.Group=="Pre"),GAS.mean_Cor_S(GAS.Group=="Post")],'-','lineWidth',0.6,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
%plot([0.5 1.25],[GAS.mean_Cor_M(GAS.Group=="Pre"),GAS.mean_Cor_M(GAS.Group=="Post")],'-','lineWidth',1,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
%plot([0.5 1.25],[GAS.mean_Cor_L(GAS.Group=="Pre"),GAS.mean_Cor_L(GAS.Group=="Post")],'-','lineWidth',1.5,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
%plot([1.5 2.25],[GAS.mean_Pre_S(GAS.Group=="Pre"),GAS.mean_Pre_S(GAS.Group=="Post")],'-','lineWidth',0.6,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
%plot([1.5 2.25],[GAS.mean_Pre_M(GAS.Group=="Pre"),GAS.mean_Pre_M(GAS.Group=="Post")],'-','lineWidth',1,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
%plot([1.5 2.25],[GAS.mean_Pre_L(GAS.Group=="Pre"),GAS.mean_Pre_L(GAS.Group=="Post")],'-','lineWidth',1.5,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
%hl1 = plot([2.5 3.25],[GAS.mean_Late_S(GAS.Group=="Pre"),GAS.mean_Late_S(GAS.Group=="Post")],'-','lineWidth',0.6,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');
%hl2 = plot([2.5 3.25],[GAS.mean_Late_M(GAS.Group=="Pre"),GAS.mean_Late_M(GAS.Group=="Post")],'-','lineWidth',1,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');
%hl3 = plot([2.5 3.25],[GAS.mean_Late_L(GAS.Group=="Pre"),GAS.mean_Late_L(GAS.Group=="Post")],'-','lineWidth',1.5,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');
hl1 = plot([4.25 5],[GAS_trial.mean_Cor_S(GAS_trial.Group=="Pre"),GAS_trial.mean_Cor_S(GAS_trial.Group=="Post")],'-','lineWidth',0.6,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
hl2 = plot([4.25 5],[GAS_trial.mean_Cor_M(GAS_trial.Group=="Pre"),GAS_trial.mean_Cor_M(GAS_trial.Group=="Post")],'-','lineWidth',1,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
hl3 = plot([4.25 5],[GAS_trial.mean_Cor_L(GAS_trial.Group=="Pre"),GAS_trial.mean_Cor_L(GAS_trial.Group=="Post")],'-','lineWidth',1.5,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
plot([5.25 6],[GAS_trial.mean_Pre_S(GAS_trial.Group=="Pre"),GAS_trial.mean_Pre_S(GAS_trial.Group=="Post")],'-','lineWidth',0.6,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
plot([5.25 6],[GAS_trial.mean_Pre_M(GAS_trial.Group=="Pre"),GAS_trial.mean_Pre_M(GAS_trial.Group=="Post")],'-','lineWidth',1,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
plot([5.25 6],[GAS_trial.mean_Pre_L(GAS_trial.Group=="Pre"),GAS_trial.mean_Pre_L(GAS_trial.Group=="Post")],'-','lineWidth',1.5,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
plot([6.25 7],[GAS_trial.mean_Late_S(GAS_trial.Group=="Pre"),GAS_trial.mean_Late_S(GAS_trial.Group=="Post")],'-','lineWidth',0.6,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');
plot([6.25 7],[GAS_trial.mean_Late_M(GAS_trial.Group=="Pre"),GAS_trial.mean_Late_M(GAS_trial.Group=="Post")],'-','lineWidth',1,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');
plot([6.25 7],[GAS_trial.mean_Late_L(GAS_trial.Group=="Pre"),GAS_trial.mean_Late_L(GAS_trial.Group=="Post")],'-','lineWidth',1.5,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');


text([4.625,5.625,6.625],repelem(-0.125,3),{'Cor','Pre','Late'},'HorizontalAlignment','center','fontsize',7,'FontName','Arial');
ylabel('Percentage (%)','Fontsize',7,'FontName','Arial');

le1 = legend([hl1,hl2,hl3],{'FP 0.5 s','FP 1.0 s','FP 1.5 s'},'Fontsize',7,'FontName','Arial','units','centimeters',...
    'Position',[xs2(3)+size5(1)+0.55,ys(1)+0.1,1,1]);% [4.7,2.8,1,1]
le1.ItemTokenSize = [12,22];
le1.Position = le1.Position + [0.025 0.045 0 0];
legend('boxoff');


%le3 = legend([hl1,hl2,hl3],{'FP 0.5 s','FP 1.0 s','FP 1.5 s'},...
%    'Fontsize',8,'FontName','Arial','units','centimeters',...
%    'Position',[xs2(3)+size1(1)+0.55,ys(1)+0.1,1,1]); % [10.7,8.9,1,1]
%legend('boxoff');
%le3.ItemTokenSize = [12,22];
%le3.Position = le3.Position + [0.025 0.045 0 0];

% le3 = legend({'Short','Medium','Long'},'fontsize',6,'units','centimeter','Position',[16.1,2.8,1,1]);
% legend('boxoff');
% le3.ItemTokenSize = [12,22];
% le3.Position = le3.Position + [0.025 0.045 0 0];


% subplot4: SessionGroup, RelT distribution, legend: Pre/Post
%ha4 = axes;
%set(ha4, 'units', 'centimeters', 'position', [xs2(1) ys(2) size1], 'nextplot', 'add','tickDir', 'out',...
%    'xlim',[xedges.edges_RelT(1),xedges.edges_RelT(end)],'xtick',[xedges.edges_RelT(1):0.25:xedges.edges_RelT(end)],...
%    'xticklabel',cellstr(string((xedges.edges_RelT(1):0.25:xedges.edges_RelT(end))*1000)),'fontsize',7,'fontname','Arial',...
%    'ylim',[0 0.17],'ytick', [0:0.05:0.15], 'yticklabel', {'0', '5', '10', '15'}, 'ticklength', [0.02 0.025]);
%fill([0,0.6,0.6,0],[0,0,1,1],cGreen,'EdgeColor','none','FaceAlpha',0.25);
% edges_RT = 0:0.05:0.6;
%f41 = shadedErrorBar(xedges.RelT,GAS.mean_RelTdist(GAS.Group=="Pre",:),...
%    GAS.sem_RelTdist(GAS.Group=="Pre",:),...
%    'lineProps',{'-','lineWidth',1.5,'color','k'},'patchSaturation',0.2);
%f42 = shadedErrorBar(xedges.RelT,GAS.mean_RelTdist(GAS.Group=="Post",:),...
%    GAS.sem_RelTdist(GAS.Group=="Post",:),...
%    'lineProps',{'-','lineWidth',1.5,'color',cOrange},'patchSaturation',0.2);
%xlabel('Release time (ms)','Fontsize',7,'FontName','Arial');
%ylabel('Probability (%)','Fontsize',7,'FontName','Arial');

%le4 = legend({'Short','Medium','Long'},...
%    'Fontsize',7,'fontname','Arial','units','centimeters','Position',[11.95,9.3,1,1]);
%legend('boxoff');
%le4.ItemTokenSize = [12,22];
%le4.Position = le4.Position + [0.025 0.045 0 0];

%xlabel('Reaction time (ms)','Fontsize',8);
%ylabel('Probability','Fontsize',8);
%title('Session','Fontsize',8);

% subplot5: TrialGroup, RelT distribution, legend: Pre/Post
ha5 = axes;
set(ha5, 'units', 'centimeters', 'position', [xs2(1) ys(3) size1], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[xedges.edges_RelT(1),xedges.edges_RelT(end)],'xtick',[xedges.edges_RelT(1):0.25:xedges.edges_RelT(end)],...
    'xticklabel',cellstr(string((xedges.edges_RelT(1):0.25:xedges.edges_RelT(end))*1000)),'fontsize',7,'fontname','Arial',...
    'ylim',[0 0.17],'ytick', [0:0.05:0.15], 'yticklabel', {'0', '5', '10', '15'}, 'ticklength', [0.02 0.025]);
fill([0,0.6,0.6,0],[0,0,1,1],cGreen,'EdgeColor','none','FaceAlpha',0.25);
f41 = shadedErrorBar(xedges.RelT,GAS_trial.mean_RelTdist(GAS_trial.Group=="Pre",:),...
    GAS_trial.sem_RelTdist(GAS_trial.Group=="Pre",:),...
    'lineProps',{'-','lineWidth',1.5,'color','k'},'patchSaturation',0.2);
f42 = shadedErrorBar(xedges.RelT,GAS_trial.mean_RelTdist(GAS_trial.Group=="Post",:),...
    GAS_trial.sem_RelTdist(GAS_trial.Group=="Post",:),...
    'lineProps',{'-','lineWidth',1.5,'color',cOrange},'patchSaturation',0.2);
xlabel('Release time (ms)','Fontsize',7,'FontName','Arial');
ylabel('Probability (%)','Fontsize',7,'FontName','Arial');
% edges_RT = 0:0.05:0.6;
%plot(0.025:0.05:0.575,GAS_trial.mean_RTdist_S(GAS_trial.Group=="Pre",:),':','lineWidth',2,'color',cGray);
%plot(0.025:0.05:0.575,GAS_trial.mean_RTdist_M(GAS_trial.Group=="Pre",:),'--','lineWidth',1.7,'color',cGray);
%plot(0.025:0.05:0.575,GAS_trial.mean_RTdist_L(GAS_trial.Group=="Pre",:),'-','lineWidth',1.5,'color',cGray);
%plot(0.025:0.05:0.575,GAS_trial.mean_RTdist_S(GAS_trial.Group=="Post",:),':','lineWidth',2,'color',cOrange);
%plot(0.025:0.05:0.575,GAS_trial.mean_RTdist_M(GAS_trial.Group=="Post",:),'--','lineWidth',1.7,'color',cOrange);
%plot(0.025:0.05:0.575,GAS_trial.mean_RTdist_L(GAS_trial.Group=="Post",:),'-','lineWidth',1.5,'color',cOrange);

% le5 = legend({'Short','Medium','Long'},...
%     'Fontsize',6,'units','centimeters','Position',[16,5.7,1,1]);%[10.3,6.8,1,1]
% legend('boxoff');
% le5.ItemTokenSize = [12,22];
% le5.Position = le5.Position + [0.025 0.045 0 0];

%xlabel('Reaction time (ms)','Fontsize',8);
%ylabel('Probability','Fontsize',8);
%title('Trial','Fontsize',8);

% subplot6: x:RT-Pre, y:RT-Post, marker:Session/Trial, each point a subjects
ha6 = axes;
set(ha6, 'units', 'centimeters', 'position', [xs2(2) ys(1) size5], 'nextplot', 'add','tickDir', 'out',...
   'fontsize',7,'fontname','Arial','xlim',[0.14 0.46],'ylim',[0.14 0.46],'xtick',[0.15,0.25,0.35,0.45],'xticklabel',{'150','250','350','450'},...
   'ytick',[0.15,0.25,0.35,0.45],'yticklabel',{'150','250','350','450'}, 'ticklength', [0.02 0.025]); %

plot([0,0.6],[0,0.6],':k','lineWidth',0.6)
%h61 = errorbar(EES.TBT.RT(EES.TBT.Group=="Pre"),EES.TBT.RT(EES.TBT.Group=="Post"),...
%    EES.TBT.RT_SE(EES.TBT.Group=="Post"),EES.TBT.RT_SE(EES.TBT.Group=="Post"),...
%    EES.TBT.RT_SE(EES.TBT.Group=="Pre"),EES.TBT.RT_SE(EES.TBT.Group=="Pre"),...
%    '.','MarkerSize',4,'MarkerEdgeColor',cBlue2,'color',cBlue2,'lineWidth',1,'CapSize',1);
h62 = errorbar(EES_trial.TBT.RT(EES_trial.TBT.Group=="Pre"),EES_trial.TBT.RT(EES_trial.TBT.Group=="Post"),...
    EES_trial.TBT.RT_SE(EES_trial.TBT.Group=="Post"),EES_trial.TBT.RT_SE(EES_trial.TBT.Group=="Post"),...
    EES_trial.TBT.RT_SE(EES_trial.TBT.Group=="Pre"),EES_trial.TBT.RT_SE(EES_trial.TBT.Group=="Pre"),...
    '.','MarkerSize',4,'MarkerEdgeColor',cBlue2,'color',cBlue2,'lineWidth',1,'CapSize',1);


%le6 = legend([h61,h62],{'Session','Trial'},'Fontsize',7,'FontName','Arial','units','centimeters','Position',[xs2(2)+3.5*0.7+0.25,ys(1)+0.7,1,1]);
%legend('boxoff');
%le6.ItemTokenSize = [12,22];
%le6.Position = le6.Position + [0.025 0.045 0 0];

xlabel('Pre RT (ms)','Fontsize',7,'FontName','Arial');
ylabel('Post RT (ms)','Fontsize',7,'FontName','Arial');

% subplot1: x:RelT-Pre, y:RelT-Post, marker:Session/Trial, each point a subjects
ha1 = axes;
set(ha1, 'units', 'centimeters', 'position', [xs2(1) ys(1) size5], 'nextplot', 'add','tickDir', 'out',...
   'fontsize',7,'fontname','Arial','xlim',[0.14 0.56],'ylim',[0.14 0.56],'xtick',[0.15,0.25,0.35,0.45,0.55],'xticklabel',{'150','250','350','450','550'},...
   'ytick',[0.15,0.25,0.35,0.45,0.55],'yticklabel',{'150','250','350','450','550'}, 'ticklength', [0.02 0.025],'xticklabelRotation',0); %

plot([0,1.5],[0,1.5],':k','lineWidth',0.6)
%errorbar(EES.TBT.RelT(EES.TBT.Group=="Pre"),EES.TBT.RelT(EES.TBT.Group=="Post"),...
%    EES.TBT.RelT_SE(EES.TBT.Group=="Post"),EES.TBT.RelT_SE(EES.TBT.Group=="Post"),...
%    EES.TBT.RelT_SE(EES.TBT.Group=="Pre"),EES.TBT.RelT_SE(EES.TBT.Group=="Pre"),...
%    '.','MarkerSize',4,'MarkerEdgeColor',cBlue,'color',cBlue,'lineWidth',1,'CapSize',1);
errorbar(EES_trial.TBT.RelT(EES_trial.TBT.Group=="Pre"),EES_trial.TBT.RelT(EES_trial.TBT.Group=="Post"),...
    EES_trial.TBT.RelT_SE(EES_trial.TBT.Group=="Post"),EES_trial.TBT.RelT_SE(EES_trial.TBT.Group=="Post"),...
    EES_trial.TBT.RelT_SE(EES_trial.TBT.Group=="Pre"),EES_trial.TBT.RelT_SE(EES_trial.TBT.Group=="Pre"),...
    '.','MarkerSize',4,'MarkerEdgeColor',cBlue,'color',cBlue,'lineWidth',1,'CapSize',1);

xlabel('Pre RelT (ms)','Fontsize',7,'FontName','Arial');
ylabel('Post RelT (ms)','Fontsize',7,'FontName','Arial');

% subplot10: x:Cor-Pre, y:Cor-Post, marker:Session/Trial, each point a subjects
ha10 = axes;
set(ha10, 'units', 'centimeters', 'position', [xs2(1) ys(4) size5], 'nextplot', 'add','tickDir', 'out',...
   'fontsize',7,'fontname','Arial','xlim',[0 1],'ylim',[0 1],'xtick',[0,0.25,0.50,0.75,1.00],'xticklabel',{'0','25','50','75','100'},'xticklabelRotation',0,...
   'ytick',[0,0.25,0.50,0.75,1.00],'yticklabel',{'0','25','50','75','100'}, 'ticklength', [0.02 0.025]); %

plot([0,1],[0,1],':k','lineWidth',0.6)
%h101 = plot(EES.TBT.Cor(EES.TBT.Group=="Pre"),EES.TBT.Cor(EES.TBT.Group=="Post"),...
%    'o','MarkerSize',4,'MarkerEdgeColor',cGreen);
h102 = plot(EES_trial.TBT.Cor(EES_trial.TBT.Group=="Pre"),EES_trial.TBT.Cor(EES_trial.TBT.Group=="Post"),...
    'o','MarkerSize',4,'MarkerEdgeColor',cGreen);
xlabel('Pre Cor (%)','Fontsize',7,'FontName','Arial');
ylabel('Post Cor (%)','Fontsize',7,'FontName','Arial');

% subplot11: x:Pre-Pre, y:Pre-Post, marker:Session/Trial, each point a subjects
ha11 = axes;
set(ha11, 'units', 'centimeters', 'position', [xs2(2) ys(4) size5], 'nextplot', 'add','tickDir', 'out',...
   'fontsize',7,'fontname','Arial','xlim',[0 1],'ylim',[0 1],'xtick',[0,0.25,0.50,0.75,1.00],'xticklabel',{'0','25','50','75','100'},'xticklabelRotation',0,...
   'ytick',[0,0.25,0.50,0.75,1.00],'yticklabel',{'0','25','50','75','100'}, 'ticklength', [0.02 0.025]); %

plot([0,1],[0,1],':k','lineWidth',0.6)
%h111 = plot(EES.TBT.Pre(EES.TBT.Group=="Pre"),EES.TBT.Pre(EES.TBT.Group=="Post"),...
%    'o','MarkerSize',4,'MarkerEdgeColor',cRed);
h112 = plot(EES_trial.TBT.Pre(EES_trial.TBT.Group=="Pre"),EES_trial.TBT.Pre(EES_trial.TBT.Group=="Post"),...
    'o','MarkerSize',4,'MarkerEdgeColor',cRed);
xlabel('Pre Pre (%)','Fontsize',7,'FontName','Arial');
ylabel('Post Pre (%)','Fontsize',7,'FontName','Arial');

% subplot12: x:Late-Pre, y:Late-Post, marker:Session/Trial, each point a subjects
ha12 = axes;
set(ha12, 'units', 'centimeters', 'position', [xs2(3) ys(4) size5], 'nextplot', 'add','tickDir', 'out',...
   'fontsize',7,'fontname','Arial','xlim',[0 1],'ylim',[0 1],'xtick',[0,0.25,0.50,0.75,1.00],'xticklabel',{'0','25','50','75','100'},'xticklabelRotation',0,...
   'ytick',[0,0.25,0.50,0.75,1.00],'yticklabel',{'0','25','50','75','100'}, 'ticklength', [0.02 0.025]); %

plot([0,1],[0,1],':k','lineWidth',0.6)
%h121 = plot(EES.TBT.Late(EES.TBT.Group=="Pre"),EES.TBT.Late(EES.TBT.Group=="Post"),...
%    'o','MarkerSize',4,'MarkerEdgeColor',cGray);
h122 = plot(EES_trial.TBT.Late(EES_trial.TBT.Group=="Pre"),EES_trial.TBT.Late(EES_trial.TBT.Group=="Post"),...
    'o','MarkerSize',4,'MarkerEdgeColor',cGray);
xlabel('Pre Late (%)','Fontsize',7,'FontName','Arial');
ylabel('Post Late (%)','Fontsize',7,'FontName','Arial');

%le12 = legend([h121,h122],{'Session','Trial'},'Fontsize',7,'FontName','Arial','units','centimeters','Position',[xs2(3)+3.5*0.7+0.25,ys(4)+0.7,1,1]);
%legend('boxoff');
%le12.ItemTokenSize = [12,22];
%le12.Position = le12.Position + [0.025 0.045 0 0];

% subplot7: SessionGroup, RelT distribution CDF, legend: Pre/Post
%ha7 = axes;
%set(ha7, 'units', 'centimeters', 'position', [xs2(2) ys(2) size1], 'nextplot', 'add','tickDir', 'out',...
%    'xlim',[xedges.edges_RelT(1),xedges.edges_RelT(end)],'xtick',[xedges.edges_RelT(1):0.25:xedges.edges_RelT(end)],...
%    'xticklabel',cellstr(string((xedges.edges_RelT(1):0.25:xedges.edges_RelT(end))*1000)),'fontsize',7,'fontname','Arial',...
%    'ylim',[0 1],'ytick', [0:0.25:1], 'yticklabel', {'0', '25', '50', '75','100'}, 'ticklength', [0.02 0.025]);
%fill([0,0.6,0.6,0],[0,0,1,1],cGreen,'EdgeColor','none','FaceAlpha',0.25);
%f41 = shadedErrorBar(xedges.RelT,GAS.mean_RelTdist_CDF(GAS.Group=="Pre",:),...
%    GAS.sem_RelTdist_CDF(GAS.Group=="Pre",:),...
%    'lineProps',{'-','lineWidth',1.5,'color','k'},'patchSaturation',0.2);
%f42 = shadedErrorBar(xedges.RelT,GAS.mean_RelTdist_CDF(GAS.Group=="Post",:),...
%    GAS.sem_RelTdist_CDF(GAS.Group=="Post",:),...
%    'lineProps',{'-','lineWidth',1.5,'color',cOrange},'patchSaturation',0.2);
%xlabel('Release time (ms)','Fontsize',7,'FontName','Arial');
%ylabel('CDF (%)','Fontsize',7,'FontName','Arial');
% edges_RT = 0:0.05:0.6;
%plot(0.025:0.05:2.475,GAS.mean_HTdist_S(GAS.Group=="Pre",:),':','lineWidth',2,'color',cGray);
%plot(0.025:0.05:2.475,GAS.mean_HTdist_M(GAS.Group=="Pre",:),'--','lineWidth',1.7,'color',cGray);
%h71 = plot(0.025:0.05:2.475,GAS.mean_HTdist_L(GAS.Group=="Pre",:),'-','lineWidth',1.5,'color',cGray);
%plot(0.025:0.05:2.475,GAS.mean_HTdist_S(GAS.Group=="Post",:),':','lineWidth',2,'color',cOrange);
%plot(0.025:0.05:2.475,GAS.mean_HTdist_M(GAS.Group=="Post",:),'--','lineWidth',1.7,'color',cOrange);
%h72 = plot(0.025:0.05:2.475,GAS.mean_HTdist_L(GAS.Group=="Post",:),'-','lineWidth',1.5,'color',cOrange);

%le7 = legend([h71 h72],{'Pre','Post'},'fontsize',6,'units','centimeter','Position',[11.8,10.6,1,1]); % [4.6,6.9,1,1]
%legend('boxoff');
%le7.ItemTokenSize = [12,22];
%le7.Position = le7.Position + [0.025 0.045 0 0];

%xlabel('Press duration (ms)','Fontsize',8);
%ylabel('Probability','Fontsize',8);
%title('Session','Fontsize',8);

% subplot8: TrialGroup, RelT distribution CDF, legend: Pre/Post
ha8 = axes;
set(ha8, 'units', 'centimeters', 'position', [xs2(2) ys(3) size1], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[xedges.edges_RelT(1),xedges.edges_RelT(end)],'xtick',[xedges.edges_RelT(1):0.25:xedges.edges_RelT(end)],...
    'xticklabel',cellstr(string((xedges.edges_RelT(1):0.25:xedges.edges_RelT(end))*1000)),'fontsize',7,'fontname','Arial',...
    'ylim',[0 1],'ytick', [0:0.25:1], 'yticklabel', {'0', '25', '50', '75','100'}, 'ticklength', [0.02 0.025]);
fill([0,0.6,0.6,0],[0,0,1,1],cGreen,'EdgeColor','none','FaceAlpha',0.25);
f41 = shadedErrorBar(xedges.RelT,GAS_trial.mean_RelTdist_CDF(GAS_trial.Group=="Pre",:),...
    GAS_trial.sem_RelTdist_CDF(GAS_trial.Group=="Pre",:),...
    'lineProps',{'-','lineWidth',1.5,'color','k'},'patchSaturation',0.2);
f42 = shadedErrorBar(xedges.RelT,GAS_trial.mean_RelTdist_CDF(GAS_trial.Group=="Post",:),...
    GAS_trial.sem_RelTdist_CDF(GAS_trial.Group=="Post",:),...
    'lineProps',{'-','lineWidth',1.5,'color',cOrange},'patchSaturation',0.2);
xlabel('Release time (ms)','Fontsize',7,'FontName','Arial');
ylabel('CDF (%)','Fontsize',7,'FontName','Arial');
% edges_RT = 0:0.05:0.6;
%plot(0.025:0.05:2.475,GAS_trial.mean_HTdist_S(GAS_trial.Group=="Pre",:),':','lineWidth',2,'color',cGray);
%plot(0.025:0.05:2.475,GAS_trial.mean_HTdist_M(GAS_trial.Group=="Pre",:),'--','lineWidth',1.7,'color',cGray);
%plot(0.025:0.05:2.475,GAS_trial.mean_HTdist_L(GAS_trial.Group=="Pre",:),'-','lineWidth',1.5,'color',cGray);
%plot(0.025:0.05:2.475,GAS_trial.mean_HTdist_S(GAS_trial.Group=="Post",:),':','lineWidth',2,'color',cOrange);
%plot(0.025:0.05:2.475,GAS_trial.mean_HTdist_M(GAS_trial.Group=="Post",:),'--','lineWidth',1.7,'color',cOrange);
%plot(0.025:0.05:2.475,GAS_trial.mean_HTdist_L(GAS_trial.Group=="Post",:),'-','lineWidth',1.5,'color',cOrange);

%xlabel('Press duration (ms)','Fontsize',8);
%ylabel('Probability','Fontsize',8);
%title('Trial','Fontsize',8);

%text(3.7,0.1,{"Trials_{pre}  = "+num2str(nSession_pre,'%.2f')+"sessions","Trials_{post} = "+num2str(nSession_post,'%.2f')+"sessions"},'FontSize',8);

% PLOT SessionGroup, 3FPs-RelT, Pre/Post
% S
%RelTMedPreMean(1) = GAS.mean_RelT_S(GAS.Group=="Pre");
%RelTMedPreSEM(1) = GAS.sem_RelT_S(GAS.Group=="Pre");
%RelTMedPostMean(1) = GAS.mean_RelT_S(GAS.Group=="Post");
%RelTMedPostSEM(1) = GAS.sem_RelT_S(GAS.Group=="Post");
% M
%RelTMedPreMean(2) = GAS.mean_RelT_M(GAS.Group=="Pre");
%RelTMedPreSEM(2) = GAS.sem_RelT_M(GAS.Group=="Pre");
%RelTMedPostMean(2) = GAS.mean_RelT_M(GAS.Group=="Post");
%RelTMedPostSEM(2) = GAS.sem_RelT_M(GAS.Group=="Post");
% L
%RelTMedPreMean(3) = GAS.mean_RelT_L(GAS.Group=="Pre");
%RelTMedPreSEM(3) = GAS.sem_RelT_L(GAS.Group=="Pre");
%RelTMedPostMean(3) = GAS.mean_RelT_L(GAS.Group=="Post");
%RelTMedPostSEM(3) = GAS.sem_RelT_L(GAS.Group=="Post");
% plot
%ha23 = axes;
%set(ha23, 'units', 'centimeters', 'position', [xs2(3) ys(2) size1], 'nextplot', 'add','tickDir', 'out',...
%    'xlim',[fplist(1)-0.25,fplist(end)+0.25],'xtick',fplist,'xticklabel',cellstr(string(fplist.*1000)),'fontsize',7,'fontname','Arial', ...
%    'ylim',[0.2 0.6],'ytick',0.2:0.1:0.5,'yticklabels',{'200','300','400','500'},'ticklength', [0.02 0.025]);
%plot(fplist, RelTMedPreMean, 'o-', 'linewidth', 1, 'color', 'k', 'markerfacecolor', 'k', 'markeredgecolor','w', 'markersize', 5)
%line([fplist; fplist], [RelTMedPreMean-RelTMedPreSEM; RelTMedPreMean+RelTMedPreSEM], ...
%    'color','k', 'linewidth', 1)

%plot(fplist, RelTMedPostMean, 'o-', 'linewidth', 1, 'color', cOrange, 'markerfacecolor', cOrange, 'markeredgecolor','w', 'markersize', 5)
%line([fplist; fplist], [RelTMedPostMean-RelTMedPostSEM; RelTMedPostMean+RelTMedPostSEM], ...
%    'color',cOrange, 'linewidth', 1)

%xlabel('Foreperiod (ms)','Fontsize',7,'FontName','Arial')
%ylabel('Release time (ms)','Fontsize',7,'FontName','Arial')

%xlm = xlim; xsep = (xlm(2)-xlm(1))./15; xtext = xlm(2)+xsep;
%ylm = ylim; ysep = (ylm(2)-ylm(1))./5; ytext = mean(ylm)+ysep;
%text(xtext,ytext,"Session",...
%    'HorizontalAlignment','left','VerticalAlignment','middle','rotation',-90,...
%    'fontsize',8,'FontName','Arial','fontweight','bold');

% PLOT TrialGroup, 3FPs-RelT, Pre/Post
% S
RelTMedPreMean(1) = GAS_trial.mean_RelT_S(GAS_trial.Group=="Pre");
RelTMedPreSEM(1) = GAS_trial.sem_RelT_S(GAS_trial.Group=="Pre");
RelTMedPostMean(1) = GAS_trial.mean_RelT_S(GAS_trial.Group=="Post");
RelTMedPostSEM(1) = GAS_trial.sem_RelT_S(GAS_trial.Group=="Post");
% M
RelTMedPreMean(2) = GAS_trial.mean_RelT_M(GAS_trial.Group=="Pre");
RelTMedPreSEM(2) = GAS_trial.sem_RelT_M(GAS_trial.Group=="Pre");
RelTMedPostMean(2) = GAS_trial.mean_RelT_M(GAS_trial.Group=="Post");
RelTMedPostSEM(2) = GAS_trial.sem_RelT_M(GAS_trial.Group=="Post");
% L
RelTMedPreMean(3) = GAS_trial.mean_RelT_L(GAS_trial.Group=="Pre");
RelTMedPreSEM(3) = GAS_trial.sem_RelT_L(GAS_trial.Group=="Pre");
RelTMedPostMean(3) = GAS_trial.mean_RelT_L(GAS_trial.Group=="Post");
RelTMedPostSEM(3) = GAS_trial.sem_RelT_L(GAS_trial.Group=="Post");
% plot
ha33 = axes;
set(ha33, 'units', 'centimeters', 'position', [xs2(3) ys(3) size1], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[fplist(1)-0.25,fplist(end)+0.25],'xtick',fplist,'xticklabel',cellstr(string(fplist.*1000)),'fontsize',7,'fontname','Arial', ...
    'ylim',[0.2 0.6],'ytick',0.2:0.1:0.5,'yticklabels',{'200','300','400','500'},'ticklength', [0.02 0.025]);
lpre = plot(fplist, RelTMedPreMean, 'o-', 'linewidth', 1, 'color', 'k', 'markerfacecolor', 'k', 'markeredgecolor','w', 'markersize', 5);
line([fplist; fplist], [RelTMedPreMean-RelTMedPreSEM; RelTMedPreMean+RelTMedPreSEM], ...
    'color','k', 'linewidth', 1)
lpost = plot(fplist, RelTMedPostMean, 'o-', 'linewidth', 1, 'color', cOrange, 'markerfacecolor', cOrange, 'markeredgecolor','w', 'markersize', 5);
line([fplist; fplist], [RelTMedPostMean-RelTMedPostSEM; RelTMedPostMean+RelTMedPostSEM], ...
    'color',cOrange, 'linewidth', 1)
xlabel('Foreperiod (ms)','Fontsize',7,'FontName','Arial')
ylabel('Release time (ms)','Fontsize',7,'FontName','Arial')

%xlm = xlim; xsep = (xlm(2)-xlm(1))./15; xtext = xlm(2)+xsep;
%ylm = ylim; ysep = (ylm(2)-ylm(1))./5; ytext = mean(ylm)+ysep;
%text(xtext,ytext,"Trial",...
%    'HorizontalAlignment','left','VerticalAlignment','middle','rotation',-90,...
%    'fontsize',8,'FontName','Arial','fontweight','bold');

le9 = legend([lpre,lpost],{'Pre','Post'},...
    'Fontsize',7,'FontName','Arial','units','centimeters',...
    'Position',[xs2(3)+size1(1)+0.6,ys(3)-0.2,1,1]); % [10.7,8.9,1,1]
legend('boxoff');
le9.ItemTokenSize = [12,22];
le9.Position = le9.Position + [0.025 0.045 0 0];
%% Save
savename = fullfile(pwd,'MiceLesion_L2V');
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
        [HT_rmvd, HT_idxrmv] = rmoutliers_custome(tdata.PressDur);
        t.HT = mean(HT_rmvd,'omitnan');
        %t.RT = mean(tdata(tdata.Type==1,:).RT);
        %RT = rmoutliers(tdata(tdata.Type==1,:).RT,'mean');
        %t.RT = mean(RT(RT>=0.1),'omitnan');
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
