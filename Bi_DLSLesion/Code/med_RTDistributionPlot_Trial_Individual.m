function med_RTDistributionPlot_Trial_Individual(btAll2d,trialNearBound_LP)
%%
altMethod = {'mean','median','geomean'};
cenMethod = altMethod{1};
errMethod = 'sem';
%% Data processing
global lesionBoundary;

EES = struct;
[EES.SBS,TBT] = packData(btAll2d,cenMethod);
TBT = TBT(TBT.Type~=0,:);

[TBT_fp05_plot, trial_idx05, session_idx05, missing_session_idx05] = plotRangeTrial(0.5,TBT);
[TBT_fp10_plot, trial_idx10, session_idx10, missing_session_idx10] = plotRangeTrial(1.0,TBT);
[TBT_fp15_plot, trial_idx15, session_idx15, missing_session_idx15] = plotRangeTrial(1.5,TBT);
[TBT_3fp_plot, trial_idx3, session_idx3, missing_session_idx3] = plotRangeTrial([0.5,1.0,1.5],TBT);

TBT_fp05_use_pre = TBT_fp05_plot(find(TBT_fp05_plot.Date<=lesionBoundary,trialNearBound_LP,'last'),:);
TBT_fp05_use_post = TBT_fp05_plot(find(TBT_fp05_plot.Date>lesionBoundary,trialNearBound_LP,'first'),:);
TBT_fp05_use = [TBT_fp05_use_pre;TBT_fp05_use_post];
TBT_fp05_use_RT = [TBT_fp05_use.RT;NaN(trialNearBound_LP*2-length(TBT_fp05_use.RT),1)];
median05 = median(reshape(TBT_fp05_use_RT,10,[]),'omitnan');
quantile05_low = quantile(reshape(TBT_fp05_use_RT,10,[]), 0.25,1);
errbar05_min = median05-quantile05_low;
quantile05_high = quantile(reshape(TBT_fp05_use_RT,10,[]), 0.75,1);
errbar05_max = quantile05_high-median05;
TBT_fp10_use_pre = TBT_fp10_plot(find(TBT_fp10_plot.Date<=lesionBoundary,trialNearBound_LP,'last'),:);
TBT_fp10_use_post = TBT_fp10_plot(find(TBT_fp10_plot.Date>lesionBoundary,trialNearBound_LP,'first'),:);
TBT_fp10_use = [TBT_fp10_use_pre;TBT_fp10_use_post];
TBT_fp10_use_RT = [TBT_fp10_use.RT;NaN(trialNearBound_LP*2-length(TBT_fp10_use.RT),1)];
median10 = median(reshape(TBT_fp10_use_RT,10,[]),'omitnan');
quantile10_low = quantile(reshape(TBT_fp10_use_RT,10,[]), 0.25,1);
errbar10_min = median10-quantile10_low;
quantile10_high = quantile(reshape(TBT_fp10_use_RT,10,[]), 0.75,1);
errbar10_max = quantile10_high-median10;
TBT_fp15_use_pre = TBT_fp15_plot(find(TBT_fp15_plot.Date<=lesionBoundary,trialNearBound_LP,'last'),:);
TBT_fp15_use_post = TBT_fp15_plot(find(TBT_fp15_plot.Date>lesionBoundary,trialNearBound_LP,'first'),:);
TBT_fp15_use = [TBT_fp15_use_pre;TBT_fp15_use_post];
TBT_fp15_use_RT = [TBT_fp15_use.RT;NaN(trialNearBound_LP*2-length(TBT_fp15_use.RT),1)];
median15 = median(reshape(TBT_fp15_use_RT,10,[]),'omitnan');
quantile15_low = quantile(reshape(TBT_fp15_use_RT,10,[]), 0.25,1);
errbar15_min = median15-quantile15_low;
quantile15_high = quantile(reshape(TBT_fp15_use_RT,10,[]), 0.75,1);
errbar15_max = quantile15_high-median15;


save('VarsToPlot_RTD.mat','TBT','TBT_fp05_plot','TBT_fp10_plot','TBT_fp15_plot','TBT_3fp_plot');
%% Plot
load('VarsToPlot_RTD.mat');

cDarkGray = [0.2,0.2,0.2];
cGreen = [0.4660 0.6740 0.1880];
cRed = [0.6350 0.0780 0.1840];
cYellow = [0.9290 0.6940 0.1250];
cBlue = [0,0.6902,0.9412];
cOrange = [0.929,0.49,0.192];
cGray = [0.4 0.4 0.4];
colorlist = [cDarkGray;cGray;cGreen;cRed;cYellow;cBlue];
color_FP = [cGray;mean([cOrange;cGray]);cOrange];


hf = figure(44); clf(hf,'reset');
set(hf, 'name', 'RT distribution_trials', 'units', 'centimeters', 'position', [35 3 18.542 12.5],...
    'PaperPositionMode', 'auto');



ax1 = subplot(5,2,[1,2]);
set(ax1, 'nextplot', 'add','tickDir', 'out', 'ylim',[0 0.7],...
    'fontsize',6);
sc1 = scatter(trial_idx3(TBT_3fp_plot.FP==0.5,:),TBT_3fp_plot.RT(TBT_3fp_plot.FP==0.5,:),3,color_FP(1,:),'filled');
sc2 = scatter(trial_idx3(TBT_3fp_plot.FP==1.0,:),TBT_3fp_plot.RT(TBT_3fp_plot.FP==1.0,:),3,color_FP(2,:),'filled');
sc3 = scatter(trial_idx3(TBT_3fp_plot.FP==1.5,:),TBT_3fp_plot.RT(TBT_3fp_plot.FP==1.5,:),3,color_FP(3,:),'filled');
xline(0,'-.');
plot([session_idx3,session_idx3],[0,0.1],'-','color',cBlue);
if isempty(missing_session_idx3) ~= 1
plot([missing_session_idx3,missing_session_idx3],[0,0.1],'-','color','r');
end
ylabel('RT(s)','Fontsize',8);
title(TBT.Subject(1) + ": RT Distribution",'Fontsize',8);

le = legend([sc1,sc2,sc3],{'FP=0.5s','FP=1.0s','FP=1.5s'},...
    'Fontsize',6,'Position',[0.89 0.79 0.04 0.05]);
legend('boxoff');
le.ItemTokenSize = [12,22];
le.Position = le.Position + [0.025 0.045 0 0];

ax2 = subplot(5,2,[3,4]);
set(ax2, 'nextplot', 'add','tickDir', 'out', 'ylim',[0 0.7],...
    'fontsize',6);
scatter(trial_idx05,TBT_fp05_plot.RT,3,color_FP(1,:),'filled');
xline(0,'-.');
plot([session_idx05,session_idx05],[0,0.1],'-','color',cBlue);
if isempty(missing_session_idx05) ~= 1
plot([missing_session_idx05,missing_session_idx05],[0,0.1],'-','color','r');
end
ylabel('RT_{FP0.5}(s)','Fontsize',8);

ax3 = subplot(5,2,[5,6]);
set(ax3, 'nextplot', 'add','tickDir', 'out', 'ylim',[0 0.7],...
    'fontsize',6);
scatter(trial_idx10,TBT_fp10_plot.RT,3,color_FP(2,:),'filled');
xline(0,'-.');
plot([session_idx10,session_idx10],[0,0.1],'-','color',cBlue);
if isempty(missing_session_idx10) ~= 1
plot([missing_session_idx10,missing_session_idx10],[0,0.1],'-','color','r');
end
ylabel('RT_{FP1.0}(s)','Fontsize',8);

ax4 = subplot(5,2,[7,8]);
set(ax4, 'nextplot', 'add','tickDir', 'out', 'ylim',[0 0.7],...
    'fontsize',6);
scatter(trial_idx15,TBT_fp15_plot.RT,3,color_FP(3,:),'filled');
xline(0,'-.');
plot([session_idx15,session_idx15],[0,0.1],'-','color',cBlue);
if isempty(missing_session_idx15) ~= 1
plot([missing_session_idx15,missing_session_idx15],[0,0.1],'-','color','r');
end
xlabel('# of trial','Fontsize',8);
ylabel('RT_{FP1.5}(s)','Fontsize',8);

ax5 = subplot(5,2,9);
set(ax5, 'nextplot', 'add','tickDir', 'out', 'ylim',[0 0.7],...
    'fontsize',6);
scatter([-50:-1,1:50]',TBT_fp05_use_RT,3,color_FP(1,:),'filled','MarkerFaceAlpha',0.5);
scatter([-50:-1,1:50]',TBT_fp10_use_RT,3,color_FP(2,:),'filled','MarkerFaceAlpha',0.5);
scatter([-50:-1,1:50]',TBT_fp15_use_RT,3,color_FP(3,:),'filled','MarkerFaceAlpha',0.5);
errorbar([-50:10:-10,10:10:50]-2,median05,errbar05_min,errbar05_max,'color',color_FP(1,:),'lineWidth',1,'CapSize',3);
errorbar([-50:10:-10,10:10:50],median10,errbar10_min,errbar10_max,'color',color_FP(2,:),'lineWidth',1,'CapSize',3);
errorbar([-50:10:-10,10:10:50]+2,median15,errbar15_min,errbar15_max,'color',color_FP(3,:),'lineWidth',1,'CapSize',3);
xline(0,'-.');
xlabel('# of trial','Fontsize',8);
ylabel('RT(s)','Fontsize',8);
%% Save
savename = fullfile(pwd,"RTDistribution_"+TBT.Subject(1));
saveas(hf,savename,'fig');
print(hf,'-dpng',savename);
print(hf,'-depsc2',savename);
end
%% Functions
function [SBS,TBT] = packData(btAll2d,cenMethod)
SBS = table; % session by session data
TBT = table; % trial by trial data
for i=1:size(btAll2d,1)
    for j=1:size(btAll2d,2)
        T = btAll2d{i,j};
        SBS = [SBS;estSBS(T,j,cenMethod)];
        
        nrow = size(T,1);
        if nrow>1
            %tempT = addvars(T,repelem(grp(i),nrow)','Before','Date','NewVariableNames','Group');
            tempTT = addvars(T,repelem(j,nrow)','After','Date','NewVariableNames','Session');
            TBT = [TBT;tempTT];
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
        t.HT = mean(tdata.PressDur);
        t.RT = mean(tdata(tdata.Type==1,:).RT);
    case 'median'
        t.HT = median(tdata.PressDur);
        t.RT = median(tdata(tdata.Type==1,:).RT);
    case 'geomean'
        t.HT = geomean(tdata.PressDur);
        t.RT = geomean(tdata(tdata.Type==1 & tdata.RT>0,:).RT);
end
outT = [outT;struct2table(t)];

end

function outT = estTBT_3FPs(TBT)
fplist = [0.5,1.0,1.5];

outT = table;
sbjlist = unique(TBT.Subject);

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

    t.RT = mean(tdata.RT(idxCor),'omitnan');
    t.RT_sem = std(tdata.RT(idxCor),'omitnan')/sqrt(sum(~isnan(tdata.RT(idxCor))));
    t.RT_S = mean(tdata.RT(idxCor&idxFPS));
    t.RT_M = mean(tdata.RT(idxCor&idxFPM));
    t.RT_L = mean(tdata.RT(idxCor&idxFPL));
%         t.RT_3FPs = [mean(tdata.RT(idxCor&idxFPS)),mean(tdata.RT(idxCor&idxFPM)),mean(tdata.RT(idxCor&idxFPL))];
    edges_RT = 0:0.05:0.6;
    t.RTdist_S = smooth(histcounts(tdata.RT(idxCor&idxFPS),edges_RT,'Normalization','probability'),3)';
    t.RTdist_M = smooth(histcounts(tdata.RT(idxCor&idxFPM),edges_RT,'Normalization','probability'),3)';
    t.RTdist_L = smooth(histcounts(tdata.RT(idxCor&idxFPL),edges_RT,'Normalization','probability'),3)';

    edges_HT = 0:0.05:2.5;
    t.HTdist = histcounts(tdata.PressDur,edges_HT,'Normalization','probability');
    t.HTdist_S = smooth(histcounts(tdata.PressDur(idxFPS),edges_HT,'Normalization','probability'),3)';
    t.HTdist_M = smooth(histcounts(tdata.PressDur(idxFPM),edges_HT,'Normalization','probability'),3)';
    t.HTdist_L = smooth(histcounts(tdata.PressDur(idxFPL),edges_HT,'Normalization','probability'),3)';

    outT = [outT;struct2table(t)];
end

end
function [TBT_fp_cor_noguessing, trial_idx, session_idx_, missing_session_idx] = plotRangeTrial(FP,TBT)
global lesionBoundary;
idx_fp = logical(sum(TBT.FP == FP,2));
TBT_fp_cor = TBT(idx_fp & TBT.Type==1,:);
TBT_fp_cor_noguessing = TBT_fp_cor(TBT_fp_cor.RT > 0.1,:);
lesionday_trial = find(TBT_fp_cor_noguessing.Date==lesionBoundary,1,'last');
prelesion_idx = (1:lesionday_trial)-(lesionday_trial+1);
postlesion_idx = (lesionday_trial+1:length(TBT_fp_cor_noguessing.RT))-lesionday_trial;
trial_idx = [prelesion_idx,postlesion_idx]';

[whichsession,idx_session_boundary,~] = unique(TBT_fp_cor_noguessing.Session,'stable');
session_idx = trial_idx(idx_session_boundary(2:end),:)-0.5;
session_idx(session_idx==0.5) = 0;
session_idx_ = setdiff(session_idx,0,'stable');

session_interval = diff(whichsession);
missing_session = session_interval(find(session_interval~=1))-1;
missing_session_idx = repelem(session_idx(find(session_interval~=1)),missing_session);
missing_session_idx_ = [];
for i = missing_session'
    missing_session_idx_ = [missing_session_idx_;(1:i)'];
end
missing_session_idx = missing_session_idx+missing_session_idx_*0.5;
if whichsession(1) ~= 1
    missing_session_idx_add = repelem(trial_idx(1),whichsession(1)-1)'-(1:(whichsession(1)-1))'*0.5;
    missing_session_idx = [missing_session_idx_add;missing_session_idx];
end
if whichsession(end) ~= TBT.Session(end)
    missing_session_idx_add = repelem(trial_idx(end),TBT.Session(end)-whichsession(end))'+(1:(TBT.Session(end)-whichsession(end)))'*0.5;
    missing_session_idx = [missing_session_idx;missing_session_idx_add];
end
end