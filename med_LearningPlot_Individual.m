function med_LearningPlot_Individual(btAll2d,plotRange)
% _________________________________________________________________________
% File:               med_LearningPlot_Individual.m
% Created on:         Sept 27, 2021
% Created by:         Yu Chen
% Last revised on:    Oct 18, 2021
% Last revised by:    Yu Chen
% _________________________________________________________________________
% Required Packages:
% 'gramm' by Pierre Morel
% _________________________________________________________________________
%% Data packaging
% central tendency & dispersion parameter
tend_disp_approachs = {'mean-std','mean-bootci','mean-sem','quartile',...
    'geomean-geomad','geomean-bootci'}; 
RT_stat = tend_disp_approachs{4};
altMethod = {'mean','median','geomean'};
cenMethod = altMethod{2};
grandCen = altMethod{1};
grandErr = 'sem';

global lesionBoundary;
btAll2d_use = {};
missSess = 0;
for i=1:size(plotRange,1)
    if find(cellfun(@(x)~isempty(x),btAll2d(i,:)),1,'last') < plotRange{i,1}(end)
        temp_missSess = plotRange{i,1}(end)-find(cellfun(@(x)~isempty(x),btAll2d(i,:)),1,'last');
        if temp_missSess > missSess
            missSess = temp_missSess;
        end
    end
end
for i=1:size(plotRange,1)
    plotRange{i,1} = plotRange{i,1}(1):plotRange{i,1}(end)-missSess;
    btAll2d_use(i,1:length(btAll2d(i,plotRange{i,1}))) = btAll2d(i,plotRange{i,1});
end

[SBS,temp_TBT] = packData(btAll2d_use,cenMethod);
TBT = temp_TBT(temp_TBT.Type~=0,:);
TBT = removevars(TBT,{'Group','PressTime','ReleaseTime','ToneTime','RW'});
TBT = movevars(TBT,{'FP','Type','RT'},'Before','PressDur');
save SBS SBS;
save TBT TBT;

SBS_sess = grpstats(removevars(SBS,{'Subject','Period','Date','Task'}),'Session',{grandCen,grandErr});


%for i=1:length(bAll)
    %b = bAll(i);
    %SBS.Name   = [SBS.Name;      b.Metadata.SubjectName];
    %SBS.Date   = [SBS.Date;      b.Metadata.Date(end-3:end)];
    %SBS.Dateincludeyear   = [SBS.Dateincludeyear;      b.Metadata.Date];
    %idx_taskname  = strfind(b.Metadata.ProtocolName,'_');
    %taskname = replace(b.Metadata.ProtocolName(idx_taskname(end)+1:end),{'Bpod','Three'},{'','3'});
    %SBS.Task   = [SBS.Task;      taskname];
    %numTrial   = length(b.Correct) + length(b.Premature) + length(b.Late);
    %SBS.nTrial = [SBS.nTrial;    numTrial];
    %SBS.rTrial = [SBS.rTrial;    numTrial./length(b.PressTime)];
    %SBS.Cor    = [SBS.Cor;       length(b.Correct)./numTrial];
    %SBS.Pre    = [SBS.Pre;       length(b.Premature)./numTrial];
    %SBS.Late   = [SBS.Late;      length(b.Late)./numTrial];
    %rt_summ  = compute_stat_summary(b.ReactionTime(b.ReactionTime>0),RT_stat);
    %SBS.RT     = [SBS.RT;        rt_summ(1)];
    %SBS.RT_var = [SBS.RT_var;    (rt_summ(3)-rt_summ(2))./2];
    
    %bt = btAll{i}(btAll{i}.Type~=0,:);
    %date_str = num2str(bt.Date);
    %TBTs.Name     = [TBTs.Name;     bt.Subject];
    %TBTs.Date     = [TBTs.Date;     string(date_str(:,end-3:end))];
    %TBTs.Task     = [TBTs.Task;     bt.Task];
    %TBTs.iTrial   = [TBTs.iTrial;   bt.iTrial];
    %TBTs.FP       = [TBTs.FP;       bt.FP];
    %TBTs.Outcome  = [TBTs.Outcome;  bt.Type];
    %TBTs.RT       = [TBTs.RT;       bt.RT];
    %TBTs.PressDur = [TBTs.PressDur; bt.PressDur];
    
    %if taskname=='Wait1'
        %t2c_thissession=find(abs(bt.FP-1.5)<1E-14,1);
        %if isempty(t2c_thissession)
            %t2c=[t2c,NaN];
        %else
            %t2c=[t2c,t2c_thissession];
        %end
    %else
        %t2c_thissession=find(abs(bt.FP-1.5)<1E-14 & abs(bt.RW-0.6)<1E-14,1);
        %if isempty(t2c_thissession)
            %t2c=[t2c,NaN];
        %else
            %t2c=[t2c,t2c_thissession];
        %end
    %end
%end

% group index for sessions before or after lesion
%SBS.Dateincludeyear=cellstr(SBS.Dateincludeyear);
%t=datetime(SBS.Dateincludeyear,'InputFormat','yyyyMMdd');
%maxtimeduration=max(caldays(caldiff(t)));
%lesionday=find(caldays(caldiff(t))==maxtimeduration,1,'last');
%groupindex=[zeros(lesionday,1);ones(length(SBS.Date)-lesionday,1)];

% trial-by-trial table
varNames = {'Subject','Date','Session','Task','iTrial','FP',...
    'Type','RT','PressDur'};
%TBT = table(TBTs.Name,TBTs.Date,TBTs.Task,TBTs.iTrial,TBTs.FP,...
%    TBTs.Outcome,TBTs.RT,TBTs.PressDur,'VariableNames',varNames);

ind_Cor = find(TBT.Type==1);

isExist3FPs = find(unique(TBT.Task) == "ThreeFPsMixedBpod");
%if isExist3FPs
%    idx_3fp = TBT.FP == 0.5 | TBT.FP == 1.0 | TBT.FP == 1.5;
%    TBT_3fp_cor = TBT(idx_3fp & TBT.Type==1,:);
%    TBT_3fp = TBT(idx_3fp,:);
%    % some rat may even fail to warm up → bug in tirals for each FP
%    subject_use = unique(SBS.Subject,'stable');
%    TBT_3fp_fake.Subject = repelem(subject_use,cellfun(@(x)length(x),plotRange)*3);
%    TBT_3fp_fake.Date = repelem(SBS.Date,3);
%    TBT_3fp_fake.Session = repelem(SBS.Session,3);
%    TBT_3fp_fake.Task = repelem(TBT_3fp.Task(1),length(SBS.Date)*3,1);
%    TBT_3fp_fake.iTrial = NaN(length(SBS.Date)*3,1);
%    TBT_3fp_fake.FP = repmat([0.5;1;1.5],length(SBS.Date),1);
%    TBT_3fp_fake.Type = NaN(length(SBS.Date)*3,1);
%    TBT_3fp_fake.RT = NaN(length(SBS.Date)*3,1);
%    TBT_3fp_fake.PressDur = NaN(length(SBS.Date)*3,1);
%    TBT_3fp_fake = table(TBT_3fp_fake.Subject,TBT_3fp_fake.Date,TBT_3fp_fake.Session,TBT_3fp_fake.Task,TBT_3fp_fake.iTrial,TBT_3fp_fake.FP,...
%    TBT_3fp_fake.Type,TBT_3fp_fake.RT,TBT_3fp_fake.PressDur,'VariableNames',varNames);
%    TBT_3fp_and_fake = [TBT_3fp;TBT_3fp_fake];

    %TBT_fp05_plot = plotRangeTrial(0.5,TBT,SBS,lesionday,trialNearBound_LP);
    %TBT_fp10_plot = plotRangeTrial(1.0,TBT,SBS,lesionday,trialNearBound_LP);
    %TBT_fp15_plot = plotRangeTrial(1.5,TBT,SBS,lesionday,trialNearBound_LP);
%end
%% Plot
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

lesionday = SBS.Session(find(SBS.Period=="Pre",1,'last'));

hf = figure(43); clf(hf,'reset');
set(hf, 'name', 'Lesion effect', 'units', 'centimeters', 'position', [1 1 13 10],...
    'PaperPositionMode', 'auto','renderer','painter'); % 生科论文要求版面裁掉边距还剩，宽度15.8cm,高度24.2cm

size1 = [4 2.5];
size2 = [4 2];
ys = [6.5 3.75 1]; % yStart
xs1 = [1.5 7]; % xStart

session_use = unique(SBS.Session);

% Trial Number
ha11 = axes;
set(ha11, 'units', 'centimeters', 'position', [xs1(1) ys(1) size1], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0.5,SBS.Session(end)+0.5],'xtick',1:SBS.Session(end),'xticklabel',{},'fontsize',7, ...
    'ylim',[0 400],'ytick',0:100:400,'yticklabels',{'0','100','200','300','400'},'ticklength', [0.02 0.025]);
nTrial = plot(SBS_sess.Session, SBS_sess.mean_nTrial, 'linewidth', 2, 'color', 'k', 'markerfacecolor', 'k', 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_nTrial'-SBS_sess.sem_nTrial'; SBS_sess.mean_nTrial'+SBS_sess.sem_nTrial'], ...
    'color','k', 'linewidth', 0.6)
nTrial_S = plot(SBS_sess.Session, SBS_sess.mean_nTrial_S,  'linewidth', 0.6, 'color', 'k', 'markerfacecolor', 'k', 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_nTrial_S'-SBS_sess.sem_nTrial_S'; SBS_sess.mean_nTrial_S'+SBS_sess.sem_nTrial_S'], ...
    'color','k', 'linewidth', 0.6)
nTrial_M = plot(SBS_sess.Session, SBS_sess.mean_nTrial_M, 'linewidth', 1, 'color', 'k', 'markerfacecolor', 'k', 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_nTrial_M'-SBS_sess.sem_nTrial_M'; SBS_sess.mean_nTrial_M'+SBS_sess.sem_nTrial_M'], ...
    'color','k', 'linewidth', 0.6)
nTrial_L = plot(SBS_sess.Session, SBS_sess.mean_nTrial_L,  'linewidth', 1.5, 'color', 'k', 'markerfacecolor', 'k', 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_nTrial_L'-SBS_sess.sem_nTrial_L'; SBS_sess.mean_nTrial_L'+SBS_sess.sem_nTrial_L'], ...
    'color','k', 'linewidth', 0.6)
plot([lesionday+0.5,lesionday+0.5],[0,400],':k','linewidth',0.6);
ylabel('Trials','Fontsize',7,'FontName','Arial')

le1 = legend([nTrial,nTrial_S,nTrial_M,nTrial_L],{'Total','FP 0.5 s','FP 1.0 s','FP 1.5 s'},...
    'Fontsize',7,'FontName','Arial','units','centimeters',...
    'Position',[11.5,7.5,1,1]); % [10.7,8.9,1,1]
legend('boxoff');
le1.ItemTokenSize = [12,22];
le1.Position = le1.Position + [0.025 0.045 0 0];

ha21 = axes;
set(ha21, 'units', 'centimeters', 'position', [xs1(1) ys(2) size1], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0.5,SBS.Session(end)+0.5],'xtick',1:SBS.Session(end),'xticklabel',{},'fontsize',7, ...
    'ylim',[0.6 1],'ytick',0.6:0.1:1,'yticklabels',{'0.6','0.7','0.8','0.9','1.0'},'ticklength', [0.02 0.025]);
rTrial = plot(SBS_sess.Session, SBS_sess.mean_rTrial, 'linewidth', 2, 'color', cGray2, 'markerfacecolor', cGray2, 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_rTrial'-SBS_sess.sem_rTrial'; SBS_sess.mean_rTrial'+SBS_sess.sem_rTrial'], ...
    'color',cGray2, 'linewidth', 0.6)
plot([lesionday+0.5,lesionday+0.5],[0.6,1],':k','linewidth',0.6);
ylabel('Trial-press Ratio','Fontsize',7,'FontName','Arial')

ha31 = axes;
set(ha31, 'units', 'centimeters', 'position', [xs1(1) ys(3) size1], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0.5,SBS.Session(end)+0.5],'xtick',1:SBS.Session(end),'xticklabel',cellstr(string([(-lesionday:-1),(1:session_use(end)-lesionday)])),'fontsize',7,'xticklabelRotation',0, ...
    'ylim',[0 1],'ytick',0:0.2:1,'yticklabels',{'0','20','40','60','80','100'},'ticklength', [0.02 0.025]);
Cor = plot(SBS_sess.Session, SBS_sess.mean_Cor, 'linewidth', 2, 'color',cGreen, 'markerfacecolor', cGreen, 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_Cor'-SBS_sess.sem_Cor'; SBS_sess.mean_Cor'+SBS_sess.sem_Cor'], ...
    'color',cGreen, 'linewidth', 0.6)
Cor_S = plot(SBS_sess.Session, SBS_sess.mean_Cor_S,  'linewidth', 0.6, 'color', cGreen, 'markerfacecolor', cGreen, 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_Cor_S'-SBS_sess.sem_Cor_S'; SBS_sess.mean_Cor_S'+SBS_sess.sem_Cor_S'], ...
    'color',cGreen, 'linewidth', 0.6)
Cor_M = plot(SBS_sess.Session, SBS_sess.mean_Cor_M, 'linewidth', 1, 'color', cGreen, 'markerfacecolor', cGreen, 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_Cor_M'-SBS_sess.sem_Cor_M'; SBS_sess.mean_Cor_M'+SBS_sess.sem_Cor_M'], ...
    'color',cGreen, 'linewidth', 0.6)
Cor_L = plot(SBS_sess.Session, SBS_sess.mean_Cor_L,  'linewidth', 1.5, 'color', cGreen, 'markerfacecolor', cGreen, 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_Cor_L'-SBS_sess.sem_Cor_L'; SBS_sess.mean_Cor_L'+SBS_sess.sem_Cor_L'], ...
    'color',cGreen, 'linewidth', 0.6)
plot([lesionday+0.5,lesionday+0.5],[0,1],':k','linewidth',0.6);
ylabel('Correct(%)','Fontsize',7,'FontName','Arial')
xlabel('# of session','Fontsize',7,'FontName','Arial')

ha22 = axes;
set(ha22, 'units', 'centimeters', 'position', [xs1(2) ys(2) size1], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0.5,SBS.Session(end)+0.5],'xtick',1:SBS.Session(end),'xticklabel',{},'fontsize',7, ...
    'ylim',[0 1],'ytick',0:0.2:1,'yticklabels',{'0','20','40','60','80','100'},'ticklength', [0.02 0.025]);
Pre = plot(SBS_sess.Session, SBS_sess.mean_Pre, 'linewidth', 2, 'color',cRed, 'markerfacecolor', cRed, 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_Pre'-SBS_sess.sem_Pre'; SBS_sess.mean_Pre'+SBS_sess.sem_Pre'], ...
    'color',cRed, 'linewidth', 0.6)
Pre_S = plot(SBS_sess.Session, SBS_sess.mean_Pre_S,  'linewidth', 0.6, 'color', cRed, 'markerfacecolor', cRed, 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_Pre_S'-SBS_sess.sem_Pre_S'; SBS_sess.mean_Pre_S'+SBS_sess.sem_Pre_S'], ...
    'color',cRed, 'linewidth', 0.6)
Pre_M = plot(SBS_sess.Session, SBS_sess.mean_Pre_M, 'linewidth', 1, 'color', cRed, 'markerfacecolor', cRed, 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_Pre_M'-SBS_sess.sem_Pre_M'; SBS_sess.mean_Pre_M'+SBS_sess.sem_Pre_M'], ...
    'color',cRed, 'linewidth', 0.6)
Pre_L = plot(SBS_sess.Session, SBS_sess.mean_Pre_L,  'linewidth', 1.5, 'color', cRed, 'markerfacecolor', cRed, 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_Pre_L'-SBS_sess.sem_Pre_L'; SBS_sess.mean_Pre_L'+SBS_sess.sem_Pre_L'], ...
    'color',cRed, 'linewidth', 0.6)
plot([lesionday+0.5,lesionday+0.5],[0,1],':k','linewidth',0.6);
ylabel('Premature(%)','Fontsize',7,'FontName','Arial')

ha12 = axes;
set(ha12, 'units', 'centimeters', 'position', [xs1(2) ys(1) size1], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0.5,SBS.Session(end)+0.5],'xtick',1:SBS.Session(end),'xticklabel',{},'fontsize',7, ...
    'ylim',[0.1 0.5],'ytick',0.1:0.1:0.5,'yticklabels',{'0.1','0.2','0.3','0.4','0.5'},'ticklength', [0.02 0.025]);
RT = plot(SBS_sess.Session, SBS_sess.mean_RT, 'linewidth', 2, 'color', cBlue, 'markerfacecolor', cBlue, 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_RT'-SBS_sess.sem_RT'; SBS_sess.mean_RT'+SBS_sess.sem_RT'], ...
    'color',cBlue, 'linewidth', 0.6)
plot([lesionday+0.5,lesionday+0.5],[0.1,0.5],':k','linewidth',0.6);
ylabel('RT(s)','Fontsize',7,'FontName','Arial')

ha32 = axes;
set(ha32, 'units', 'centimeters', 'position', [xs1(2) ys(3) size1], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0.5,SBS.Session(end)+0.5],'xtick',1:SBS.Session(end),'xticklabel',cellstr(string([(-lesionday:-1),(1:session_use(end)-lesionday)])),'fontsize',7,'xticklabelRotation',0, ...
    'ylim',[0 1],'ytick',0:0.2:1,'yticklabels',{'0','20','40','60','80','100'},'ticklength', [0.02 0.025]);
Late = plot(SBS_sess.Session, SBS_sess.mean_Late, 'linewidth', 2, 'color',cGray, 'markerfacecolor', cGray, 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_Late'-SBS_sess.sem_Late'; SBS_sess.mean_Late'+SBS_sess.sem_Late'], ...
    'color',cGray, 'linewidth', 0.6)
Late_S = plot(SBS_sess.Session, SBS_sess.mean_Late_S,  'linewidth', 0.6, 'color', cGray, 'markerfacecolor', cGray, 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_Late_S'-SBS_sess.sem_Late_S'; SBS_sess.mean_Late_S'+SBS_sess.sem_Late_S'], ...
    'color',cGray, 'linewidth', 0.6)
Late_M = plot(SBS_sess.Session, SBS_sess.mean_Late_M, 'linewidth', 1, 'color', cGray, 'markerfacecolor', cGray, 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_Late_M'-SBS_sess.sem_Late_M'; SBS_sess.mean_Late_M'+SBS_sess.sem_Late_M'], ...
    'color',cGray, 'linewidth', 0.6)
Late_L = plot(SBS_sess.Session, SBS_sess.mean_Late_L,  'linewidth', 1.5, 'color', cGray, 'markerfacecolor', cGray, 'markeredgecolor','w', 'markersize', 5);
line([SBS_sess.Session'; SBS_sess.Session'], [SBS_sess.mean_Late_L'-SBS_sess.sem_Late_L'; SBS_sess.mean_Late_L'+SBS_sess.sem_Late_L'], ...
    'color',cGray, 'linewidth', 0.6)
plot([lesionday+0.5,lesionday+0.5],[0,1],':k','linewidth',0.6);
ylabel('Late(%)','Fontsize',7,'FontName','Arial')
xlabel('# of session','Fontsize',7,'FontName','Arial')


%figure(3);clf(3);
%set(gcf, 'Name', 'Learning curve', 'unit', 'normalized', 'position',[0 0 1 1], 'paperpositionmode', 'auto' )
%g.set_title(subject_use(1)+"2"+ subject_use(end)+ ": Learning Curve"+"  ("+cenMethod+" & "+grandErr+")");
%g.draw();
%figName = ['LearningProgress_', char(subject_use(1)),'2', char(subject_use(end))];
%% Save
subject_use = unique(SBS.Subject,'stable');
figName = ['MiceLesion_SBS_Curve_', char(subject_use(1)),'2', char(subject_use(end))];
savename = fullfile(pwd,figName);
saveas(hf,savename,'fig');
print(hf,'-dpng',savename);
print(hf,'-depsc2',savename);

end
%% Functions
function TBT_fp_plot = plotRangeTrial(FP,TBT,SBS,lesionday,trialNearBound)
idx_fp = TBT.FP == FP;
TBT_fp_cor = TBT(idx_fp & TBT.Outcome==1,:);
lessionday_trial = find(TBT_fp_cor.Date==string(SBS.Date(lesionday)),1,'last');
plotRange_trial = lessionday_trial-trialNearBound+1:lessionday_trial+trialNearBound;
TBT_fp_plot = TBT_fp_cor(plotRange_trial,:);
end

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
t.nTrial_S = sum(abs(tdata.FP-0.5)<1e-4);
t.nTrial_M = sum(abs(tdata.FP-1.0)<1e-4);
t.nTrial_L = sum(abs(tdata.FP-1.5)<1e-4);
t.Dark   = sum(data.Type==0)./size(data,1);
t.rTrial = 1-t.Dark;
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
        t.RT_S = mean(tdata(tdata.Type==1 & abs(tdata.FP-0.5)<1e-4,:).RT);
        t.RT_M = mean(tdata(tdata.Type==1 & abs(tdata.FP-1.0)<1e-4,:).RT);
        t.RT_L = mean(tdata(tdata.Type==1 & abs(tdata.FP-1.5)<1e-4,:).RT);
    case 'median'
        t.HT = median(tdata.PressDur);
        t.RT = median(tdata(tdata.Type==1,:).RT);
        t.RT_S = median(tdata(tdata.Type==1 & abs(tdata.FP-0.5)<1e-4,:).RT);
        t.RT_M = median(tdata(tdata.Type==1 & abs(tdata.FP-1.0)<1e-4,:).RT);
        t.RT_L = median(tdata(tdata.Type==1 & abs(tdata.FP-1.5)<1e-4,:).RT);
    case 'geomean'
        t.HT = geomean(tdata.PressDur);
        t.RT = geomean(tdata(tdata.Type==1 & tdata.RT>0,:).RT);
        t.RT_S = geomean(tdata(tdata.Type==1 & tdata.RT>0 & abs(tdata.FP-0.5)<1e-4,:).RT);
        t.RT_M = geomean(tdata(tdata.Type==1 & tdata.RT>0 & abs(tdata.FP-1.0)<1e-4,:).RT);
        t.RT_L = geomean(tdata(tdata.Type==1 & tdata.RT>0 & abs(tdata.FP-1.5)<1e-4,:).RT);
end
outT = [outT;struct2table(t)];

end