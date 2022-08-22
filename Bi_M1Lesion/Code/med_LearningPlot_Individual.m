function med_LearningPlot_Individual(bAll,btAll,trialNearBound_LP)
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
% data session-by-session
SBS.Name = {};
SBS.Date = {};
SBS.Task = {};
SBS.nTrial = [];
SBS.rTrial = [];
SBS.Cor = [];
SBS.Pre = [];
SBS.Late = [];
SBS.RT = [];
SBS.RT_var = [];

% data trial-by-trial struct
TBTs.Name = {};
TBTs.Date = {};
TBTs.Task = {};
TBTs.iTrial = [];
TBTs.FP = [];
TBTs.Outcome = [];
TBTs.RT = [];
TBTs.PressDur = [];
t2c = [];
SBS.Dateincludeyear = [];

% central tendency & dispersion parameter
tend_disp_approachs = {'mean-std','mean-bootci','mean-sem','quartile',...
    'geomean-geomad','geomean-bootci'}; 
RT_stat = tend_disp_approachs{4};

for i=1:length(bAll)
    b = bAll(i);
    SBS.Name   = [SBS.Name;      b.Metadata.SubjectName];
    SBS.Date   = [SBS.Date;      b.Metadata.Date(end-3:end)];
    SBS.Dateincludeyear   = [SBS.Dateincludeyear;      b.Metadata.Date];
    idx_taskname  = strfind(b.Metadata.ProtocolName,'_');
    taskname = replace(b.Metadata.ProtocolName(idx_taskname(end)+1:end),{'Bpod','Three'},{'','3'});
    SBS.Task   = [SBS.Task;      taskname];
    numTrial   = length(b.Correct) + length(b.Premature) + length(b.Late);
    SBS.nTrial = [SBS.nTrial;    numTrial];
    SBS.rTrial = [SBS.rTrial;    numTrial./length(b.PressTime)];
    SBS.Cor    = [SBS.Cor;       length(b.Correct)./numTrial];
    SBS.Pre    = [SBS.Pre;       length(b.Premature)./numTrial];
    SBS.Late   = [SBS.Late;      length(b.Late)./numTrial];
    rt_summ  = compute_stat_summary(b.ReactionTime(b.ReactionTime>0),RT_stat);
    SBS.RT     = [SBS.RT;        rt_summ(1)];
    SBS.RT_var = [SBS.RT_var;    (rt_summ(3)-rt_summ(2))./2];
    
    bt = btAll{i}(btAll{i}.Type~=0,:);
    date_str = num2str(bt.Date);
    TBTs.Name     = [TBTs.Name;     bt.Subject];
    TBTs.Date     = [TBTs.Date;     string(date_str(:,end-3:end))];
    TBTs.Task     = [TBTs.Task;     bt.Task];
    TBTs.iTrial   = [TBTs.iTrial;   bt.iTrial];
    TBTs.FP       = [TBTs.FP;       bt.FP];
    TBTs.Outcome  = [TBTs.Outcome;  bt.Type];
    TBTs.RT       = [TBTs.RT;       bt.RT];
    TBTs.PressDur = [TBTs.PressDur; bt.PressDur];
    
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
end

% group index for sessions before or after lesion
SBS.Dateincludeyear=cellstr(SBS.Dateincludeyear);
t=datetime(SBS.Dateincludeyear,'InputFormat','yyyyMMdd');
maxtimeduration=max(caldays(caldiff(t)));
lesionday=find(caldays(caldiff(t))==maxtimeduration,1,'last');
groupindex=[zeros(lesionday,1);ones(length(SBS.Date)-lesionday,1)];

% trial-by-trial table
varNames = {'Name','Date','Task','iTrial','FP',...
    'Outcome','RT','PressDur'};
TBT = table(TBTs.Name,TBTs.Date,TBTs.Task,TBTs.iTrial,TBTs.FP,...
    TBTs.Outcome,TBTs.RT,TBTs.PressDur,'VariableNames',varNames);

ind_Cor = find(TBT.Outcome==1);

isExist3FPs = find(unique(TBT.Task) == "ThreeFPsMixedBpod");
if isExist3FPs
    idx_3fp = TBT.FP == 0.5 | TBT.FP == 1.0 | TBT.FP == 1.5;
    TBT_3fp_cor = TBT(idx_3fp & TBT.Outcome==1,:);
    TBT_3fp = TBT(idx_3fp,:);
    % some rat may even fail to warm up → bug in tirals for each FP
    TBT_3fp_fake.Name = repelem(TBT_3fp.Name(1),length(SBS.Date)*3,1);
    TBT_3fp_fake.Date = repelem(string(SBS.Date),3);
    TBT_3fp_fake.Task = repelem(TBT_3fp.Task(1),length(SBS.Date)*3,1);
    TBT_3fp_fake.iTrial = NaN(length(SBS.Date)*3,1);
    TBT_3fp_fake.FP = repmat([0.5;1;1.5],length(SBS.Date),1);
    TBT_3fp_fake.Outcome = NaN(length(SBS.Date)*3,1);
    TBT_3fp_fake.RT = NaN(length(SBS.Date)*3,1);
    TBT_3fp_fake.PressDur = NaN(length(SBS.Date)*3,1);
    TBT_3fp_fake = table(TBT_3fp_fake.Name,TBT_3fp_fake.Date,TBT_3fp_fake.Task,TBT_3fp_fake.iTrial,TBT_3fp_fake.FP,...
    TBT_3fp_fake.Outcome,TBT_3fp_fake.RT,TBT_3fp_fake.PressDur,'VariableNames',varNames);
    TBT_3fp_and_fake = [TBT_3fp;TBT_3fp_fake];

    %TBT_fp05_plot = plotRangeTrial(0.5,TBT,SBS,lesionday,trialNearBound_LP);
    %TBT_fp10_plot = plotRangeTrial(1.0,TBT,SBS,lesionday,trialNearBound_LP);
    %TBT_fp15_plot = plotRangeTrial(1.5,TBT,SBS,lesionday,trialNearBound_LP);
end
%% Plot
cDarkGray = [0.2,0.2,0.2];
cGreen = [0.4660 0.6740 0.1880];
cRed = [0.6350 0.0780 0.1840];
cYellow = [0.9290 0.6940 0.1250];
cBlue = [0,0.6902,0.9412];
cOrange = [0.929,0.49,0.192];
cGray = [0.4 0.4 0.4];
colorlist = [cDarkGray;cGray;cGreen;cRed;cYellow;cBlue];
color_FP = [cGray;mean([cOrange;cGray]);cOrange];

% Trial Number
g(1,1) = gramm('X', SBS.Date, 'Y',SBS.nTrial);
g(1,1).facet_grid([], SBS.Task, 'scale', 'free_x','space','free_x');
g(1,1).geom_point(); g(1,1).set_point_options('base_size',7);
g(1,1).geom_line();  g(1,1).set_line_options('base_size',2);

g(1,1).axe_property('ylim', [0 400], 'xticklabels', {}, 'XGrid', 'on', 'YGrid', 'on');
g(1,1).geom_vline('xintercept',lesionday+0.5,'style','k-.');
g(1,1).set_names('x',{}, 'y', 'Trials','column', '');
g(1,1).set_color_options('map',colorlist(1,:),'n_color',1,'n_lightness',1);
g(1,1).set_order_options('column',0,'X',0);

% Trial Press Ratio
g(2,1) = gramm('X', SBS.Date, 'Y',SBS.rTrial);
g(2,1).facet_grid([], cellstr(SBS.Task), 'scale', 'free_x','space','free_x', 'column_labels', false);
g(2,1).geom_point(); g(2,1).set_point_options('base_size',7);
g(2,1).geom_line();  g(2,1).set_line_options('base_size',2);
g(2,1).geom_vline('xintercept',lesionday+0.5,'style','k-.');
g(2,1).axe_property('ylim', [0 1], 'xticklabels', {}, 'XGrid', 'on', 'YGrid', 'on');
g(2,1).set_names('x',{}, 'y', 'Trial-press Ratio');
g(2,1).set_color_options('map',colorlist(2,:),'n_color',1,'n_lightness',1);
g(2,1).set_order_options('column',0,'X',0);

% Performance
g(3,1) = gramm('X',repmat(SBS.Date,3,1), 'Y',[SBS.Cor;SBS.Pre;SBS.Late],'color',repelem({'Correct';'Premature';'Late'},length(SBS.Date)));
g(3,1).facet_grid([], repmat(SBS.Task,3,1), 'scale', 'free_x','space','free_x');
g(3,1).geom_point(); g(3,1).set_point_options('base_size',7);
g(3,1).geom_line();  g(3,1).set_line_options('base_size',2);
g(3,1).geom_hline('yintercept',0.7,'style','k--');
g(3,1).geom_vline('xintercept',lesionday+0.5,'style','k-.');
g(3,1).axe_property('ylim', [0 1], 'XGrid', 'on', 'YGrid', 'on','XTickLabelRotation',90);
g(3,1).set_names('x','Date', 'y', 'Performance','column', '','color','');
g(3,1).set_color_options('map',colorlist([3,4,5],:),'n_color',3,'n_lightness',1);
g(3,1).set_order_options('column',0,'X',0,'color',{'Correct','Premature','Late'});
g(3,1).set_layout_options('legend_position',[0.05 0.22 0.08 0.11]);

%t2c
%g(2,2) = gramm('X',SBS.Date, 'Y',t2c);
%g(2,2).facet_grid([], cellstr(SBS.Task), 'scale', 'free_x','space','free_x', 'column_labels', false);
%g(2,2).geom_point(); g(2,2).set_point_options('base_size',7);
%g(2,2).geom_line();  g(2,2).set_line_options('base_size',2);
%g(2,2).axe_property('ylim', [0 300], 'XGrid', 'on', 'YGrid', 'on');
%g(2,2).set_names('x','Date', 'y', 'Trial to Criteria');
%g(2,2).set_color_options('map',cOrange,'n_color',1,'n_lightness',1);
%g(2,2).set_order_options('column',0,'X',0);

% Premature
%g(4,1) = gramm('X', SBS.Date, 'Y',SBS.Pre);
%g(4,1).facet_grid([], cellstr(SBS.Task), 'scale', 'free_x','space','free_x', 'column_labels', false);
%g(4,1).geom_point(); g(4,1).set_point_options('base_size',7);
%g(4,1).geom_line();  g(4,1).set_line_options('base_size',2);
%g(4,1).axe_property('ylim', [0 0.8], 'xticklabels', {}, 'XGrid', 'on', 'YGrid', 'on');
%g(4,1).set_names('x',{}, 'y', 'Premature');
%g(4,1).set_color_options('map',colorlist(4,:),'n_color',1,'n_lightness',1);
%g(4,1).set_order_options('column',0);

% Late
%g(5,1) = gramm('X', SBS.Date, 'Y',SBS.Late);
%g(5,1).facet_grid([], cellstr(SBS.Task), 'scale', 'free_x','space','free_x', 'column_labels', false);
%g(5,1).geom_point(); g(5,1).set_point_options('base_size',7);
%g(5,1).geom_line();  g(5,1).set_line_options('base_size',2);
%g(5,1).axe_property('ylim', [0 0.3], 'xticklabels', {}, 'XGrid', 'on', 'YGrid', 'on');
%g(5,1).set_names('x',{}, 'y', 'Late');
%g(5,1).set_color_options('map',colorlist(5,:),'n_color',1,'n_lightness',1);
%g(5,1).set_order_options('column',0);

% RT
g(1,2) = gramm('X',cellstr(TBT.Date(ind_Cor)),'Y',TBT.RT(ind_Cor));
g(1,2).facet_grid([], cellstr(TBT.Task(ind_Cor)), 'scale', 'free_x','space','free_x', 'column_labels', false);
g(1,2).stat_summary('type', @(x)compute_stat_summary(x,RT_stat),'geom',{'black_errorbar','point','line'});
g(1,2).geom_vline('xintercept',lesionday+0.5,'style','k-.');
g(1,2).axe_property('ylim', [0 0.7], 'XGrid', 'on', 'YGrid', 'on');
g(1,2).set_names('x','Date', 'y', 'RT(s)');
g(1,2).set_color_options('map',colorlist(6,:),'n_color',1,'n_lightness',1);
g(1,2).set_order_options('column',0,'X',0);

if isExist3FPs
    g(1,2).axe_property('xticklabels',{});
    g(1,2).set_names('x',{}, 'y', 'RT(s)');    
    
    g(2,2) = gramm('X',cellstr(TBT_3fp_cor.Date),'Y',TBT_3fp_cor.RT,'color',TBT_3fp_cor.FP);
    g(2,2).facet_grid([], cellstr(TBT_3fp_cor.Task), 'scale', 'free_x','space','free_x', 'column_labels', false);
    g(2,2).stat_summary('type', @(x)compute_stat_summary(x,RT_stat),'geom',{'errorbar','point'},'dodge',0.5);
    g(2,2).geom_vline('xintercept',lesionday+0.5,'style','k-.');
    g(2,2).axe_property('ylim', [0 0.7], 'XGrid', 'on', 'YGrid', 'on','XTickLabelRotation',90);
    g(2,2).set_names('x','Date', 'y', 'RT_3FPs(s)','color','FP(s)');
    g(2,2).set_color_options('map',color_FP,'n_color',3,'n_lightness',1);
    g(2,2).set_order_options('column',0,'X',0);
    g(2,2).no_legend();

    SBS_3fp_statarray = grpstats(removevars(TBT_3fp_and_fake,{'Name','Task'}),{'Date','FP'});
    SBS_3fp_nTrial = reshape(SBS_3fp_statarray.GroupCount-1,3,[])';
    g(1,1) = gramm('X',repmat(SBS.Date,4,1), 'Y',[SBS_3fp_nTrial(:,1);SBS_3fp_nTrial(:,2);SBS_3fp_nTrial(:,3);SBS.nTrial],'color',repelem({['FP=',mat2str(SBS_3fp_statarray.FP(1)),'s'];['FP=',mat2str(SBS_3fp_statarray.FP(2)),'s'];['FP=',mat2str(SBS_3fp_statarray.FP(3)),'s'];'Total'},length(SBS.Date)));
    g(1,1).facet_grid([], repmat(SBS.Task,4,1), 'scale', 'free_x','space','free_x');
    g(1,1).geom_point(); g(1,1).set_point_options('base_size',7);
    g(1,1).geom_line();  g(1,1).set_line_options('base_size',2);
    g(1,1).axe_property('ylim', [0 400], 'xticklabels', {}, 'XGrid', 'on', 'YGrid', 'on');
    g(1,1).geom_vline('xintercept',lesionday+0.5,'style','k-.');
    g(1,1).set_names('x',{}, 'y', 'Trials','column', '','color','');
    g(1,1).set_color_options('map',[color_FP;cDarkGray],'n_color',4,'n_lightness',1);
    g(1,1).set_order_options('column',0,'X',0,'color',{['FP=',mat2str(SBS_3fp_statarray.FP(1)),'s'];['FP=',mat2str(SBS_3fp_statarray.FP(2)),'s'];['FP=',mat2str(SBS_3fp_statarray.FP(3)),'s'];'Total'});
    g(1,1).set_layout_options('legend_position',[0.05 0.84 0.08 0.11]);

   %g(3,2) = gramm('X',repmat([-50:-1,1:50]',3,1), 'Y',[TBT_fp05_plot.RT;TBT_fp10_plot.RT;TBT_fp15_plot.RT],'color',repelem({'FP=0.5s';'FP=1s';'FP=1.5s'},length(TBT_fp05_plot.RT)));
   %g(3,2).facet_grid([], repmat(SBS.Task,3,1), 'scale', 'free_x','space','free_x');
   %g(3,2).geom_point(); g(3,2).set_point_options('base_size',3);
   %g(3,2).geom_vline('xintercept',0,'style','k-.');
   %g(3,2).axe_property('ylim', [0 0.7], 'XGrid', 'on', 'YGrid', 'on');
   %g(3,2).set_names('x','# of trial balanced', 'y', 'RT(s)','column', '','color','');
   %g(3,2).set_color_options('map',color_FP,'n_color',3,'n_lightness',1);
   %g(3,2).set_order_options('column',0,'X',0,'color',{'FP=0.5s','FP=1s','FP=1.5s'});
   %g(3,2).set_layout_options('legend_position',[0.55 0.22 0.08 0.11]);
end

% % RT (alternatives)
% g(5,1) = gramm('X', SBS.Date, 'Y',SBS.RT);
% g(5,1).facet_grid([], cellstr(SBS.Task), 'scale', 'free_x','space','free_x', 'column_labels', false);
% g(5,1).geom_point(); g(5,1).set_point_options('base_size',7);
% g(5,1).geom_line();  g(5,1).set_line_options('base_size',2);
% g(5,1).axe_property('ylim', [100 500], 'xticklabels', {}, 'XGrid', 'on', 'YGrid', 'on');
% g(5,1).set_names('x',{}, 'y', 'RT(ms)');
% % g(5,1).set_title('RT');
% g(5,1).set_color_options('map',colorlist(5,:),'n_color',1,'n_lightness',1);
% 
% % RT variance
% g(6,1) = gramm('X', SBS.Date, 'Y',SBS.RT_var);
% g(6,1).facet_grid([], cellstr(SBS.Task), 'scale', 'free_x','space','free_x', 'column_labels', false);
% g(6,1).geom_point(); g(6,1).set_point_options('base_size',7);
% g(6,1).geom_line();  g(6,1).set_line_options('base_size',2);
% g(6,1).axe_property('ylim', [0 400], 'XGrid', 'on', 'YGrid', 'on');
% g(6,1).set_names('x','Date', 'y', 'RT Variance(ms)');
% % g(6,1).set_title('RT variance');
% g(6,1).set_color_options('map',colorlist(6,:),'n_color',1,'n_lightness',1);

figure(3);clf(3);
set(gcf, 'Name', 'Learning curve', 'unit', 'normalized', 'position',[0 0 1 1], 'paperpositionmode', 'auto' )
g.set_title(b.Metadata.SubjectName + ": Learning Curve"+"  ("+RT_stat+")");
g.draw();
figName = ['LearningProgress_', b.Metadata.SubjectName];
%% Save
figPath = fullfile(pwd,'AnalysisFig');
if ~exist(figPath,'dir')
    mkdir(figPath);
end
figFile = fullfile(figPath,figName);
saveas(gcf, figFile, 'png');
saveas(gcf, figFile, 'fig');

end
%% Function
function TBT_fp_plot = plotRangeTrial(FP,TBT,SBS,lesionday,trialNearBound)
idx_fp = TBT.FP == FP;
TBT_fp_cor = TBT(idx_fp & TBT.Outcome==1,:);
lessionday_trial = find(TBT_fp_cor.Date==string(SBS.Date(lesionday)),1,'last');
plotRange_trial = lessionday_trial-trialNearBound+1:lessionday_trial+trialNearBound;
TBT_fp_plot = TBT_fp_cor(plotRange_trial,:);
end