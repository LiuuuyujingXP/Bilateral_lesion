
%%Data
Animalfiles = arrayfun(@(x)x.name,dir('*.mat'),'UniformOutput',false);
Perf_saline = cell(length(Animalfiles),1);
Perf_DCZ = cell(length(Animalfiles),1);
PrDur_cued_saline = cell(length(Animalfiles),1);
PrDur_uncued_saline = cell(length(Animalfiles),1);
PrDur_cued_DCZ = cell(length(Animalfiles),1);
PrDur_uncued_DCZ = cell(length(Animalfiles),1);
PrDur_distr_cued_saline = cell(length(Animalfiles),1);
PrDur_distr_uncued_saline = cell(length(Animalfiles),1);
PrDur_distr_cued_DCZ = cell(length(Animalfiles),1);
PrDur_distr_uncued_DCZ = cell(length(Animalfiles),1);

AfterFP_cued_saline = cell(length(Animalfiles),1);
AfterFP_uncued_saline = cell(length(Animalfiles),1);
AfterFP_cued_DCZ = cell(length(Animalfiles),1);
AfterFP_uncued_DCZ = cell(length(Animalfiles),1);
AfterFP_distr_cued_saline = cell(length(Animalfiles),1);
AfterFP_distr_uncued_saline = cell(length(Animalfiles),1);
AfterFP_distr_cued_DCZ = cell(length(Animalfiles),1);
AfterFP_distr_uncued_DCZ = cell(length(Animalfiles),1);
edges = [0:50:5000];
edge_centers = mean(edges(1:end-1), edges(2:end));
FP = 1500;
Animalnames = cell(1,length(Animalfiles));
for n = 1:length(Animalfiles) 
    b_used = load(Animalfiles{n}); 
    Animalnames{n} =b_used.bmix(1).Metadata.SubjectName; 
    drugind=[];
    Perf_saline{n} = zeros(2,4);
    Perf_DCZ{n} = zeros(2,4);
    for i= 1:length(b_used.bmix)
        DCZfind=strfind(b_used.bmix(i).Experiment,'DCZ');
        drugind=[drugind ~isempty( DCZfind)];
    end
    DCZind=find(drugind==1);
    Salineind=find(drugind==0);
    
    for i = 1:length(Salineind)
        Perf_saline{n} = Perf_saline{n} + b_used.bmix(Salineind(i)).Performance(:,1:4);
        PrDur_cued_saline{n} = [PrDur_cued_saline{n},b_used.bmix(Salineind(i)).PrDurCued];
        PrDur_uncued_saline{n} = [PrDur_uncued_saline{n},b_used.bmix(Salineind(i)).PrDurUncued];
    end

    for i = 1:length(DCZind)
        Perf_DCZ{n} = Perf_DCZ{n} + b_used.bmix(DCZind(i)).Performance(:,1:4);
        PrDur_cued_DCZ{n} = [PrDur_cued_DCZ{n},b_used.bmix(DCZind(i)).PrDurCued];
        PrDur_uncued_DCZ{n} = [PrDur_uncued_DCZ{n},b_used.bmix(DCZind(i)).PrDurUncued];
    end
    
    Releasetime_CS = PrDur_cued_saline{n}-FP;
    AfterFP_cued_saline{n} = Releasetime_CS(Releasetime_CS>0);
    Releasetime_US = PrDur_uncued_saline{n}-FP;
    AfterFP_uncued_saline{n} = Releasetime_US(Releasetime_US>0);
    Releasetime_CD = PrDur_cued_DCZ{n}-FP;
    AfterFP_cued_DCZ{n} = Releasetime_CD(Releasetime_CD>0);
    Releasetime_UD = PrDur_uncued_DCZ{n}-FP;
    AfterFP_uncued_DCZ{n} = Releasetime_UD(Releasetime_UD>0);
    
    AfterFP_distr_cued_saline{n} = histcounts(AfterFP_cued_saline{n}, edges, 'Normalization','probability');
    AfterFP_distr_cued_saline{n} = smoothdata(AfterFP_distr_cued_saline{n}, 'gaussian' , 5);
    AfterFP_distr_uncued_saline{n} = histcounts(AfterFP_uncued_saline{n}, edges, 'Normalization','probability');
    AfterFP_distr_uncued_saline{n} = smoothdata(AfterFP_distr_uncued_saline{n}, 'gaussian' , 5);
    AfterFP_distr_cued_DCZ{n} = histcounts(AfterFP_cued_DCZ{n}, edges, 'Normalization','probability');
    AfterFP_distr_cued_DCZ{n} = smoothdata(AfterFP_distr_cued_DCZ{n}, 'gaussian' , 5);
    AfterFP_distr_uncued_DCZ{n} = histcounts(AfterFP_uncued_DCZ{n}, edges, 'Normalization','probability');
    AfterFP_distr_uncued_DCZ{n} = smoothdata(AfterFP_distr_uncued_DCZ{n}, 'gaussian' , 5);

    PrDur_distr_cued_saline{n} = histcounts(PrDur_cued_saline{n}, edges, 'Normalization','probability');
    PrDur_distr_cued_saline{n} = smoothdata(PrDur_distr_cued_saline{n}, 'gaussian' , 5);
    PrDur_distr_uncued_saline{n} = histcounts(PrDur_uncued_saline{n}, edges, 'Normalization','probability');
    PrDur_distr_uncued_saline{n} = smoothdata(PrDur_distr_uncued_saline{n}, 'gaussian' , 5);
    
    PrDur_distr_cued_DCZ{n} = histcounts(PrDur_cued_DCZ{n}, edges, 'Normalization','probability');
    PrDur_distr_cued_DCZ{n} = smoothdata(PrDur_distr_cued_DCZ{n}, 'gaussian' , 5);
    PrDur_distr_uncued_DCZ{n} = histcounts(PrDur_uncued_DCZ{n}, edges, 'Normalization','probability');
    PrDur_distr_uncued_DCZ{n} = smoothdata(PrDur_distr_uncued_DCZ{n}, 'gaussian' , 5);
    
end
PrDur_distr_cued_Mat=[];
PrDur_distr_uncued_Mat=[];
AfterFP_distr_cued_Mat=[];
AfterFP_distr_uncued_Mat=[];
for n= 1:length(Animalfiles)
    PrDur_distr_cued_Mat       =   [PrDur_distr_cued_Mat;     PrDur_distr_cued_saline{n}; PrDur_distr_cued_DCZ{n};];
    PrDur_distr_uncued_Mat       =   [PrDur_distr_uncued_Mat;      PrDur_distr_uncued_saline{n}; PrDur_distr_uncued_DCZ{n};];
    AfterFP_distr_cued_Mat       =   [AfterFP_distr_cued_Mat;      AfterFP_distr_cued_saline{n}; AfterFP_distr_cued_DCZ{n}];
    AfterFP_distr_uncued_Mat       =   [AfterFP_distr_uncued_Mat;      AfterFP_distr_uncued_saline{n}; AfterFP_distr_uncued_DCZ{n}];
end

PrDur_median_cued_saline=[];
PrDur_median_uncued_saline=[];
PrDur_median_cued_DCZ=[];
PrDur_median_uncued_DCZ=[];
AfterFP_median_cued_saline=[];
AfterFP_median_uncued_saline=[];
AfterFP_median_cued_DCZ=[];
AfterFP_median_uncued_DCZ=[];
PrDur_cued_Ci_saline=[];
PrDur_uncued_Ci_saline=[];
PrDur_cued_Ci_DCZ=[];
PrDur_uncued_Ci_DCZ=[];
AfterFP_cued_Ci_saline=[];
AfterFP_uncued_Ci_saline=[];
AfterFP_cued_Ci_DCZ=[];
AfterFP_uncued_Ci_DCZ=[];
for n =1:length(Animalfiles)
    PrDur_median_cued_saline=[PrDur_median_cued_saline; median(PrDur_cued_saline{n})];
    PrDur_median_uncued_saline=[PrDur_median_uncued_saline; median(PrDur_uncued_saline{n})];
    PrDur_median_cued_DCZ=[ PrDur_median_cued_DCZ; median(PrDur_cued_DCZ{n})];
    PrDur_median_uncued_DCZ=[PrDur_median_uncued_DCZ; median(PrDur_uncued_DCZ{n})];
    AfterFP_median_cued_saline=[AfterFP_median_cued_saline; median(AfterFP_cued_saline{n})];
    AfterFP_median_uncued_saline=[AfterFP_median_uncued_saline; median(AfterFP_uncued_saline{n})];
    AfterFP_median_cued_DCZ=[AfterFP_median_cued_DCZ; median(AfterFP_cued_DCZ{n})];
    AfterFP_median_uncued_DCZ=[AfterFP_median_uncued_DCZ; median(AfterFP_uncued_DCZ{n})];
    PrDur_cued_Ci_saline(n,:)   = cell2mat(cellfun(@(x)prctile(bootstrp(1000, @median, x), [2.5, 97.5]),  PrDur_cued_saline(n), 'UniformOutput', false)');
    PrDur_uncued_Ci_saline(n,:)   = cell2mat(cellfun(@(x)prctile(bootstrp(1000, @median, x), [2.5, 97.5]),  PrDur_uncued_saline(n), 'UniformOutput', false)');
    PrDur_cued_Ci_DCZ(n,:)   = cell2mat(cellfun(@(x)prctile(bootstrp(1000, @median, x), [2.5, 97.5]),  PrDur_cued_DCZ(n), 'UniformOutput', false)');
    PrDur_uncued_Ci_DCZ(n,:)   = cell2mat(cellfun(@(x)prctile(bootstrp(1000, @median, x), [2.5, 97.5]),  PrDur_uncued_DCZ(n), 'UniformOutput', false)');
    AfterFP_cued_Ci_saline(n,:)   = cell2mat(cellfun(@(x)prctile(bootstrp(1000, @median, x), [2.5, 97.5]),  AfterFP_cued_saline(n), 'UniformOutput', false)');
    AfterFP_uncued_Ci_saline(n,:)   = cell2mat(cellfun(@(x)prctile(bootstrp(1000, @median, x), [2.5, 97.5]),  AfterFP_uncued_saline(n), 'UniformOutput', false)');
    AfterFP_cued_Ci_DCZ(n,:)   = cell2mat(cellfun(@(x)prctile(bootstrp(1000, @median, x), [2.5, 97.5]),  AfterFP_cued_DCZ(n), 'UniformOutput', false)');
    AfterFP_uncued_Ci_DCZ(n,:)   = cell2mat(cellfun(@(x)prctile(bootstrp(1000, @median, x), [2.5, 97.5]),  AfterFP_uncued_DCZ(n), 'UniformOutput', false)');
end
%% Plot
Salinecolor=[120 120 120]/255;
DCZcolor=[253,166,0]/255;

figure(815); clf;
set(gcf, 'unit', 'centimeters', 'position',[2 2 24 15], 'paperpositionmode', 'auto' );

% Press duration
ha1d=axes('unit', 'centimeters', 'position', [2 2 4 5], 'nextplot', 'add',...
    'xlim', [0 4000], 'ylim', [0.5 size(PrDur_distr_cued_Mat, 1)+0.5],'ytick', [0.5:size(PrDur_distr_cued_Mat, 1)+0.5]+0.5, 'yticklabel', {'Saline', 'DCZ'}, 'xtick', [0:1000:4000], 'ticklength', [.02 .025])

% do a little upsamnpling x10 
NewCenters_cued = linspace(edge_centers(1), edge_centers(end), 10*length(edge_centers));
NewDistribution_cued = interp1(edge_centers, PrDur_distr_cued_Mat', NewCenters_cued, 'v5cubic');
NewDistribution_cued = NewDistribution_cued';
himagePressDur = imagesc(NewCenters_cued, [1:size(NewDistribution_cued, 1)], NewDistribution_cued, [0 max(NewDistribution_cued(:))]);
colormap(ha1d, fake_parula);

for i = 1:length(Animalfiles)
    line(FP*[1 1], 0.5+[0 2]+2*(i-1), 'color', 'w', 'linewidth', 2, 'linestyle', ':');
    line([0 4000], 0.5+2*(i-1)*[1 1], 'color', 'r', 'linewidth', 1, 'linestyle', '-')
%     text(2300, 2*i-0.5, Animalnames{i}, 'fontsize', 10, 'color', 'w');
end
xlabel('Press duration (ms)');
title('Cued');

ha2d=axes('unit', 'centimeters', 'position', [6.5 2 4 5], 'nextplot', 'add',...
    'xlim', [0 4000], 'ylim', [0.5 size(PrDur_distr_uncued_Mat, 1)+0.5],'ytick', [], 'xtick', [0:1000:4000], 'ticklength', [.02 .025])

% do a little upsamnpling x10 
NewCenters_uncued = linspace(edge_centers(1), edge_centers(end), 10*length(edge_centers));
NewDistribution_uncued = interp1(edge_centers, PrDur_distr_uncued_Mat', NewCenters_uncued, 'v5cubic');
NewDistribution_uncued = NewDistribution_uncued';

himagePressDur = imagesc(NewCenters_uncued, [1:size(NewDistribution_uncued, 1)], NewDistribution_uncued, [0 max(NewDistribution_uncued(:))]);
colormap(ha2d, fake_parula)

for i = 1:length(Animalfiles)
    line(FP*[1 1], 0.5+[0 2]+2*(i-1), 'color', 'w', 'linewidth', 2, 'linestyle', ':');
    line([0 4000], 0.5+2*(i-1)*[1 1], 'color', 'r', 'linewidth', 1, 'linestyle', '-')
    text(2800, 2*i-0.5, Animalnames{i}, 'fontsize', 10, 'color', 'w');
end
xlabel('Press duration (ms)');
title('Uncued');

% RT include late
ha3d=axes('unit', 'centimeters', 'position', [13 2 4 5], 'nextplot', 'add',...
    'xlim', [0 3000], 'ylim', [0.5 size(AfterFP_distr_cued_Mat, 1)+0.5],'ytick', [0.5:size(AfterFP_distr_cued_Mat, 1)+0.5]+0.5, 'yticklabel', {'Saline', 'DCZ'}, 'xtick', [0:1000:3000], 'ticklength', [.02 .025])

% do a little upsamnpling x10 
NewCenters_cued = linspace(edge_centers(1), edge_centers(end), 10*length(edge_centers));
NewDistribution_cued = interp1(edge_centers, AfterFP_distr_cued_Mat', NewCenters_cued, 'v5cubic');
NewDistribution_cued = NewDistribution_cued';
himagePressDur = imagesc(NewCenters_cued, [1:size(NewDistribution_cued, 1)], NewDistribution_cued, [0 max(NewDistribution_cued(:))]);
colormap(ha3d, fake_parula);

for i = 1:length(Animalfiles)
    line([0 3000], 0.5+2*(i-1)*[1 1], 'color', 'r', 'linewidth', 1, 'linestyle', '-')
%     text(2300, 2*i-0.5, Animalnames{i}, 'fontsize', 10, 'color', 'w');
end
xlabel('RT include late (ms)');
title('Cued');

ha4d=axes('unit', 'centimeters', 'position', [17.5 2 4 5], 'nextplot', 'add',...
    'xlim', [0 3000], 'ylim', [0.5 size(AfterFP_distr_uncued_Mat, 1)+0.5],'ytick', [], 'xtick', [0:1000:3000], 'ticklength', [.02 .025])

% do a little upsamnpling x10 
NewCenters_uncued = linspace(edge_centers(1), edge_centers(end), 10*length(edge_centers));
NewDistribution_uncued = interp1(edge_centers, AfterFP_distr_uncued_Mat', NewCenters_uncued, 'v5cubic');
NewDistribution_uncued = NewDistribution_uncued';

himagePressDur = imagesc(NewCenters_uncued, [1:size(NewDistribution_uncued, 1)], NewDistribution_uncued, [0 max(NewDistribution_uncued(:))]);
colormap(ha4d, fake_parula)

for i = 1:length(Animalfiles)
    line([0 3000], 0.5+2*(i-1)*[1 1], 'color', 'r', 'linewidth', 1, 'linestyle', '-')
    text(2100, 2*i-0.5, Animalnames{i}, 'fontsize', 10, 'color', 'w');
end
xlabel('RT include late (ms)');
title('Uncued');

ha5d=axes('unit', 'centimeters', 'position', [2 9 4 4], 'nextplot', 'add',...
    'xlim', [0 3], 'ylim', [1250 2500],'ytick', [1000:250:2500], 'xtick', [1,2],'xticklabel',{'Saline','DCZ'}, 'ticklength', [.02 .025])
b = bar([1 2],[mean(PrDur_median_cued_saline),mean(PrDur_median_cued_DCZ)],'FaceColor','flat','LineStyle','none');
b.CData(1,:) = Salinecolor;
b.CData(2,:) = DCZcolor;
for n = 1:length(Animalfiles)
    errorbar(1,PrDur_median_cued_saline(n),abs(PrDur_cued_Ci_saline(n,1)-PrDur_median_cued_saline(n)),abs(PrDur_cued_Ci_saline(n,2)-PrDur_median_cued_saline(n)),'k.','MarkerSize',12,'lineWidth',1,'CapSize',3);
    errorbar(2,PrDur_median_cued_DCZ(n),abs(PrDur_cued_Ci_DCZ(n,1)-PrDur_median_cued_DCZ(n)),abs(PrDur_cued_Ci_DCZ(n,2)-PrDur_median_cued_DCZ(n)),'k.','MarkerSize',12,'lineWidth',1,'CapSize',3);
    line([1 2],[PrDur_median_cued_saline(n) PrDur_median_cued_DCZ(n)],'color','k')
end
[p, h]=ranksum(PrDur_median_cued_saline, PrDur_median_cued_DCZ);
text(1, 1350, sprintf('p=%2.5f', p))
ylabel('Press duration (ms)');
title('Cued');

ha6d=axes('unit', 'centimeters', 'position', [6.3 9 4 4], 'nextplot', 'add',...
    'xlim', [0 3], 'ylim', [1250 2500],'ytick', [1000:250:2500],'yticklabel',{}, 'xtick', [1,2],'xticklabel',{'Saline','DCZ'}, 'ticklength', [.02 .025])
b = bar([1 2],[mean(PrDur_median_uncued_saline),mean(PrDur_median_uncued_DCZ)],'FaceColor','flat','LineStyle','none');
b.CData(1,:) = Salinecolor;
b.CData(2,:) = DCZcolor;
for n = 1:length(Animalfiles)
    errorbar(1,PrDur_median_uncued_saline(n),abs(PrDur_uncued_Ci_saline(n,1)-PrDur_median_uncued_saline(n)),abs(PrDur_uncued_Ci_saline(n,2)-PrDur_median_uncued_saline(n)),'k.','MarkerSize',12,'lineWidth',1,'CapSize',3);
    errorbar(2,PrDur_median_uncued_DCZ(n),abs(PrDur_uncued_Ci_DCZ(n,1)-PrDur_median_uncued_DCZ(n)),abs(PrDur_uncued_Ci_DCZ(n,2)-PrDur_median_uncued_DCZ(n)),'k.','MarkerSize',12,'lineWidth',1,'CapSize',3);
    line([1 2],[PrDur_median_uncued_saline(n) PrDur_median_uncued_DCZ(n)],'color','k')
end
[p, h]=ranksum(PrDur_median_uncued_saline, PrDur_median_uncued_DCZ);
text(1, 1350, sprintf('p=%2.5f', p))
title('Uncued');

ha7d=axes('unit', 'centimeters', 'position', [13 9 4 4], 'nextplot', 'add',...
    'xlim', [0 3], 'ylim', [0 1250],'ytick', [0:250:1250], 'xtick', [1,2],'xticklabel',{'Saline','DCZ'}, 'ticklength', [.02 .025])
b = bar([1 2],[mean(AfterFP_median_cued_saline),mean(AfterFP_median_cued_DCZ)],'FaceColor','flat','LineStyle','none');
b.CData(1,:) = Salinecolor;
b.CData(2,:) = DCZcolor;
for n = 1:length(Animalfiles)
    errorbar(1,AfterFP_median_cued_saline(n),abs(AfterFP_cued_Ci_saline(n,1)-AfterFP_median_cued_saline(n)),abs(AfterFP_cued_Ci_saline(n,2)-AfterFP_median_cued_saline(n)),'k.','MarkerSize',12,'lineWidth',1,'CapSize',3);
    errorbar(2,AfterFP_median_cued_DCZ(n),abs(AfterFP_cued_Ci_DCZ(n,1)-AfterFP_median_cued_DCZ(n)),abs(AfterFP_cued_Ci_DCZ(n,2)-AfterFP_median_cued_DCZ(n)),'k.','MarkerSize',12,'lineWidth',1,'CapSize',3);
    line([1 2],[AfterFP_median_cued_saline(n) AfterFP_median_cued_DCZ(n)],'color','k')
end
[p, h]=ranksum(AfterFP_median_cued_saline, AfterFP_median_cued_DCZ);
text(1, 100, sprintf('p=%2.5f', p))
ylabel('RT include late (ms)');
title('Cued');

ha8d=axes('unit', 'centimeters', 'position', [17.3 9 4 4], 'nextplot', 'add',...
    'xlim', [0 3], 'ylim', [0 1250],'ytick', [0:250:1250],'yticklabel',{}, 'xtick', [1,2],'xticklabel',{'Saline','DCZ'}, 'ticklength', [.02 .025])
b = bar([1 2],[mean(AfterFP_median_uncued_saline),mean(AfterFP_median_uncued_DCZ)],'FaceColor','flat','LineStyle','none');
b.CData(1,:) = Salinecolor;
b.CData(2,:) = DCZcolor;
for n = 1:length(Animalfiles)
    errorbar(1,AfterFP_median_uncued_saline(n),abs(AfterFP_uncued_Ci_saline(n,1)-AfterFP_median_uncued_saline(n)),abs(AfterFP_uncued_Ci_saline(n,2)-AfterFP_median_uncued_saline(n)),'k.','MarkerSize',12,'lineWidth',1,'CapSize',3);
    errorbar(2,AfterFP_median_uncued_DCZ(n),abs(AfterFP_uncued_Ci_DCZ(n,1)-AfterFP_median_uncued_DCZ(n)),abs(AfterFP_uncued_Ci_DCZ(n,2)-AfterFP_median_uncued_DCZ(n)),'k.','MarkerSize',12,'lineWidth',1,'CapSize',3);
    line([1 2],[AfterFP_median_uncued_saline(n) AfterFP_median_uncued_DCZ(n)],'color','k')
end
[p, h]=ranksum(AfterFP_median_uncued_saline, AfterFP_median_uncued_DCZ);
text(1, 100, sprintf('p=%2.5f', p))
title('Uncued');
%% Save
savename = ['Kornblum_Heatmap_Group'];
mkdir('Analysis');
savename = fullfile(pwd, 'Analysis', savename)                                                           
saveas(gcf, savename, 'png')