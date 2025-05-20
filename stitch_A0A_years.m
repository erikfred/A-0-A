% stitch_A0A_years.m
%
% starting point for attempting to create multi-year pressure time series
% from discrete deployments
%

clear; close all

% name mapping between years
name1={'POBS-01','POBS-03','POBS-15','POBS-12','POBS-05','POBS-07','POBS-11',...
    'POBS-04','POBS-16','POBS-02'};
name2={'POBS-08','POBS-02','POBS-12','POBS-05','POBS-14','POBS-13','POBS-16',...
    'POBS-03','POBS-11','POBS-01'};

g1=load('../pressure_data/geometry.mat');
g2=load('../pressure_data_Y2/geometry.mat');

% check actual distance between 'repeat' sites
for i=1:length(name1)
    i1=find(strcmp(name1{i},g1.staname));
    lat1=g1.stalat(i1);
    lon1=g1.stalon(i1);
    dep1=g1.stadepth(i1);

    i2=find(strcmp(name2{i},g2.staname));
    lat2=g2.stalat(i2);
    lon2=g2.stalon(i2);
    dep2=g2.stadepth(i2);

    disp([name1{i} ' vs. ' name2{i} ':'])
    disp(['Separation = ' num2str(round(lldistkm([lat1 lon1],[lat2 lon2])*1000)) ' m'])
    disp(['Depth difference = ' num2str(abs(dep1-dep2)) ' m'])
    disp(newline)
end

% what do the data actually look like side-by-side?
for i=1:length(name1)
    i1=find(strcmp(name1{i},g1.staname));
    i2=find(strcmp(name2{i},g2.staname));

    if length(g1.filename{i1})>10
        Y1=load(['../stitched_data/drift_corrected/' g1.filename{i1}(1:5) '.mat'],'dataf');
    else
        Y1=load(['../stitched_data/drift_corrected/' g1.filename{i1} '.mat'],'dataf');
    end
    Y2=load(['../stitched_data_Y2/drift_corrected/' g2.filename{i2} '.mat'],'dataf');

    ta=Y1.dataf.tf;
    pa1=Y1.dataf.p1_dcor;
    pa2=Y1.dataf.p2_dcor;

    tb=Y2.dataf.tf;
    pb1=Y2.dataf.p1_dcor;
    pb2=Y2.dataf.p2_dcor;

    figure(55); clf
    subplot(211); hold on
    plot(ta,pa1-mean(pa1),'linewidth',1)
    plot(tb,pb1-mean(pb1),'linewidth',1)
    datetick('x',3)
    box on; grid on
    subplot(212); hold on
    plot(ta,pa2-mean(pa2),'linewidth',1)
    plot(tb,pb2-mean(pb2),'linewidth',1)
    datetick('x',3)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print(['../figures_Y2/combine_years/' name1{i} '&' name2{i}],'-dpng','-r300')

end