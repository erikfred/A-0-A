% best_gauge_Y2.m
%
% For gauges that do not agree, make plots/comparisons that help logically
% justify which of the two is to be trusted.
%

clear; close all

G=load('../pressure_data_Y2/geometry.mat');

badname={'POBS01';'POBS03';'POBS04'};

POBS_dir=dir('../stitched_data_Y2/drift_corrected');
POBS_list={POBS_dir.name}';
file_check=cellfun(@(v)v(1),POBS_list);
i_list=find(eq(file_check,'P')); % subselects only files begining with 'P'
POBS_files=POBS_list(i_list);

for i=1:length(badname)
    % find and load corresponding matfile
    ibad=find(strcmp([badname{i} '.mat'],POBS_files));
    load([POBS_dir(1).folder '/' POBS_files{ibad}])
    Pbad=dataf;

    % extract relevant geometry info
    ig=find(strcmp(badname{i},G.filename));
    dbad=G.stadepth(ig);

    % find nearest-depth station(s)
    [dtest,idt]=sort(abs(dbad-G.stadepth));
    n=2;
    id=idt(n); % this will be the nearest-depth other than self

    % compare against 4 nearest-depth stations, excluding other bads
    count=0;
    while count<4
        % check for other bad stations and skip
        if any(strcmp(G.filename{id},badname))
            n=n+1;
            id=idt(n);
            continue
        end

        % load comp station's pressure data
        load([POBS_dir(1).folder '/' G.filename{id} '.mat'])
        Pgood=dataf;

        % extract temporal overlap
        [tplot,ia,ib]=intersect(round(Pbad.tf,6),round(Pgood.tf,6));
        pb1=Pbad.p1_dcor(ia);
        pb2=Pbad.p2_dcor(ia);
        pg1=Pgood.p1_dcor(ib);
        pg2=Pgood.p2_dcor(ib);

        % deal with NaNs as necessary
        if any(isnan(pg1))
            pg1(isnan(pg1))=pg2(isnan(pg1))-pg2(1)+pg1(1);
        elseif any(isnan(pg2))
            pg2(isnan(pg2))=pg1(isnan(pg2))-pg1(1)+pg2(1);
        end

        % calculate trend lines for each difference
        pdif11=(pb1-mean(pb1))-(pg1-mean(pg1));
        m11=polyfit((1:length(tplot))',pdif11,1);
        fit11=polyval(m11,(1:length(tplot))');
        pdif12=(pb1-mean(pb1))-(pg2-mean(pg2));
        m12=polyfit((1:length(tplot))',pdif12,1);
        fit12=polyval(m12,(1:length(tplot))');
        pdif21=(pb2-mean(pb2))-(pg1-mean(pg1));
        m21=polyfit((1:length(tplot))',pdif21,1);
        fit21=polyval(m21,(1:length(tplot))');
        pdif22=(pb2-mean(pb2))-(pg2-mean(pg2));
        m22=polyfit((1:length(tplot))',pdif22,1);
        fit22=polyval(m22,(1:length(tplot))');

        % some display parameters
        lx=round(length(tplot)*2/3);
        ly=-2*std(pb1);

        % 4 plots for the 4 gauge combinations
        figure(4); clf
        % 1 v 1
        subplot(221); hold on
        plot(tplot,pb1-mean(pb1),'linewidth',1)
        plot(tplot,pg1-mean(pg1),'linewidth',1)
        plot(tplot,pdif11,'k','linewidth',2)
        plot(tplot,fit11,'k','linewidth',2)
        text(tplot(lx),ly,[num2str(round(m11(1)*24*365,2)) ' hPa/yr'],'fontsize',12)
        datetick('x',3,'keeplimits')
        title('BAD 1 - GOOD 1')
        box on; grid on
        % 1 v 2
        subplot(222); hold on
        plot(tplot,pb1-mean(pb1),'linewidth',1)
        plot(tplot,pg2-mean(pg2),'linewidth',1)
        plot(tplot,pdif12,'k','linewidth',2)
        plot(tplot,fit12,'k','linewidth',2)
        text(tplot(lx),ly,[num2str(round(m12(1)*24*365,2)) ' hPa/yr'],'fontsize',12)
        datetick('x',3,'keeplimits')
        title('BAD 1 - GOOD 2')
        box on; grid on
        % 2 v 1
        subplot(223); hold on
        plot(tplot,pb2-mean(pb2),'linewidth',1)
        plot(tplot,pg1-mean(pg1),'linewidth',1)
        plot(tplot,pdif21,'k','linewidth',2)
        plot(tplot,fit21,'k','linewidth',2)
        text(tplot(lx),ly,[num2str(round(m21(1)*24*365,2)) ' hPa/yr'],'fontsize',12)
        datetick('x',3,'keeplimits')
        title('BAD 2 - GOOD 1')
        box on; grid on
        % 2 v 2
        subplot(224); hold on
        plot(tplot,pb2-mean(pb2),'linewidth',1)
        plot(tplot,pg2-mean(pg2),'linewidth',1)
        plot(tplot,pdif22,'k','linewidth',2)
        plot(tplot,fit22,'k','linewidth',2)
        text(tplot(lx),ly,[num2str(round(m22(1)*24*365,2)) ' hPa/yr'],'fontsize',12)
        datetick('x',3,'keeplimits')
        title('BAD 2 - GOOD 2')
        box on; grid on

        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];
        print(['../figures_Y2/gauge_assessments/' badname{i} '_vs_' G.filename{id}],'-dpng','-r300')

        n=n+1; id=idt(n);
        count=count+1;
    end
end