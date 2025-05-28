% stitch_A0A_years.m
%
% starting point for attempting to create multi-year pressure time series
% from discrete deployments
%

clear; close all

meancomp=false; % assume two parts have same mean
pobs9comp=false; % assume alignment with POBS-09
glocomp=false; % assume alignment with GLORYS output
cpiescomp=true; % assume alignment with nearest CPIES

cm=lines;

% name mapping between years
name1={'POBS-04','POBS-03','POBS-02','POBS-01','POBS-07','POBS-15','POBS-11',...
    'POBS-12','POBS-05','POBS-16'};
name2={'POBS-03','POBS-02','POBS-01','POBS-08','POBS-13','POBS-12','POBS-16',...
    'POBS-05','POBS-14','POBS-11'};

g1=load('../pressure_data/geometry.mat');
g2=load('../pressure_data_Y2/geometry.mat');

if pobs9comp
    p09=load('../stitched_data/drift_corrected/POBS09.mat','dataf');
end
if glocomp
    load('../../hikurangi/ocean_data/GLORYS/glorys_pressure_Y2','p_glo')
end
if cpiescomp
    load('../CPIES_telemetry_2024/CPIES_formatted.mat');
end

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

    % remove first 2 weeks, when calibrations aren't great
    ta=ta(24*14:end); pa1=pa1(24*14:end); pa2=pa2(24*14:end);
    tb=tb(24*14:end); pb1=pb1(24*14:end); pb2=pb2(24*14:end);

    %--- de-mean and plot side by side
    if meancomp
        % individual plots
        figure(55); clf
        subplot(211); hold on
        plot(ta,pa1-mean(pa1),'linewidth',1)
        plot(tb,pb1-mean(pb1),'linewidth',1)
        title('Gauge 1')
        datetick('x',3)
        box on; grid on
        subplot(212); hold on
        plot(ta,pa2-mean(pa2),'linewidth',1)
        plot(tb,pb2-mean(pb2),'linewidth',1)
        title('Gauge 2')
        datetick('x',3)
        box on; grid on

        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];
        % print(['../figures_Y2/combine_years/demean/' name1{i} '&' name2{i}],'-dpng','-r300')

        % Gauge 1 stack
        figure(56); hold on
        plot(ta,pa1-mean(pa1)+(i-1)*10,'color',cm(i,:),'linewidth',1)
        text(ta(1)-30,pa1(1)-mean(pa1)+(i-1)*10,name1{i})
        plot(tb,pb1-mean(pb1)+(i-1)*10,'color',cm(i,:),'linewidth',1)
        text(tb(end)+10,pb1(end)-mean(pb1)+(i-1)*10,name2{i})
        datetick('x',3)
        box on; grid on

        % Gauge 2 stack
        figure(57); hold on
        plot(ta,pa2-mean(pa2)+(i-1)*10,'color',cm(i,:),'linewidth',1)
        text(ta(1)-30,pa2(1)-mean(pa2)+(i-1)*10,name1{i})
        plot(tb,pb2-mean(pb2)+(i-1)*10,'color',cm(i,:),'linewidth',1)
        text(tb(end)+10,pb2(end)-mean(pb2)+(i-1)*10,name2{i})
        datetick('x',3)
        box on; grid on
    end

    %--- align to POBS-09 (spans two deployments)
    if pobs9comp
        tc=p09.dataf.tf;
        pc1=p09.dataf.p1_dcor; pc1=pc1-mean(pc1);
        pc2=p09.dataf.p2_dcor; pc2=pc2-mean(pc2);

        % find early overlap
        [tearly,ia,ib]=intersect(round(ta,6),round(tc,6));
        ppa1=pa1(ia); ppc1=pc1(ib);
        ppa2=pa2(ia); ppc2=pc2(ib);

        % trim to remove first 2 weeks of drift
        tearly=tearly(24*14:end);
        ppa1=ppa1(24*14:end); ppc1=ppc1(24*14:end);
        ppa2=ppa2(24*14:end); ppc2=ppc2(24*14:end);

        % mean of difference is optimal offset between them
        ofst_a1=mean(ppa1-ppc1);
        ofst_a2=mean(ppa2-ppc2);

        % find late overlap
        [tlate,ia,ib]=intersect(round(tb,6),round(tc,6));
        ppb1=pb1(ia); ppc1=pc1(ib);
        ppb2=pb2(ia); ppc2=pc2(ib);

        % do not need to trim late segment for drift

        % mean of difference is optimal offset between them
        ofst_b1=mean(ppb1-ppc1);
        ofst_b2=mean(ppb2-ppc2);

        % find linear fit to combined time series
        tinv=([ta;tb]-ta(1))/tb(end);
        pinv1=[pa1-ofst_a1;pb1-ofst_b1];
        pinv2=[pa2-ofst_a2;pb2-ofst_b2];
        G=[tinv,ones(size(tinv))];
        m1=inv(G'*G)*G'*pinv1; plin1=G*m1;
        m2=inv(G'*G)*G'*pinv2; plin2=G*m2;
        p9rate1(i)=(plin1(end)-plin1(1))/(tb(end)-ta(1))*365;
        p9rate2(i)=(plin2(end)-plin2(1))/(tb(end)-ta(1))*365;

        % individual plots
        figure(25); clf;
        subplot(211); hold on
        plot(tc,pc1,'k','linewidth',1)
        plot(ta,pa1-ofst_a1,'linewidth',2)
        text(ta(1)-30,pa1(1)-ofst_a1,name1{i})
        plot(tb,pb1-ofst_b1,'linewidth',2)
        text(tb(end)+10,pb1(end)-ofst_b1,name2{i})
        plot([ta;tb],plin1,'k')
        title('Gauge 1')
        datetick('x',3)
        box on; grid on
        subplot(212); hold on
        plot(tc,pc2,'k','linewidth',1)
        plot(ta,pa2-ofst_a2,'linewidth',2)
        text(ta(1)-30,pa2(1)-ofst_a2,name1{i})
        plot(tb,pb2-ofst_b2,'linewidth',2)
        text(tb(end)+10,pb2(end)-ofst_b2,name2{i})
        plot([ta;tb],plin2,'k')
        title('Gauge 2')
        datetick('x',3)
        box on; grid on

        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];
        print(['../figures_Y2/combine_years/pobs09/' name1{i} '&' name2{i}],'-dpng','-r300')
        
        % Gauge 1 stack
        figure(26); hold on
        plot(tc,pc1+(i-1)*10,'k','linewidth',1)
        plot(ta,pa1-ofst_a1+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(ta(1)-30,pa1(1)-ofst_a1+(i-1)*10,name1{i})
        plot(tb,pb1-ofst_b1+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(tb(end)+10,pb1(end)-ofst_b1+(i-1)*10,name2{i})
        plot([ta;tb],plin1+(i-1)*10,'--k')
        datetick('x',3)
        box on; grid on

        % Gauge 2 stack
        figure(27); hold on
        plot(tc,pc2+(i-1)*10,'k','linewidth',1)
        plot(ta,pa2-ofst_a2+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(ta(1)-30,pa2(1)-ofst_a2+(i-1)*10,name1{i})
        plot(tb,pb2-ofst_b2+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(tb(end)+10,pb2(end)-ofst_b2+(i-1)*10,name2{i})
        plot([ta;tb],plin2+(i-1)*10,'--k')
        datetick('x',3)
        box on; grid on
    end

    %--- align to GLORYS
    if glocomp
        i3=find(strcmp(name1{i},p_glo.nearest_sta));

        tc=(p_glo.t{i3}(1):1/24:p_glo.t{i3}(end))';
        pc=interp1(p_glo.t{i3},p_glo.bpa{i3}',tc); pc=pc-mean(pc);

        % find early overlap
        [tearly,ia,ib]=intersect(round(ta,6),round(tc,6));
        ppa1=pa1(ia); ppa2=pa2(ia); ppc=pc(ib);

        % trim to remove first 2 weeks of drift
        ppa1=ppa1(24*14:end); ppa2=ppa2(24*14:end); ppc=ppc(24*14:end);

        % mean of difference is optimal offset between them
        ofst_a1=mean(ppa1-ppc);
        ofst_a2=mean(ppa2-ppc);

        % find late overlap
        [tlate,ia,ib]=intersect(round(tb,6),round(tc,6));
        ppb1=pb1(ia); ppb2=pb2(ia); ppc=pc(ib);

        % do not need to trim late segment for drift

        % mean of difference is optimal offset between them
        ofst_b1=mean(ppb1-ppc);
        ofst_b2=mean(ppb2-ppc);

        % find linear fit to combined time series
        tinv=([ta;tb]-ta(1))/tb(end);
        pinv1=[pa1-ofst_a1;pb1-ofst_b1];
        pinv2=[pa2-ofst_a2;pb2-ofst_b2];
        G=[tinv,ones(size(tinv))];
        m1=inv(G'*G)*G'*pinv1; plin1=G*m1;
        m2=inv(G'*G)*G'*pinv2; plin2=G*m2;
        glorate1(i)=(plin1(end)-plin1(1))/(tb(end)-ta(1))*365;
        glorate2(i)=(plin2(end)-plin2(1))/(tb(end)-ta(1))*365;

        % individual plots
        figure(15); clf;
        subplot(211); hold on
        plot(tc,pc,'k','linewidth',1)
        plot(ta,pa1-ofst_a1,'linewidth',2)
        text(ta(1)-30,pa1(1)-ofst_a1,name1{i})
        plot(tb,pb1-ofst_b1,'linewidth',2)
        text(tb(end)+10,pb1(end)-ofst_b1,name2{i})
        plot([ta;tb],plin1,'k')
        title('Gauge 1')
        datetick('x',3)
        box on; grid on
        subplot(212); hold on
        plot(tc,pc,'k','linewidth',1)
        plot(ta,pa2-ofst_a2,'linewidth',2)
        text(ta(1)-30,pa2(1)-ofst_a2,name1{i})
        plot(tb,pb2-ofst_b2,'linewidth',2)
        text(tb(end)+10,pb2(end)-ofst_b2,name2{i})
        plot([ta;tb],plin2,'k')
        title('Gauge 2')
        datetick('x',3)
        box on; grid on

        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];
        print(['../figures_Y2/combine_years/glorys/' name1{i} '&' name2{i}],'-dpng','-r300')
        
        % Gauge 1 stack
        figure(16); hold on
        plot(tc,pc+(i-1)*10,'k','linewidth',1)
        plot(ta,pa1-ofst_a1+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(ta(1)-30,pa1(1)-ofst_a1+(i-1)*10,name1{i})
        plot(tb,pb1-ofst_b1+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(tb(end)+10,pb1(end)-ofst_b1+(i-1)*10,name2{i})
        plot([ta;tb],plin1+(i-1)*10,'k')
        datetick('x',3)
        box on; grid on

        % Gauge 2 stack
        figure(17); hold on
        plot(tc,pc+(i-1)*10,'k','linewidth',1)
        plot(ta,pa2-ofst_a2+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(ta(1)-30,pa2(1)-ofst_a2+(i-1)*10,name1{i})
        plot(tb,pb2-ofst_b2+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(tb(end)+10,pb2(end)-ofst_b2+(i-1)*10,name2{i})
        plot([ta;tb],plin2+(i-1)*10,'k')
        datetick('x',3)
        box on; grid on

        % what does it look like if I subtract out the model?
        tcomb=[ta;tb];
        pcomb1=[pa1-ofst_a1;pb1-ofst_b1];
        pcomb2=[pa2-ofst_a2;pb2-ofst_b2];
        [tdif,ia,ib]=intersect(round(tc,6),round(tcomb,6));
        pcor1=pcomb1(ib)-pc(ia);
        pcor2=pcomb2(ib)-pc(ia);

        % find linear fit to corrected time series
        tinv=(tdif-tdif(1))/tdif(end);
        G=[tinv,ones(size(tinv))];
        m1=inv(G'*G)*G'*pcor1; plin1=G*m1;
        m2=inv(G'*G)*G'*pcor2; plin2=G*m2;
        gcorrate1(i)=(plin1(end)-plin1(1))/(tdif(end)-tdif(1))*365;
        gcorrate2(i)=(plin2(end)-plin2(1))/(tdif(end)-tdif(1))*365;

        ijump=find(diff(tdif)>0.05);

        % Gauge 1 stack
        figure(18); hold on
        plot(tdif(1:ijump),pcor1(1:ijump)+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(tdif(1)-30,pcor1(1)+(i-1)*10,name1{i})
        plot(tdif(ijump+1:end),pcor1(ijump+1:end)+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(tdif(end)+10,pcor1(end)+(i-1)*10,name2{i})
        plot(tdif,plin1+(i-1)*10,'k','linewidth',1)
        datetick('x',3)
        box on; grid on

        % Gauge 2 stack
        figure(19); hold on
        plot(tdif(1:ijump),pcor2(1:ijump)+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(tdif(1)-30,pcor2(1)+(i-1)*10,name1{i})
        plot(tdif(ijump+1:end),pcor2(ijump+1:end)+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(tdif(end)+10,pcor2(end)+(i-1)*10,name2{i})
        plot(tdif,plin2+(i-1)*10,'k','linewidth',1)
        datetick('x',3)
        box on; grid on
    end

    %--- align to CPIES
    if cpiescomp
        % find depth-matched CPIES station
        itemp=find(strcmp(name1{i},g1.staname));
        [~,i3]=min(abs(cpies.z-g1.stadepth(itemp)));

        if strcmp('POBS-07',g1.staname{itemp})
            i3=1;
        end

        tc=(cpies.t{i3}(1):1/24:cpies.t{i3}(end))';
        pc=interp1(cpies.t{i3},cpies.bpa{i3},tc); pc=pc-mean(pc);

        % find early overlap
        [tearly,ia,ib]=intersect(round(ta,6),round(tc,6));
        ppa1=pa1(ia); ppa2=pa2(ia); ppc=pc(ib);

        % trim to remove first 2 weeks of drift
        ppa1=ppa1(24*14:end); ppa2=ppa2(24*14:end); ppc=ppc(24*14:end);

        % mean of difference is optimal offset between them
        ofst_a1=mean(ppa1-ppc);
        ofst_a2=mean(ppa2-ppc);

        % find late overlap
        [tlate,ia,ib]=intersect(round(tb,6),round(tc,6));
        ppb1=pb1(ia); ppb2=pb2(ia); ppc=pc(ib);

        % do not need to trim late segment for drift

        % mean of difference is optimal offset between them
        ofst_b1=mean(ppb1-ppc);
        ofst_b2=mean(ppb2-ppc);

        % find linear fit to combined time series
        tinv=([ta;tb]-ta(1))/tb(end);
        pinv1=[pa1-ofst_a1;pb1-ofst_b1];
        pinv2=[pa2-ofst_a2;pb2-ofst_b2];
        G=[tinv,ones(size(tinv))];
        m1=inv(G'*G)*G'*pinv1; plin1=G*m1;
        m2=inv(G'*G)*G'*pinv2; plin2=G*m2;
        cpierate1(i)=(plin1(end)-plin1(1))/(tb(end)-ta(1))*365;
        cpierate2(i)=(plin2(end)-plin2(1))/(tb(end)-ta(1))*365;

        % individual plots
        figure(35); clf;
        subplot(211); hold on
        plot(tc,pc,'k','linewidth',1)
        plot(ta,pa1-ofst_a1,'linewidth',2)
        text(ta(1)-30,pa1(1)-ofst_a1,name1{i})
        plot(tb,pb1-ofst_b1,'linewidth',2)
        text(tb(end)+10,pb1(end)-ofst_b1,name2{i})
        plot([ta;tb],plin1,'k')
        title('Gauge 1')
        datetick('x',3)
        box on; grid on
        subplot(212); hold on
        plot(tc,pc,'k','linewidth',1)
        plot(ta,pa2-ofst_a2,'linewidth',2)
        text(ta(1)-30,pa2(1)-ofst_a2,name1{i})
        plot(tb,pb2-ofst_b2,'linewidth',2)
        text(tb(end)+10,pb2(end)-ofst_b2,name2{i})
        plot([ta;tb],plin2,'k')
        title('Gauge 2')
        datetick('x',3)
        box on; grid on

        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];
        print(['../figures_Y2/combine_years/cpies/' name1{i} '&' name2{i}],'-dpng','-r300')
        
        % Gauge 1 stack
        figure(36); hold on
        plot(tc,pc+(i-1)*10,'k','linewidth',1)
        plot(ta,pa1-ofst_a1+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(ta(1)-30,pa1(1)-ofst_a1+(i-1)*10,name1{i})
        plot(tb,pb1-ofst_b1+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(tb(end)+10,pb1(end)-ofst_b1+(i-1)*10,name2{i})
        plot([ta;tb],plin1+(i-1)*10,'k')
        datetick('x',3)
        box on; grid on

        % Gauge 2 stack
        figure(37); hold on
        plot(tc,pc+(i-1)*10,'k','linewidth',1)
        plot(ta,pa2-ofst_a2+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(ta(1)-30,pa2(1)-ofst_a2+(i-1)*10,name1{i})
        plot(tb,pb2-ofst_b2+(i-1)*10,'color',cm(i,:),'linewidth',2)
        text(tb(end)+10,pb2(end)-ofst_b2+(i-1)*10,name2{i})
        plot([ta;tb],plin2+(i-1)*10,'k')
        datetick('x',3)
        box on; grid on
    end
end

if meancomp
    figure(56)
    ylim([-20 120])
    title('Gauge 1')
    set(gca,'fontsize',14)
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../figures_Y2/combine_years/demean/stack_G1','-dpng','-r300')
    print('../figures_Y2/combine_years/demean/stack_G1','-depsc','-vector')

    figure(57)
    ylim([-20 120])
    title('Gauge 2')
    set(gca,'fontsize',14)
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../figures_Y2/combine_years/demean/stack_G2','-dpng','-r300')
    print('../figures_Y2/combine_years/demean/stack_G2','-depsc','-vector')
end

if pobs9comp
    figure(26)
    ylim([-20 120])
    title('Gauge 1')
    set(gca,'fontsize',14)
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../figures_Y2/combine_years/pobs09/stack_G1','-dpng','-r300')
    print('../figures_Y2/combine_years/pobs09/stack_G1','-depsc','-vector')

    figure(27)
    ylim([-20 120])
    title('Gauge 2')
    set(gca,'fontsize',14)
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../figures_Y2/combine_years/pobs09/stack_G2','-dpng','-r300')
    print('../figures_Y2/combine_years/pobs09/stack_G2','-depsc','-vector')
end

if glocomp
    figure(16)
    ylim([-20 120])
    title('Gauge 1')
    set(gca,'fontsize',14)
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../figures_Y2/combine_years/glorys/stack_G1','-dpng','-r300')
    print('../figures_Y2/combine_years/glorys/stack_G1','-depsc','-vector')

    figure(17)
    ylim([-20 120])
    title('Gauge 2')
    set(gca,'fontsize',14)
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../figures_Y2/combine_years/glorys/stack_G2','-dpng','-r300')
    print('../figures_Y2/combine_years/glorys/stack_G2','-depsc','-vector')

    figure(18)
    ylim([-20 120])
    title('Gauge 1')
    set(gca,'fontsize',14)
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../figures_Y2/combine_years/glorys/corstack_G1','-dpng','-r300')
    print('../figures_Y2/combine_years/glorys/corstack_G1','-depsc','-vector')

    figure(19)
    ylim([-20 120])
    title('Gauge 2')
    set(gca,'fontsize',14)
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../figures_Y2/combine_years/glorys/corstack_G2','-dpng','-r300')
    print('../figures_Y2/combine_years/glorys/corstack_G2','-depsc','-vector')
end

if cpiescomp
    figure(36)
    ylim([-20 120])
    title('Gauge 1')
    set(gca,'fontsize',14)
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../figures_Y2/combine_years/cpies/stack_G1','-dpng','-r300')
    print('../figures_Y2/combine_years/cpies/stack_G1','-depsc','-vector')

    figure(37)
    ylim([-20 120])
    title('Gauge 2')
    set(gca,'fontsize',14)
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../figures_Y2/combine_years/cpies/stack_G2','-dpng','-r300')
    print('../figures_Y2/combine_years/cpies/stack_G2','-depsc','-vector')
end