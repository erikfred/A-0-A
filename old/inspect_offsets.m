% inspect_offsets.m

clear; close all

POBS_dir=dir('../pressure_data');
POBS_list={POBS_dir.name}';
file_check=cellfun(@(v)v(1),POBS_list);
i_list=find(eq(file_check,'P')); % subselects only files begining with 'P'
i_list([4,5,13])=[];

% zoom in on specific anomalous offsets
for i=1:length(i_list)
    k=i_list(i);

    load([POBS_dir(k).folder '/' POBS_dir(k).name])

    t=data.ta(data.lNormState1);
    p1=data.a1(data.lNormState1);
    p2=data.a2(data.lNormState1);
    T1=data.Ta1(data.lNormState1);
    T2=data.Ta2(data.lNormState1);

    t_0=datenum(2023,10,19,16,10,0);
    t_f=datenum(2023,10,19,16,30,0);
    keyboard % change t-limits to desired range

    % trim data
    it=t>=t_0 & t<=t_f;
    tt=t(it); pp1=p1(it); pp2=p2(it);

    figure(20); clf; hold on
    plot(tt,detrend(pp1)+0.5,'linewidth',1)
    plot(tt,detrend(pp2)-0.5,'linewidth',1)
    plot(tt,(pp1-mean(pp1))-(pp2-mean(pp2)),'k','linewidth',2)
    ylabel('P (hPa)')
    legend('Gauge 1','Gauge 2',['Gauge 1' char(8212) 'Gauge 2'],'location','northwest')
    datetick('x')
    set(gca,'fontsize',14)
    box on; grid on

    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../figures/spurious_offsets/POBS4/Oct19','-dpng','-r300')
end