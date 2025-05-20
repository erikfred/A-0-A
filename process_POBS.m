% process_POBS.m
%
% Script to process 2022-2023 Hikurangi A-0-A data
%

clear; close all

%% Calibration parameters

% Temperature sensitivity parameters (yet to be determined)
p.dadT1 = 0;
p.TRef1 = 0;
p.dadT2 = 0;
p.TRef2 = 0;
p.dadTb = 0;
p.TRefb = 0;

% Find calibration parameters
p.pThresh = 1100;               % 1100 hPa threshold for a calibration (housing pressure ~= 1 atm)
p.minTime4Cal = 290;            % Minimum duration in seconds for a calibration to be counted
p.complexRange80 = 100;           % If 90% value - 10% value is larger than this then not a simple calibration
p.pThreshNorm = 150;            % Threshold for normal external pressure (as dev from median)
p.tBufferNorm = 120;            % Make non-normal any sample within this of a non-normal (calibration) interval
p.nMadNorm = 6;                 % If pressure in normal (external) mode is this many MADs from the median, deem it non-normal

% Process calibration parameters
p.daMax = 0.5;                 % During a calibration successive samples will not change by more than this
p.tCalLim = [190 250];            % Time limits for calibration in seconds since start of stable output

%% Load Data

POBS_dir=dir('../../hikurangi/processed_POBS/');
POBS_list={POBS_dir.name}';
file_check=cellfun(@(v)v(1),POBS_list);
i_list=find(~eq(file_check,'.'));
for i=i_list'
    % if ~strcmp(POBS_dir(i).name,'POBS09') && ~strcmp(POBS_dir(i).name,'POBS7_co-located-w-QA15')
    %     continue
    % end
    
    % load POBS files
    data=load([POBS_dir(i).folder '/' POBS_dir(i).name '/' POBS_dir(i).name '.mat']);
    t_data=load([POBS_dir(i).folder '/' POBS_dir(i).name '/' POBS_dir(i).name '_triax.mat']);

    if any(isnan(data.a1))
        inan=find(isnan(data.a1));
        if inan(end)<86400 % first day already unusable b/c of deck time and descent
            data.Ta1(1:inan(end))=[]; data.Ta2(1:inan(end))=[];
            data.a1(1:inan(end))=[]; data.a2(1:inan(end))=[];
            data.ta(1:inan(end))=[];
        elseif inan(1)>(length(data.ta)-86400) % last day likewise unusable
            data.Ta1(inan(1):end)=[]; data.Ta2(inan(1):end)=[];
            data.a1(inan(1):end)=[]; data.a2(inan(1):end)=[];
            data.ta(inan(1):end)=[];
        elseif inan(1)<86400 && inan(end)>(length(data.ta)-86400) % cut first and last
            icut=[1:inan(1)-1,inan',inan(end)+1:length(data.ta)];
            data.Ta1(icut)=[]; data.Ta2(icut)=[];
            data.a1(icut)=[]; data.a2(icut)=[];
            data.ta(icut)=[];
        elseif length(inan)==1 % if only one point, just interpolate it
            data.Ta1(inan-1:inan+1)=nanmedian(data.Ta1(inan-5:inan+5));
            data.Ta2(inan-1:inan+1)=nanmedian(data.Ta2(inan-5:inan+5));
            data.a1(inan-1:inan+1)=nanmedian(data.a1(inan-5:inan+5));
            data.a2(inan-1:inan+1)=nanmedian(data.a2(inan-5:inan+5));
        elseif length(inan)<=600 && inan(end)-inan(1)==length(inan)-1 % continuous segments shorter than 10 minutes can be interpolated
            data.Ta1(inan(1)-1:inan(end)+1)=nanmedian(data.Ta1(inan(1)-5:inan(end)+5));
            data.Ta2(inan(1)-1:inan(end)+1)=nanmedian(data.Ta2(inan(1)-5:inan(end)+5));
            data.a1(inan(1)-1:inan(end)+1)=nanmedian(data.a1(inan(1)-5:inan(end)+5));
            data.a2(inan(1)-1:inan(end)+1)=nanmedian(data.a2(inan(1)-5:inan(end)+5));
        else % special case for POBS7 (updated for clock correction)
            data.Ta1(inan(end)-1:inan(end)+1)=nanmedian(data.Ta1(inan(end)-5:inan(end)+5));
            data.Ta2(inan(end)-1:inan(end)+1)=nanmedian(data.Ta2(inan(end)-5:inan(end)+5));
            data.a1(inan(end)-1:inan(end)+1)=nanmedian(data.a1(inan(end)-5:inan(end)+5));
            data.a2(inan(end)-1:inan(end)+1)=nanmedian(data.a2(inan(end)-5:inan(end)+5));
            inan=inan(1:end-2);
            data.Ta1(1:inan(end))=[]; data.Ta2(1:inan(end))=[];
            data.a1(1:inan(end))=[]; data.a2(1:inan(end))=[];
            data.ta(1:inan(end))=[];
        end
    end

    % find calibrations using 1 Hz data
    [ta,calInfo1,calInfo2,data.lNormState1,data.lNormState2] = find_cal(data.ta,data.a1,data.a2,p);
    
    % truncate POBS file appropriately
    [~,ia,~]=intersect(data.ta,ta);
    if length(ia)~=length(ta) % if these don't match something went wrong
        keyboard
    end
    data.ta=data.ta(ia);
    data.Ta1=data.Ta1(ia);
    data.Ta2=data.Ta2(ia);
    data.a1=data.a1(ia);
    data.a2=data.a2(ia);
    
    % analyze calibrations
    if isempty(calInfo1.t0)
        fprintf(['No calibrations found in ' POBS_dir(i).name '\n\n'])
        calInfoAll1=[];
        calInfoAll2=[];
    else
        % [calInfoAll1,calInfoAll2] = analyze_cals(data,calInfo1,calInfo2,p,POBS_dir(i).name); % with plots
        [calInfoAll1,calInfoAll2] = analyze_cals(data,calInfo1,calInfo2,p); % without plots
    end

    % add triax data
    data.tt=t_data.tt;
    data.xt=t_data.xt;
    data.yt=t_data.yt;
    data.zt=t_data.zt;
    data.Tt=t_data.Tt;

    % ensure pressure and triax span same interval
    [~,ia,ib]=intersect(round(data.ta,6),round(data.tt,6));
    data.ta=data.ta(ia);
    data.Ta1=data.Ta1(ia);
    data.Ta2=data.Ta2(ia);
    data.a1=data.a1(ia);
    data.a2=data.a2(ia);
    data.lNormState1=data.lNormState1(ia);
    data.lNormState2=data.lNormState2(ia);
    data.tt=data.tt(ib);
    data.xt=data.xt(ib);
    data.yt=data.yt(ib);
    data.zt=data.zt(ib);
    data.Tt=data.Tt(ib);

    % analyze barometer (same intervals as calibrations)
    % there are a few bad times thrown in -- just cut them out
    ibad = find(diff(data.tb)>10);
    data.tb(ibad)=[];
    data.b(ibad)=[];
    data.Tb(ibad)=[];

    if isempty(calInfo1.t0)
        barInfoAll=[];
    else
        % barInfoAll = analyze_bar(data,calInfoAll1,p,POBS_dir(i).name); % with plots
        barInfoAll = analyze_bar(data,calInfoAll1,p); % without plots
    end

    % decimate data to 1 sample/minute
    [~,dataDec.Ta1,~,~] = downsample_uneven(data.ta,data.Ta1,1/24/60);
    [~,dataDec.Ta2,~,~] = downsample_uneven(data.ta,data.Ta2,1/24/60);
    [~,dataDec.a1,~,~] = downsample_uneven(data.ta,data.a1,1/24/60);
    [dataDec.ta,dataDec.a2,~,~] = downsample_uneven(data.ta,data.a2,1/24/60);
    [~,dataDec.xt,~,~] = downsample_uneven(data.tt,data.xt,1/24/60);
    [~,dataDec.yt,~,~] = downsample_uneven(data.tt,data.yt,1/24/60);
    [~,dataDec.zt,~,~] = downsample_uneven(data.tt,data.zt,1/24/60);
    [dataDec.tt,dataDec.Tt,~,~] = downsample_uneven(data.tt,data.Tt,1/24/60);

    % % write 1 Hz data to text file
    % sname=['../pressure_data/datafiles_1Hz/' POBS_dir(i).name];
    % writematrix([round(data.ta,6),round(data.a1,4),round(data.Ta1,4),round(data.a2,4),round(data.Ta2,4)],sname)

    % write 1 Hz data to text file
    sdir=['/Volumes/ExtremeSSD/decimated/year1/' POBS_dir(i).name];
    if ~exist(sdir,'dir')
        mkdir(sdir)
    end
    sname=[sdir '/' POBS_dir(i).name];

    tstr=datestr(data.tb,30);
    T=table(tstr,data.b*100,data.Tb);
    T.Properties.VariableNames(1:3)={'t (yyyymmddTHHMMSS)','Pb (Pa)','T (C)'};
    writetable(T,[sname '_BAR'])

    tstr=datestr(data.ta,30);
    T=table(tstr,data.a1*100,data.Ta1,data.a2*100,data.Ta2);
    T.Properties.VariableNames(1:5)={'t (yyyymmddTHHMMSS)','P1 (Pa)','T1 (C)','P2 (Pa)','T2 (C)'};
    writetable(T,[sname '_APG'])

    T=table(tstr,data.xt,data.yt,data.zt,data.Tt);
    T.Properties.VariableNames(1:5)={'t (yyyymmddTHHMMSS)','ax (m/s^2)','ay (m/s^2)','az (m/s^2)','T (C)'};
    writetable(T,[sname '_TAX'])

    % save(['../pressure_data/' POBS_dir(i).name],'data','dataDec','calInfoAll1','calInfoAll2','barInfoAll','-v7.3')

    % plot calibration values and barometer readings
    if ~isempty(calInfo1.t0)
        figure(3); clf; hold on
        plot(calInfoAll1.t0,calInfoAll1.pCal,'o-','markersize',10,'linewidth',1)
        plot(calInfoAll2.t0,calInfoAll2.pCal+20,'s-','markersize',10,'linewidth',1)
        ylabel('P_A_P_G (hPa)')
        yyaxis right
        plot(barInfoAll.t0p,barInfoAll.pCal,'x-k','markersize',10,'linewidth',1)
        legend('APG1','APG2','Barometer')
        ylabel('P_b_a_r (hPa)')
        % ylim([949 969])
        datetick('x',6)
        set(gca,'fontsize',14)
        box on; grid on

        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];
        % print(['../figures/calibrations/' POBS_dir(i).name '_cal_summary'],'-dpng','-r300')
        % print(['../figures/calibrations/' POBS_dir(i).name '_cal_summary'],'-depsc','-vector')

        figure(3); clf; hold on
        barcor1=calInfoAll1.pCal-barInfoAll.pCal;
        barcor2=calInfoAll2.pCal-barInfoAll.pCal;
        plot(calInfoAll1.t0,barcor1-barcor2,'-ok','markersize',10,'linewidth',1)
        ylabel('P (hPa)')
        legend('Gauge 1–Gauge2','location','northeast')
        datetick('x',6)
        set(gca,'fontsize',14)
        box on; grid on

        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];
        % print(['../figures/calibrations/' POBS_dir(i).name '_barcor_difference'],'-dpng','-r300')
        % print(['../figures/calibrations/' POBS_dir(i).name '_barcor_difference'],'-depsc','-vector')

        if strcmp(POBS_dir(i).name,'POBS2')
            figure(3); clf;
            subplot(211); hold on
            plot(calInfoAll1.t0,calInfoAll1.pCal,'o','markersize',10,'linewidth',2)
            plot(calInfoAll2.t0,calInfoAll2.pCal,'s','markersize',10,'linewidth',2)
            plot(barInfoAll.t0p,zeros(size(barInfoAll.pCal))+865,'xk','markersize',7,'linewidth',1)
            ylabel('P (hPa)')
            yyaxis right
            plot(dataDec.ta(3600:end),dataDec.Ta1(3600:end),'color',[0 0.4470 0.7410],'linewidth',0.1)
            plot(dataDec.ta(3600:end),dataDec.Ta2(3600:end),'color','r','linewidth',0.1)
            ylabel('T (C)')
            legend('Gauge 1','Gauge2','Flat Barometer','location','northwest')
            datetick('x',6)
            title([POBS_dir(i).name ' calibrations'])
            set(gca,'fontsize',14)
            box on; grid on
            subplot(212); hold on
            barcor1=calInfoAll1.pCal-zeros(size(barInfoAll.pCal));
            barcor2=calInfoAll2.pCal-zeros(size(barInfoAll.pCal));
            plot(calInfoAll1.t0,barcor1-barcor1(end),'o','markersize',10,'linewidth',2)
            plot(calInfoAll2.t0,barcor2-barcor2(end),'s','markersize',10,'linewidth',2)
            ylabel('P (hPa)')
            legend('Gauge 1','Gauge2','location','northeast')
            datetick('x',6)
            title('Flat barometer calibrations')
            set(gca,'fontsize',14)
            box on; grid on

            fh=gcf;
            fh.PaperUnits='inches';
            fh.PaperPosition=[0 0 8.5 11];
            % print(['../figures/calibrations/' POBS_dir(i).name '_cal_summary_nobarometer'],'-dpng','-r300')

            figure(3); clf;
            subplot(211); hold on
            plot(calInfoAll1.t0,calInfoAll1.pCal,'o','markersize',10,'linewidth',2)
            plot(calInfoAll2.t0,calInfoAll2.pCal,'s','markersize',10,'linewidth',2)
            plot(barInfoAll.t0p,(barInfoAll.T-mean(barInfoAll.T))*3.71+865,'xk','markersize',7,'linewidth',1)
            ylabel('P (hPa)')
            yyaxis right
            plot(dataDec.ta(3600:end),dataDec.Ta1(3600:end),'color',[0 0.4470 0.7410],'linewidth',0.1)
            plot(dataDec.ta(3600:end),dataDec.Ta2(3600:end),'color','r','linewidth',0.1)
            ylabel('T (C)')
            legend('Gauge 1','Gauge2','T-scaled','location','northwest')
            datetick('x',6)
            title([POBS_dir(i).name ' calibrations'])
            set(gca,'fontsize',14)
            box on; grid on
            subplot(212); hold on
            barcor1=calInfoAll1.pCal-(barInfoAll.T-mean(barInfoAll.T))*3.71;
            barcor2=calInfoAll2.pCal-(barInfoAll.T-mean(barInfoAll.T))*3.71;
            plot(calInfoAll1.t0,barcor1-barcor1(end),'o','markersize',10,'linewidth',2)
            plot(calInfoAll2.t0,barcor2-barcor2(end),'s','markersize',10,'linewidth',2)
            ylabel('P (hPa)')
            legend('Gauge 1','Gauge2','location','northeast')
            datetick('x',6)
            title('Scaled temperature-corrected calibrations')
            set(gca,'fontsize',14)
            box on; grid on

            fh=gcf;
            fh.PaperUnits='inches';
            fh.PaperPosition=[0 0 8.5 11];
            % print(['../figures/calibrations/' POBS_dir(i).name '_cal_summary_temperature'],'-dpng','-r300')
        end

        % useful for plotting just one set of calibrations
        % figure(5); clf;
        % subplot(211); hold on
        % p1=calInfoAll2.pCal-mean(calInfoAll2.pCal);
        % plot(calInfoAll2.t0,p1+20,'rs','markersize',10,'linewidth',1)
        % p2=barInfoAll.pCal-mean(barInfoAll.pCal);
        % plot(barInfoAll.t0p,p2+15,'rx','markersize',7,'linewidth',1)
        % plot(calInfoAll2.t0,p1-p2,'ko','markersize',10,'linewidth',2)
        % ylabel('P (hPa)')
        % legend('APG','Barometer','Corrected APG','location','northeast')
        % datetick('x',6)
        % title([POBS_dir(i).name ' calibrations'])
        % set(gca,'fontsize',14)
        % box on; grid on
    end
end