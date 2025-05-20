function [calInfoOut1,calInfoOut2] = analyze_cals(data,calInfoIn1,calInfoIn2,p,makePlots)
% Analyze calibration intervals to get discrete calibration values
%
% Usage
%   [calInfoOut1,calInfoOut2] = analyze_cals(data,calInfoIn1,calInfoIn2,p,makePlots)
%
% Input
%     data       - Data structure with fields ta, a1, a2, Ta1, Ta2
%     calInfoIn1  - Output structure from 'find_cal.m' for Gauge 1 with following fields
%                       i0 - index of first sample of each calibration
%                       i1 - index of last sample of each calibration
%                       t0 - datetime of each calibration start
%                       aMed - median pressure of calibration interval 
%                       range80 - 80th percentile range of pressure during calibration
%                       complex - indicates that the calibration was complex (warrants manual investigation)
%      calInfoIn2  - Output structure for Gauge 2
%      p         - Parameter structure with fields
%                       daMax - Calibration starts when successive samples change by less than this
%                       tCalLim - Time limits (s) of calibrations samples to use relative to start time
%      makePlots - If instrument name as string is given, plots will be made and saved
%
% Outputs
%      calInfoOut1 - As for calInfoIn1 but with additional fields
%                       i0p - index of first sample of each calibration (reduced window)
%                       i1p - index of last sample of each calibration (reduced window)
%                       length - number of samples in calibration interval (reduced window)
%                       t0p - datetime of first sample (reduced window)
%                       pCal - median pressure of calibration interval
%                       pCalTCor - median pressure corrected for temperature
%                       T - median temperature during calibration interval
%                       Tstd - standard deviation of temperature
%                       duration - calibration interval length in s
%      calInfoOut2 - As for calInfoOut1 but for Gauge 2
%

if nargin<5
  makePlots = [];
else
    figure(111); clf; hold on
    figure(112); clf; hold on
end

%% GAUGE 1

calInfoOut1 = calInfoIn1;
calInfoOut1.i0p = NaN(size(calInfoOut1.i0));
calInfoOut1.i1p = NaN(size(calInfoOut1.i0));
calInfoOut1.t0p = NaN(size(calInfoOut1.i0));
calInfoOut1.pCal = NaN(size(calInfoOut1.i0));
calInfoOut1.pCalTCor = NaN(size(calInfoOut1.i0));
calInfoOut1.length = NaN(size(calInfoOut1.i0));
calInfoOut1.T = NaN(size(calInfoOut1.i0));
calInfoOut1.Tstd = NaN(size(calInfoOut1.i0));
calInfoOut1.duration = NaN(size(calInfoOut1.i0));

%Loop through calibrations
for i = 1:length(calInfoOut1.i0)

    % separate directories by experiment year
    if data.ta(1)<datenum(2023,09,01)
        pltdir='figures/calibrations';
    elseif data.ta(1)>datenum(2023,09,01)
        pltdir='figures_Y2/calibrations';
    end

    if ~calInfoOut1.complex(i)
        % Find start and end of calibration based on changes in pressure
        index = find(abs(diff(data.a1(calInfoOut1.i0(i):calInfoOut1.i1(i)))) < p.daMax);
        while (index(2) - index(1)) ~=1
            index = index(2:end);
        end
        while (index(end) - index(end-1)) ~=1
            index = index(1:end-1);
        end
        i0p = calInfoOut1.i0(i)+index(1);
        i1p = calInfoOut1.i0(i)+index(end);

        tp = (data.ta(i0p:i1p) - data.ta(i0p))*86400;
        Tp = data.Ta1(i0p:i1p);
        p0 = data.a1(i0p:i1p);
        p0_cor = p0 - (Tp-p.TRef1)*p.dadT1;

        % temperature scaling
        %%%% DOES NOT APPRECIABLY IMPROVE CALIBRATIONS %%%%
        % G=[Tp-mean(Tp),ones(size(tp))];
        % lamda_list=linspace(1,180,1000);
        % me=[]; e_fit=[]; stds=[];
        % for jj=1:1000
        %     gexp=exp(-tp/lamda_list(jj)); gexp(gexp<10^-7)=0;
        %     Ge=[G,gexp];
        %     me(:,jj)=inv(Ge'*Ge)*Ge'*p0;
        %     e_fit(:,jj)=Ge*me(:,jj);
        %     stds(jj)=std(p0-e_fit(:,jj));
        % end
        % [~,imin]=min(stds);
        % pmod=e_fit(:,imin);
        % e_cor=p0-pmod;
        % m=[me(:,imin);lamda_list(imin)];
        % p0_cor=p0-G(:,1)*m(1);
        % disp([num2str(m(1)) ' hPa/C' newline])

        index = tp>=p.tCalLim(1) &  tp<=p.tCalLim(2);
        pCal = median(p0(index));
        pCalTCor = median(p0_cor(index));

        calInfoOut1.i0p(i) = i0p;
        calInfoOut1.i1p(i) = i1p;
        calInfoOut1.length(i) = i1p-i0p+1;
        calInfoOut1.t0p(i) = data.ta(i0p);
        calInfoOut1.pCal(i) = pCal;
        calInfoOut1.pCalTCor(i) = pCalTCor;
        calInfoOut1.T(i) = median(Tp);
        calInfoOut1.Tstd(i) = std(Tp);
        calInfoOut1.duration(i) = tp(end)-tp(1);
    end

    % Plots of pressure and temperature during calibration
    if ~isempty(makePlots)
        figure(101); clf

        subplot(211); hold on
        yyaxis left
        plot(tp,p0,'b');
        ylabel('P (hPa)')
        lim_y=ylim;
        lim_x=xlim;
        hp=patch([190 250 250 190],[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.8 0.8 0.8]);
        hp.EdgeColor='none'; hp.FaceVertexAlphaData=0.2; hp.FaceAlpha='flat';
        ylim(lim_y)
        set(gca,'fontsize',14)
        yyaxis right
        plot(tp,Tp,'r');
        ylabel('T (C)')
        title({[makePlots ', Gauge 1, calibration #' num2str(i)]; datestr(floor(data.ta(calInfoOut1.i0(i))))});
        set(gca,'fontsize',14)
        box on; grid on

        subplot(212); hold on
        plot(tp(index),p0(index),'b','linewidth',1)
        xlim(lim_x)
        yline(pCal,'k--','linewidth',1)
        ylabel('P (hPa)')
        set(gca,'fontsize',14)
        yyaxis right
        plot(tp(index),Tp(index))
        xlim(lim_x)
        ylabel('T (C)')
        xlabel('time (s)')
        set(gca,'fontsize',14)
        box on; grid on

        % temperature corrected
        % subplot(349)
        % hold on
        % plot(tg(index),p0_cor(index),'b','linewidth',1)
        % lim_x=xlim;
        % plot([lim_x(1) lim_x(2)],[gCalTCor gCalTCor]','k--','linewidth',1)
        % xlabel('Time, s')
        % ylabel('g (m/s^2)')
        % title('Corrected')
        % box on

        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 8.5 11];
        if ~exist(['../' pltdir '/' makePlots '/gauge1'],'dir')
            mkdir(['../' pltdir '/' makePlots '/gauge1'])
        end
        print(['../' pltdir '/' makePlots '/gauge1/' makePlots '_G1_cal' num2str(i)],'-dtiff','-r100')

        figure(111)
        palt=p0-pCal;
        plot(palt,'color',[0.8 0.8 0.8])
        psum1{i}=palt;
        if i==length(calInfoOut1.i0)
            l=min(cellfun(@length,psum1));
            pcat=cell2mat(cellfun(@(v)v(1:l),psum1,'UniformOutput',false));
            psum=sum(pcat,2)/length(psum1);
            pstd=std((pcat-psum)')';

            h1=plot(psum,'k','linewidth',1);
            legend(h1,'mean transient')
            set(gca,'fontsize',14)
            ylabel('P (hPa)')
            xlabel('time (s)')
            title([makePlots ' G1'])
            box on; grid on
            
            fh=gcf;
            fh.PaperUnits='inches';
            fh.PaperPosition=[0 0 11 8.5];
            print(['../' pltdir '/' makePlots '/gauge1/' makePlots '_G1_stack'],'-dpng','-r100')
            pause(1)

            ylim([-0.2 0.2])
            print(['../' pltdir '/' makePlots '/gauge1/' makePlots '_G1_stack_zoom'],'-dpng','-r100')

            figure(11); clf
            subplot(211); hold on
            plot(psum,'k','linewidth',1)
            lim_y=ylim;
            hp=patch([190 250 250 190],[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.8 0.8 0.8]);
            hp.EdgeColor='none'; hp.FaceVertexAlphaData=0.2; hp.FaceAlpha='flat';
            set(gca,'fontsize',14)
            ylabel('P (hPa)')
            xlabel('Time (s)')
            title('Gauge 1')
            yyaxis right
            plot(pstd,'r','linewidth',1)
            ylabel('Standard Deviation (hPa)')
            box on; grid on
        end
    end
    if false % demonstrative plots for manuscript (worked well with POBS-12)
        figure(5); clf
        subplot(211); hold on
        plot((data.ta(i0p-136:i1p+100) - data.ta(i0p-136))*86400,data.a1(i0p-136:i1p+100),'linewidth',1)
        plot((data.ta(i0p-136:i1p+100) - data.ta(i0p-136))*86400,data.a2(i0p-136:i1p+100),'linewidth',1)
        ylim([137600 137950])
        xlim([0 505])
        set(gca,'fontsize',14)
        set(gca,'xticklabels',[])
        box on; grid on
        subplot(212); hold on
        plot((data.ta(i0p-136:i1p+100) - data.ta(i0p-136))*86400,data.a2(i0p-136:i1p+100),'linewidth',1)
        plot((data.ta(i0p-136:i1p+100) - data.ta(i0p-136))*86400,data.a2(i0p-136:i1p+100),'linewidth',1)
        ylim([650 1000])
        xlim([0 505])
        set(gca,'fontsize',14)
        ylabel('P (hPa)')
        xlabel('Time (s)')
        box on; grid on
        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];
        print(['../' pltdir '/demo/cal_a'],'-dpng','-r300')
        print(['../' pltdir '/demo/cal_a'],'-depsc','-vector')

        figure(6); clf; hold on
        fill([327 386 386 327],[942.66 942.66 942.84 942.84],[0.95 0.95 0.95],'edgecolor','none')
        plot((data.ta(i0p+75:i1p-1) - data.ta(i0p-136))*86400,data.a1(i0p+75:i1p-1),'linewidth',1,'color',[0 114 189]/255)
        ylim([942.66 942.84])
        set(gca,'fontsize',14)
        ylabel('P1 (hPa)')
        yyaxis right
        plot((data.ta(i0p+75:i1p-1) - data.ta(i0p-136))*86400,data.a2(i0p+75:i1p-1),'r','linewidth',1)
        ylabel('P2 (hPa)')
        xlabel('Time (s)')
        xlim([205 410])
        box on; grid on

        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];
        print(['../' pltdir '/demo/cal_b'],'-dpng','-r300')
        print(['../' pltdir '/demo/cal_b'],'-depsc','-vector')
    end
end

%% GAUGE 2

calInfoOut2 = calInfoIn2;
calInfoOut2.i0p = NaN(size(calInfoOut2.i0));
calInfoOut2.i1p = NaN(size(calInfoOut2.i0));
calInfoOut2.t0p = NaN(size(calInfoOut2.i0));
calInfoOut2.pCal = NaN(size(calInfoOut2.i0));
calInfoOut2.pCalTCor = NaN(size(calInfoOut2.i0));
calInfoOut2.length = NaN(size(calInfoOut2.i0));
calInfoOut2.T = NaN(size(calInfoOut2.i0));
calInfoOut2.Tstd = NaN(size(calInfoOut2.i0));
calInfoOut2.duration = NaN(size(calInfoOut2.i0));

%Loop through calibrations
for i = 1:length(calInfoOut2.i0)
    if ~calInfoOut2.complex(i)
        % Find start and end of calibration based on changes in pressure
        index = find(abs(diff(data.a2(calInfoOut2.i0(i):calInfoOut2.i1(i)))) < p.daMax);
        while (index(2) - index(1)) ~=1
            index = index(2:end);
        end
        while (index(end) - index(end-1)) ~=1
            index = index(1:end-1);
        end
        i0p = calInfoOut2.i0(i)+index(1);
        i1p = calInfoOut2.i0(i)+index(end);

        tp = (data.ta(i0p:i1p) - data.ta(i0p))*86400;
        Tp = data.Ta2(i0p:i1p);
        p0 = data.a2(i0p:i1p);
        p0_cor = p0 - (Tp-p.TRef2)*p.dadT2;

        index = tp>=p.tCalLim(1) &  tp<=p.tCalLim(2);
        pCal = median(p0(index));
        pCalTCor = median(p0_cor(index));

        calInfoOut2.i0p(i) = i0p;
        calInfoOut2.i1p(i) = i1p;
        calInfoOut2.length(i) = i1p-i0p+1;
        calInfoOut2.t0p(i) = data.ta(i0p);
        calInfoOut2.pCal(i) = pCal;
        calInfoOut2.pCalTCor(i) = pCalTCor;
        calInfoOut2.T(i) = median(Tp);
        calInfoOut2.Tstd(i) = std(Tp);
        calInfoOut2.duration(i) = tp(end)-tp(1);
    end

    % Plots of pressure and temperature during calibration
    if ~isempty(makePlots)
        figure(102); clf

        subplot(211); hold on
        yyaxis left
        plot(tp,p0,'b');
        ylabel('P (hPa)')
        lim_y=ylim;
        lim_x=xlim;
        hp=patch([190 250 250 190],[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.8 0.8 0.8]);
        hp.EdgeColor='none'; hp.FaceVertexAlphaData=0.2; hp.FaceAlpha='flat';
        ylim(lim_y)
        set(gca,'fontsize',14)
        yyaxis right
        plot(tp,Tp,'r');
        ylabel('T (C)')
        title({[makePlots ', Gauge 2, calibration #' num2str(i)]; datestr(floor(data.ta(calInfoOut2.i0(i))))});
        set(gca,'fontsize',14)
        box on; grid on

        subplot(212); hold on
        plot(tp(index),p0(index),'b','linewidth',1)
        xlim(lim_x)
        yline(pCal,'k--','linewidth',1)
        ylabel('P (hPa)')
        set(gca,'fontsize',14)
        yyaxis right
        plot(tp(index),Tp(index))
        xlim(lim_x)
        ylabel('T (C)')
        xlabel('time (s)')
        set(gca,'fontsize',14)
        box on; grid on

        % temperature corrected
        % subplot(349)
        % hold on
        % plot(tg(index),p0_cor(index),'b','linewidth',1)
        % lim_x=xlim;
        % plot([lim_x(1) lim_x(2)],[gCalTCor gCalTCor]','k--','linewidth',1)
        % xlabel('Time, s')
        % ylabel('g (m/s^2)')
        % title('Corrected')
        % box on

        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 8.5 11];
        if ~exist(['../' pltdir '/' makePlots '/gauge2'],'dir')
            mkdir(['../' pltdir '/' makePlots '/gauge2'])
        end
        print(['../' pltdir '/' makePlots '/gauge2/' makePlots '_G2_cal' num2str(i)],'-dtiff','-r100')

        figure(112)
        palt=p0-pCal;
        plot(palt,'color',[0.8 0.8 0.8])
        psum2{i}=palt;
        if i==length(calInfoOut2.i0)
            l=min(cellfun(@length,psum2));
            pcat=cell2mat(cellfun(@(v)v(1:l),psum2,'UniformOutput',false));
            psum=sum(pcat,2)/length(psum2);
            pstd=std((pcat-psum)')';

            h2=plot(psum,'k','linewidth',1);
            legend(h2,'mean transient')
            set(gca,'fontsize',14)
            ylabel('P (hPa)')
            xlabel('time (s)')
            title([makePlots ' G2'])
            box on; grid on
            
            fh=gcf;
            fh.PaperUnits='inches';
            fh.PaperPosition=[0 0 11 8.5];
            print(['../' pltdir '/' makePlots '/gauge2/' makePlots '_G2_stack'],'-dpng','-r100')
            pause(1)

            ylim([-0.2 0.2])
            print(['../' pltdir '/' makePlots '/gauge2/' makePlots '_G2_stack_zoom'],'-dpng','-r100')

            figure(11)
            subplot(212); hold on
            plot(psum,'k','linewidth',1)
            lim_y=ylim;
            hp=patch([190 250 250 190],[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.8 0.8 0.8]);
            hp.EdgeColor='none'; hp.FaceVertexAlphaData=0.2; hp.FaceAlpha='flat';
            set(gca,'fontsize',14)
            ylabel('P (hPa)')
            xlabel('Time (s)')
            title('Gauge 2')
            yyaxis right
            plot(pstd,'r','linewidth',1)
            ylabel('Standard Deviation (hPa)')
            box on; grid on

            fh=gcf;
            fh.PaperUnits='inches';
            fh.PaperPosition=[0 0 8.5 11];
            print(['../' pltdir '/' makePlots '/G1_G2_deviation'],'-dpng','-r100')
        end
    end
end