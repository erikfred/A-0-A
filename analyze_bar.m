function [barInfoOut] = analyze_bar(data,calInfoIn,p,makePlots)
% Analyze calibration intervals to get discrete calibration values
%
% Usage
%   [barInfoOut] = analyze_bar(data,calInfoIn,p,makePlots)
%
% Input
%     data       - Data structure with fields ta, a1, a2, Ta1, Ta2
%     calInfoIn  - Output structure from 'find_cal.m' for either gauge with following fields
%                       i0 - Vector with first sample of each calibration
%                       i1 - Vector with last sample of each calibration
%                       t0 - Vector with serial data number of each calibration start
%                       aMed - Vector with median pressure of calibration interval 
%                       range80 - Vector with 80th percentile range of pressure during calibration
%                       complex - Indicates that the calibration was complex (warrants manual investigation)
%      p         - Parameter structure with fields
%                       daMax - Calibration starts when successive samples change by less than this
%                       tCalLim - Time limits (s) of calibrations samples to use relative to start time
%      makePlots - If instrument name as string is given, plots will be made and saved
%
% Outputs
%      barInfoOut - As for calInfoIn but with additional fields
%                       i0p - index of first sample of each calibration (reduced window)
%                       i1p - index of last sample of each calibration (reduced window)
%                       length - number of samples in calibration interval (reduced window)
%                       t0p - datetime of first sample (reduced window)
%                       pCal - median pressure of calibration interval
%                       pCalTCor - median pressure corrected for temperature
%                       T - median temperature during calibration interval
%                       Tstd - standard deviation of temperature
%                       duration - calibration interval length in s
%

if nargin<4
  makePlots = [];
end

barInfoOut.i0p = NaN(size(calInfoIn.i0));
barInfoOut.i1p = NaN(size(calInfoIn.i0));
barInfoOut.t0p = calInfoIn.t0p;
barInfoOut.pCal = NaN(size(calInfoIn.i0));
barInfoOut.pCalTCor = NaN(size(calInfoIn.i0));
barInfoOut.length = NaN(size(calInfoIn.i0));
barInfoOut.T = NaN(size(calInfoIn.i0));
barInfoOut.Tstd = NaN(size(calInfoIn.i0));
barInfoOut.duration = NaN(size(calInfoIn.i0));

%Loop through calibrations
for i = 1:length(barInfoOut.i0p)
    if ~calInfoIn.complex(i)
        % Find start and end of calibration interval from APG indices
        [~,i0p] = min(abs(data.tb-calInfoIn.t0p(i)));
        [~,i1p] = min(abs(data.tb-(calInfoIn.t0p(i)+calInfoIn.duration(i)/60/60/24)));

        tp = (data.tb(i0p:i1p) - data.tb(i0p))*86400;
        Tp = data.Tb(i0p:i1p);
        p0 = data.b(i0p:i1p);
        p0_cor = p0 - (Tp-p.TRefb)*p.dadTb;

        index = tp>=p.tCalLim(1) &  tp<=p.tCalLim(2);
        pCal = median(p0(index));
        pCalTCor = median(p0_cor(index));

        barInfoOut.i0p(i) = i0p;
        barInfoOut.i1p(i) = i1p;
        barInfoOut.length(i) = i1p-i0p+1;
        barInfoOut.t0p(i) = data.tb(i0p);
        barInfoOut.duration(i) = tp(end);
        barInfoOut.pCal(i) = pCal;
        barInfoOut.pCalTCor(i) = pCalTCor;
        barInfoOut.T(i) = median(Tp);
        barInfoOut.Tstd(i) = std(Tp);
    end

    % Plots of pressure during calibration
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
        title({[makePlots ', Barometer, calibration #' num2str(i)]; datestr(floor(data.tb(barInfoOut.i0p(i))))});
        set(gca,'fontsize',14)
        box on; grid on

        subplot(212); hold on
        plot(tp(index),p0(index),'b','linewidth',1)
        xlim(lim_x)
        yline(pCal,'k--','linewidth',1)
        ylabel('P (hPa)')
        set(gca,'fontsize',14)
        xlim(lim_x)
        xlabel('time (s)')
        set(gca,'fontsize',14)
        box on; grid on

        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 8.5 11];
        if ~exist(['../figures/calibrations/' makePlots '/barometer'],'dir')
            mkdir(['../figures/calibrations/' makePlots '/barometer'])
        end
        % print(['../figures/calibrations/' makePlots '/barometer/' makePlots '_cal' num2str(i)],'-dtiff','-r100')
    end
end