function [t_out,calInfo1,calInfo2,lNormState1,lNormState2] = find_cal(t,a1,a2,p)
% Finds calibrations in an A-0-A pressure time series
%
% Usage
%   [calInfo1,calInfo2,lNormalState1,lNormState2] = find_cal(t,a1,a2,p)
%
% Inputs
%   t  - time series
%   a1 - pressure 1
%   a2 - pressure 2
%   p  - parameter structure with fields
%           pThresh          - System is in calibration when pressure is below this value
%           minTime4Cal      - A calibration requires at least this much time (s) in internal sensing
%           complexRange80   - If 80% of the calibration interval time series is not within this range
%                               from median then the calibration is complex and warrants manual inspection
%           pThreshNorm      - Max deviation from median for determining that the external
%                               pressure is as expected
%           tBufferNorm      - Adjacent samples within this buffer of a non-normal (calibration)
%                               interval are likewise deemed to be non-normal
%           nMadNorm         - If pressure in normal (external) state is this many MADs from the
%                               median, deem it non-normal
%
%
% Outputs
%  t_out         -
%  calInfo1      - Structure with information about Gauge 1 calibrations with following fields
%                   normalState = The normal state of the sensor (usually 1)
%                   i0 - Vector with first sample of each calibration
%                   i1 - Vector with last sample of each calibration
%                   state - Vector with state of each calibration (1 for external, 0 for internal)
%                   t0 - Vector with serial data number of each calibration start
%                   aMed - Vector with median pressure of calibration interval
%                   range80 - Vector with 80th percentile range of pressure during the calibration
%                   complex - Indicates that the calibration was complex (warrants manual inspection)
%                   roughDuration - Rough duration of calibration
%  calInfo2      - Structure with information about Gauge 2 calibrations with same fields
%  lNormalState1 - Logical vector of length t set to true when Gauge 1 is in its normal state
%  lNormalState2 - Logical vector of length t set to true when Gauge 2 is in its normal state
%

%% GAUGE 1

% state for each sample is set to
%    0 if measuring internal (calibration)
%    1 if measuring external (non-calibration)
%    2 if other (transition or something going wrong)
state = 2*ones(size(a1));
state(abs(a1-median(a1))<=p.pThreshNorm) = 1;
state(a1<=p.pThresh) = 0;

% beginning of time series will likely include deck & descent
% likewise, end may include ascent and deck
% trim such cases
istart=find(state==1,1);
iend=find(state==1,1,'last');
t=t(istart:iend);
a1=a1(istart:iend);
a2=a2(istart:iend);
state=state(istart:iend);

% Start and End samples of calibration intervals
i0 = find(diff(state)==-2)+1;
i1 = find(diff(state)==2);
% ensure these look right
if length(i0)~=length(i1)
    keyboard
elseif any(i1-i0<0)
    keyboard
end

% No calibrations found
if isempty(i0) || isempty(i1)
    calInfo1.i0 = [];
    calInfo1.i1 = [];
    calInfo1.t0 =  [];
    calInfo1.aMed =  [];
    calInfo1.range80 =  [];
    calInfo1.complex =  [];
    calInfo1.duration = [];
    lNormState1 = ones(size(t));

    calInfo2.i0 = [];
    calInfo2.i1 = [];
    calInfo2.t0 =  [];
    calInfo2.aMed =  [];
    calInfo2.range80 =  [];
    calInfo2.complex =  [];
    calInfo2.duration = [];
    lNormState2 = ones(size(t));

    t_out=t;
    return
end

% ensure each calibration interval is as long as expected
ikeep = (t(i1)-t(i0))>=p.minTime4Cal/86400;
if any(~ikeep)
    bads = datestr(t(i0(~ikeep)));
    for i=1:size(bads,1)
        warning(['Calibration on ' bads(i,:) ' shorter than expected'])
    end
    i0 = i0(ikeep);
    i1 = i1(ikeep);
end

% Set the calibration information structure
calInfo1.i0 = i0;
calInfo1.i1 = i1;
calInfo1.t0 = t(i0);

calInfo1.aMed = zeros(size(calInfo1.i0));
calInfo1.range80 = zeros(size(calInfo1.i0));
for i=1:length(calInfo1.i0)
    calInfo1.aMed(i) = median(a1(i0(i):i1(i)));
    temp = sort(a1(i0(i):i1(i)));
    calInfo1.range80(i) = temp(round(length(temp)*0.9)) - temp(round(length(temp)*0.1));
end

calInfo1.complex = calInfo1.range80>p.complexRange80;
calInfo1.roughDuration = (t(i1) - t(i0)) * 86400;

% Create a logical that is true for samples that are in the normal orienation,
% set apart from calibrations and without extreme values
medianA = median(a1);
madA = mad(a1);
% Normal state when observed pressure not to far from median value
lNormState1 = a1>medianA-p.nMadNorm*madA & a1<medianA+p.nMadNorm*madA;
% First samples when state is not normal
on  = find(diff(lNormState1)==-1)+1;
% Last sample when state is not normal
off = find(diff(lNormState1)==1);
% ensure these look right
if length(i0)~=length(i1)
    keyboard
elseif any(i1-i0<0)
    keyboard
end
% Make anything close to a non normal state, non-normal
for i=1:length(on)
    lNormState1(t>=t(on(i))-p.tBufferNorm/86400 & t<=t(off(i))+p.tBufferNorm/86400) = false;
end

%% GAUGE 2

% state for each sample is set to
%    0 if measuring internal (calibration)
%    1 if measuring external (non-calibration)
%    2 if other (transition or something going wrong)
state = 2*ones(size(a2));
state(abs(a2-median(a2))<=p.pThreshNorm) = 1;
state(a2<=p.pThresh) = 0;

% in theory, Gauge 2 should have same start/end effects as Gauge 1
% % beginning of time series will likely include deck & descent
% % likewise, end may include ascent and deck
% % trim such cases
% istart=find(state==1,1);
% iend=find(state==1,1,'last');
% t=t(istart:iend);
% a2=a2(istart:iend);
% a2=a2(istart:iend);
% state=state(istart:iend);

% Start and End samples of calibration intervals
i0 = find(diff(state)==-2)+1;
i1 = find(diff(state)==2);
% ensure these look right
if length(i0)~=length(i1)
    keyboard
elseif any(i1-i0<0)
    keyboard
end

% No calibrations found
if isempty(i0) || isempty(i1)
    calInfo2.i0 = [];
    calInfo2.i1 = [];
    calInfo2.t0 =  [];
    calInfo2.aMed =  [];
    calInfo2.range80 =  [];
    calInfo2.complex =  [];
    calInfo2.duration = [];
    lNormState2 = ones(size(t));
    return
end

% ensure each calibration interval is as long as expected
ikeep = (t(i1)-t(i0))>=p.minTime4Cal/86400;
if any(~ikeep)
    bads = datestr(t(i0(~ikeep)));
    for i=1:size(bads,1)
        warning(['Calibration on ' bads(i,:) ' shorter than expected'])
    end
    i0 = i0(ikeep);
    i1 = i1(ikeep);
end

% Set the calibration information structure
calInfo2.i0 = i0;
calInfo2.i1 = i1;
calInfo2.t0 = t(i0);

calInfo2.aMed = zeros(size(calInfo2.i0));
calInfo2.range80 = zeros(size(calInfo2.i0));
for i=1:length(calInfo2.i0)
    calInfo2.aMed(i) = median(a2(i0(i):i1(i)));
    temp = sort(a2(i0(i):i1(i)));
    calInfo2.range80(i) = temp(round(length(temp)*0.9)) - temp(round(length(temp)*0.1));
end

calInfo2.complex = calInfo2.range80>p.complexRange80;
calInfo2.roughDuration = (t(i1) - t(i0)) * 86400;

% Create a logical that is true for samples that are in the normal orienation,
% set apart from calibrations and without extreme values
medianA = median(a2);
madA = mad(a2);
% Normal state when observed pressure not to far from median value
lNormState2 = a2>medianA-p.nMadNorm*madA & a2<medianA+p.nMadNorm*madA;
% First samples when state is not normal
on  = find(diff(lNormState2)==-1)+1;
% Last sample when state is not normal
off = find(diff(lNormState2)==1);
% ensure these look right
if length(i0)~=length(i1)
    keyboard
elseif any(i1-i0<0)
    keyboard
end
% Make anything close to a non normal state, non-normal
for i=1:length(on)
    lNormState2(t>=t(on(i))-p.tBufferNorm/86400 & t<=t(off(i))+p.tBufferNorm/86400) = false;
end

% export trimmed time series
t_out = t;