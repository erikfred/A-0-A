function [tds,xds,nUsed,nNaN] = downsample_uneven(t,x,dt,timeBasis)
% Downsamples an uneven time series by averaging within bins
% This second version does not specify output sample times but instead obtains 
% them from the input data.  It also returns the number of samples in each average
% 
% Usage
%   [tds,xds,nUsed,nNaN] = downsample_uneven(t,x,dt,timeBasis)
%
% Inputs
%   t         - Times of samples in time series
%   x         - Time series values (matrix okay)
%   dt        - Length of each averaging window (output sample interval) in days
%               Averaging intervals will always extend from N*dt to (N+1)*dt
%   timeBasis - String to indicate how times of samples are referenced to averaging interval
%                 'start'
%                 'middle' (default)
%                 'end'
%
% Outputs
%   tds   - Times of samples in output time series (time of center of window)
%   xds   - Output averaged time series (NaN's where no input samples in time interval)
%   nUsed - Number of samples averaged for each output in xds
%   nNaN  - Number of Nans in each interval
%
% This algorithm is fast and accounts for NaNs
%
% Modified May, 2020 to limit rounding errors from using cumulative sum differences to get means


% Small value compared with sample interval (in days)
small = 1e-4/86400;
% Cumulative sums of very large time series can lead to rounding when differenced to get means
% Limit non zero cumsums to this many samples
nCumSumMax = 1e6;

% Default arguments
if nargin<4
  timeBasis = 'middle';
end

% Check inputs have consistent dimensions
xDim = size(x);
if ~any(length(t(:))==xDim)
  error('downsample_uneven: Input dimensions of time and data inconsistent')
end

% Put time series into column orientation
t = t(:);
if xDim(1) == length(t)
  flip = false;
else
  flip = true;
  x = x';
  xDim = xDim([2 1]);
end

% Create time series that encloses samples intervals
t0 = floor((t(1)+small)/dt)*dt;
t1 = ceil((t(end)-small)/dt)*dt;
tds = t0:dt:t1+dt/1e6;
tds = tds(:);

% Create output time series as NaNs 
xds = NaN(length(tds)-1,xDim(2));
nUsed = zeros(length(tds)-1,xDim(2));
nNaN = zeros(length(tds)-1,xDim(2));

% For each input time determine the appropriate index in the output time series
% First and last sample may have an invalid index
if strcmpi(timeBasis,'end') 
  index = ceil((t-small-t0)/dt);
  while index(1)==0
    index = index(2:end);
    t = t(2:end);
    x = x(2:end,:);
    xDim(1) = xDim(1)-1;
  end
else
  index = floor((t+small-t0)/dt)+1;
  while index(end)==length(tds)
    index = index(1:end-1);
    t = t(1:end-1);
    x = x(1:end-1,:);
    xDim(1) = xDim(1)-1;
  end
end

% Find the indices corresponding to a change in output sample including the first and last input sample
iEnd = [true; diff(index)>0; true];
iEnd = find(iEnd);

% Determine the sample that each averaged value represents 
index = index(iEnd(1:end-1));

% Number of observations in each average
n = iEnd(2:end)-iEnd(1:end-1);

% Deal with NaNs
if any(isnan(x(:)))
  % Logical matix of NaN
  nanLogical = isnan(x);
  % Cumulative sum of NaNs with extra 0 at start
  nanLogicalCum = [zeros(1,xDim(2)); cumsum(nanLogical)];
  % Number of NaNs in each avarage
  nn = (nanLogicalCum(iEnd(2:end),:) - nanLogicalCum(iEnd(1:end-1),:));
  % Assign average values to output
  for j = 1:xDim(2)
    nNaN(index,j) = nn(:,j);
  end
  
% Zero out NaNs  
  x(nanLogical) = 0;
  
else
  nn = zeros(length(n),xDim(2));
end

% Set values with data to zero
for j = 1:xDim(2)
  xds(index,j) = 0;
end

% Loop to minimize rounding errors
nZero = 0;
loop = true;
while loop
  
  x(1:nZero,:) = 0;

  % Cumulative sum with extra 0 at start
  xc = [zeros(1,xDim(2)); cumsum(x)];

  % Averaged values (accounting for NaNs)
  xm = (xc(iEnd(2:end),:) - xc(iEnd(1:end-1),:)) ./ (repmat(n,1,xDim(2)) - nn);

  % Get rid of ±infinite values when there are samples with no non-NaN observations
  xm(isinf(xm)) = NaN;

  % Assign average values to output
  for j = 1:xDim(2)
    xds(index,j) = xds(index,j).*(iEnd(1:end-1)<=nZero) + xm(:,j).*(iEnd(1:end-1)>nZero);
  end

  nZero = nZero + nCumSumMax;
  if nZero>=xDim(1)
    loop = false;
  end
  
end

nUsed(index,:) = (repmat(n,1,xDim(2)) - nn);
nNaN(index,:) = nn;

% Sample time is average of interval
if strcmpi(timeBasis,'start')
  tds = tds(1:end-1);
elseif strcmpi(timeBasis,'end')
  tds = tds(2:end);
else
  tds = (tds(1:end-1) + tds(2:end))/2;
end

% Put outputs into same orientation as inputs
if flip
  xds = xds';
  tds = tds';
end