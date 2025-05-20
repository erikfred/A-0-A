function [smooth] = Z_godin(data)
% 3/20/2012  Parker MacCready, based on jfilt by Jonathan Lilly 2005
%
% This applies the 24-24-25 Godin Filter to a vector, or a matrix in which
% each COLUMN is a data record.  It differs from the version in
% Z_dasfilt.m in that it makes a single filter shape of length 71, instead
% of using repeated passes of a boxcar. This results in more NaN-padding.
% The other difference is it accepts matrices of column vectors.
% 
% It returns a dataset of the same size you started with, padded whth
% NaN's at the ends.
%
% ** use ONLY with hourly data! **

test = size(data);
% transpose input if necessary
if test(1)<test(2)
    data=data';
end
    
% This is the shape given in Emery and Thomson (1997) Eqn. (5.10.37)
k = [0:11];
filter = NaN * ones(71,1);
filter(36:47) = (0.5/(24*24*25))*(1200-(12-k).*(13-k)-(12+k).*(13+k));
k = [12:35];
filter(48:71) = (0.5/(24*24*25))*(36-k).*(37-k);
filter(1:35) = flipud(filter(37:71));

% % alternatively you can make the shape this way,
% it is identical but less efficient to create
% aa = NaN * ones(25,24,24);
% for ii = [-12:12]
%     aa(ii+13,:,:) = toeplitz([ii:-1:ii-23],[ii:23+ii]);
% end
% X = [-35:35]; % there are 71 values, from -35 to 35, representing indices
% % relative to 0, the offest from the filtered point
% filter = hist(aa(:),X)'; % with distribution N
% filter=filter./sum(filter);

n=length(filter);
smooth=zeros(size(data));
a=round(n/2);
N=size(data,1);
for i=1:size(data,2)
    temp=conv(data(:,i),filter);
    smooth(:,i)=temp(a:a+N-1);
    smooth(1:n,i)=nan*ones(n,1);
    smooth(N-n+1:N,i)=nan*ones(n,1);
end

% return to original shape
if test(1)<test(2)
    smooth=smooth';
end


