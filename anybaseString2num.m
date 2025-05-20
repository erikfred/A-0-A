function outnum = anybaseString2num(base,string,start)
% converts any number to base string from Ascii 33(!) to 126(~)
% start is 33 '!' by default and can be changed to any other value > 33
% For binary conversion use base 2 and start 48, for octal use base 8 and
% start 48

if ~exist('start','var')
    start = 33;
elseif(start<33 || start>126)
    warning('anyBase2string:startOutOfRange','Characters below 33 and above 126 are not printable, changing start to 33');
    start = 33;
end
if(base>(126-start+1))
    warning('anyBase2string:maxBaseExceeded','Maximum allowed base is (127-start (33 default)) = %d\n Using %d',(126-start+1),(126-start+1));
    base = (126-start+1);
end
outnum = zeros(size(string));
for kk = 1:numel(string)
number = 0;
thisCharArray = char(string{kk});
for i=numel(thisCharArray):-1:1
    number = number + (thisCharArray(i)-start)*base^(i-1);
end
outnum(kk) = number;
end
end
