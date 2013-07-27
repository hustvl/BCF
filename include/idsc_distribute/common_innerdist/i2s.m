function s = i2s(n,len)

if ~exist('len') | n>=10^len
	s = int2str(n);
else
	s = [];
	if n==0			lenn=1;
	else			lenn = floor(log10(n))+1;
	end
	for ii=lenn+1:len
		s = [s '0'];
	end
	s = [s int2str(n)];
end
	