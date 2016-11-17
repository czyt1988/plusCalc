function issame = cmpfloat(d1,d2,mind)
	if nargin < 3
		mind = 1e-6;
	end
	issame = (d1<(d2+mind) && d1>(d2-mind));
end