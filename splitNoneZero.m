function splitData = splitNoneZero( data )
%分解非0数据到一个除了它自己全都是0的数据里
%[1,0,2,0,1] =>
%
%[1,0,0,0,0]
%[0,0,2,0,0]
%[0,0,0,0,1]

index = find(data~=0);
n = length(data);
for i=1:length(index)
    splitData(:,i) = zeros(n,1);
    splitData(index(i),i) = data(index(i));
end

end

