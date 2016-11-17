figure;
hold on;
meaPoints = 11;
for i = 2:size(recorderCell,1)
    x = recorderCell{i,1};
    x = x * ones(1,meaPoints);
    y = 1:meaPoints;
    z = recorderCell{i,7};
   plot3(x,y,z);
end
view(117,44);
set(gcf,'color','w');
title('µ¥Ò»»º³å¹ÞÑ¹Á¦Âö¶¯');

figure
z = ones(meaPoints,size(recorderCell,1)-1);
for i = 2:size(recorderCell,1)
    z(1:meaPoints,i-1) = recorderCell{i,7};
end
x = cell2mat(recorderCell(2:end,1));
y = 1:meaPoints;
[X,Y] = meshgrid(x,y);
contourf(X,Y,z);
colorbar;
set(gcf,'color','w');
title('µ¥Ò»»º³å¹ÞÑ¹Á¦Âö¶¯');

figure
z = ones(meaPoints,size(recorderCell,1)-1);
for i = 2:size(recorderCell,1)
    z(1:meaPoints,i-1) = recorderCell{i,8};
end
x = cell2mat(recorderCell(2:end,1));
y = 1:meaPoints;
[X,Y] = meshgrid(x,y);
contourf(X,Y,z);
colorbar;
set(gcf,'color','w');
title('µ¥Ò»»º³å¹Þ1±¶Æµ');

figure
z = ones(meaPoints,size(recorderCell,1)-1);
for i = 2:size(recorderCell,1)
    z(1:meaPoints,i-1) = recorderCell{i,9};
end
x = cell2mat(recorderCell(2:end,1));
y = 1:meaPoints;
[X,Y] = meshgrid(x,y);
contourf(X,Y,z);
colorbar;
set(gcf,'color','w');
title('µ¥Ò»»º³å¹Þ2±¶Æµ');

figure
z = ones(meaPoints,size(recorderCell,1)-1);
for i = 2:size(recorderCell,1)
    z(1:meaPoints,i-1) = recorderCell{i,10};
end
x = cell2mat(recorderCell(2:end,1));
y = 1:meaPoints;
[X,Y] = meshgrid(x,y);
contourf(X,Y,z);
colorbar;
set(gcf,'color','w');
title('µ¥Ò»»º³å¹Þ3±¶Æµ');

figure
hold on;
x = zeros(1,size(recorderCell,1)-1);
y = zeros(1,size(recorderCell,1)-1);
for i = 2:size(recorderCell,1)
    x(:,i-1) = recorderCell{i,1};   
    y(:,i-1) = recorderCell{i,11}; 
end
plot(x,y); 
title('µ¥Ò»»º³å¹Þ¹ÞÇ°Ñ¹Á¦Âö¶¯');

figure
hold on;
x = zeros(1,size(recorderCell,1)-1);
y = zeros(1,size(recorderCell,1)-1);
for i = 2:size(recorderCell,1)
    x(:,i-1) = recorderCell{i,1};   
    y(:,i-1) = recorderCell{i,12}; 
end
plot(x,y); 
title('µ¥Ò»»º³å¹Þ¹ÞºóÑ¹Á¦Âö¶¯');
