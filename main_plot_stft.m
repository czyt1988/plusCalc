%绘制一个数据的3d频谱
close all;
%双容0.5m间距2
%'rawData';
%'saMainFreFilterStruct';
structName = 'incrementDenoisingStruct';
Fres = getFromDataStruct(dataStruct,'Fre',structName);
Mags = getFromDataStruct(dataStruct,'Mag',structName);
multFreMag = getFromDataStruct(dataStruct,'multFreMag',structName);
section = 1:15;
figure
%subplot(1,2,2)
plot_3DSpectrum(Fres,Mags,section);
hold on;
multFreMag10 = multFreMag(1,section);
plot3(ones(length(multFreMag10),1).*10,section,multFreMag10,'.-');
xlabel('frequecy(Hz)');
xlim([0,50]);
ylabel('measuring point number');
zlabel('amplitude');
grid on;
view(29,50);
set(gcf,'color','w');

figure
%subplot(1,2,1)
color_of_line = [44,75,247]./255;
color_of_marker_face = [108,57,232]./255;
plot(section,multFreMag10,'color',color_of_line,'Marker','<'...
    ,'MarkerEdgeColor',color_of_line,'MarkerFaceColor',color_of_marker_face)
grid on;
xlabel('measuring point number');
ylabel('amplitude');
d = axis;
hold on;
fill([6,7,7,6,6],[0,0,d(4),d(4),0],'r'...
    ,'FaceColor',[229,33,82]./255,'FaceAlpha',0.3 ...
    ,'EdgeColor',[250,72,85]./255)
xlim(d(1:2))
ylim(d(3:4))
set(gcf,'unit','pixels','position',[200,200,350,250]);
set(gcf,'color','w');



