function main_plotMultFre()
    clc;
    clear;
    close all;
    currentPath = fileparts(mfilename('fullpath'));
    %% 加载数据
    dvThrData = getAnalysisDataFun('dvthrfre');%双容理论
    dvSimData = getAnalysisDataFun('dvsimfre');%双容模拟
    dvExpData = getAnalysisDataFun('dvexpfre');%双容实验
    strSimData = getAnalysisDataFun('strfre');%直管模拟数据
    ovSameAsDoubleSim = getAnalysisDataFun('ovsimfre');%单容模拟数据   
    [dv05Sim.fre,dv05Sim.mag]=getSimulationData(fullfile(currentPath,'originData\26米双容间距0.5米多测点.mat'));
    [dv09Sim.fre,dv09Sim.mag]=getSimulationData(fullfile(currentPath,'originData\26米双容间距0.9米多测点.mat'));
    [dvUpSim.fre,dvUpSim.mag]=getSimulationData(fullfile(currentPath,'originData\26米双容间距0米向上叠加.mat'));
    [dvUpCrossSim.fre,dvUpCrossSim.mag]=getSimulationData(fullfile(currentPath,'originData\26米双容间距0米向上叠加交错.mat'));
    [dv09D2Sim.fre,dv09D2Sim.mag]=getSimulationData(fullfile(currentPath,'originData\26米双容间距1米D2=1.66.mat'));
    [ovSameAsOneSim.fre,ovSameAsOneSim.mag]=getSimulationData(fullfile(currentPath,'originData\26米单容.mat'));
    %[newSim.fre,newSim.mag]=getSimulationData(fullfile(currentPath,'originData\new.mat'));

    %% 分析阶数
    baseFre = 5:5:40;
    allowDis = 0.5;
    sepIndex = [14,15];
    startDrawPoint = 1;
    dvThrOrderMag = getOrderMag(dvThrData.fre,dvThrData.mag,baseFre,allowDis);
    dvThrOrderMag = dvThrOrderMag(:,startDrawPoint:end);
    dvThrRedRate = calcReduceRate(sepIndex+1-startDrawPoint,dvThrOrderMag');% 脉动抑制率
    
    
    dv0SimOrderMag = getOrderMag(dvSimData.fre,dvSimData.mag,baseFre,allowDis);
    dv0SimOrderMag = dv0SimOrderMag(:,startDrawPoint:end);
    dvSimRedRate = calcReduceRate(sepIndex+1-startDrawPoint,dv0SimOrderMag');
    
    dvExpOrderMag = getOrderMag(dvExpData.fre,dvExpData.mag,baseFre,allowDis);
    dvExpOrderMag = dvExpOrderMag(:,startDrawPoint:end);
    dvExpRedRate = calcReduceRate(sepIndex+1-startDrawPoint,dvExpOrderMag');
    
    dv05SimOrderMag = getOrderMag(dv05Sim.fre,dv05Sim.mag,baseFre,allowDis);
    dv05SimOrderMag = dv05SimOrderMag(:,startDrawPoint:end);
    dv05SimRedRate = calcReduceRate(sepIndex+1-startDrawPoint,dv05SimOrderMag');
    
    dv09SimOrderMag = getOrderMag(dv09Sim.fre,dv09Sim.mag,baseFre,allowDis);
    dv09SimOrderMag = dv09SimOrderMag(:,startDrawPoint:end);
    dv09SimRedRate = calcReduceRate(sepIndex+1-startDrawPoint,dv09SimOrderMag');
    
    dvUpSimOrderMag = getOrderMag(dvUpSim.fre,dvUpSim.mag,baseFre,allowDis);
    dvUpSimOrderMag = dvUpSimOrderMag(:,startDrawPoint:end);
    dvUpSimRedRate = calcReduceRate(sepIndex+1-startDrawPoint,dvUpSimOrderMag');
    
    dvUpCrossSimOrderMag = getOrderMag(dvUpCrossSim.fre,dvUpCrossSim.mag,baseFre,allowDis);
    dvUpCrossSimOrderMag = dvUpCrossSimOrderMag(:,startDrawPoint:end);
    dvUpCrossSimRedRate = calcReduceRate(sepIndex+1-startDrawPoint,dvUpCrossSimOrderMag');  
    
    strSimOrderMag = getOrderMag(strSimData.fre,strSimData.mag,baseFre,allowDis);
    strSimOrderMag(:,15:16) = [];
    strSimOrderMag(:,end) = 0;
    strSimOrderMag = strSimOrderMag(:,startDrawPoint:end);
    strSimRedRate = calcReduceRate(sepIndex+1-startDrawPoint,strSimOrderMag');
    
    ov2SimOrderMag = getOrderMag(ovSameAsDoubleSim.fre,ovSameAsDoubleSim.mag,baseFre,allowDis);
    ov2SimOrderMag = ov2SimOrderMag(:,startDrawPoint:end);
    ov2SimRedRate = calcReduceRate(sepIndex+1-startDrawPoint,ov2SimOrderMag');
    
    ov1SimOrderMag = getOrderMag(ovSameAsOneSim.fre,ovSameAsOneSim.mag,baseFre,allowDis);
    ov1SimOrderMag = ov1SimOrderMag(:,startDrawPoint:end);
    ov1SimRedRate = calcReduceRate(sepIndex+1-startDrawPoint,ov1SimOrderMag');
    
    dv09D2SimOrderMag = getOrderMag(dv09D2Sim.fre,dv09D2Sim.mag,baseFre,allowDis);
    dv09D2SimOrderMag = dv09D2SimOrderMag(:,startDrawPoint:end);
    dv09D2SimRedRate = calcReduceRate(sepIndex+1-startDrawPoint,dv09D2SimOrderMag');
    %newSimOrderMag = getOrderMag(newSim.fre,newSim.mag,baseFre,allowDis);
    %% 
    

    %% 绘制双容和直管的模拟的前n阶对比
%     index = 24;%index5=node1 ; index9=node4;index14 = node6;index15 = node7;index20 = node9;index24 = node11
%     d = joinVector(dv0SimOrderMag(:,index),ov1SimOrderMag(:,index),strSimOrderMag(:,index));
%     lengText = {'two-tank system','single-tank system','pipe without tanks'};
%     plotBarCmp(d,lengText);
%     ylim([0,4]);
%     xlim([0,9]);
%     set(gcf,'position',[300,300,330,230]);
   %% 绘制双容和等体积双容的模拟的前n阶对比
   if 1
        index = 24;%index5=node1 ; index9=node4;index14 = node6;index15 = node7;index20 = node9;index24 = node11
        d = joinVector(dv0SimOrderMag(:,index),ov2SimOrderMag(:,index),strSimOrderMag(:,index));
        lengText = {'two-tank unit','single-tank unit','pipe without tanks'};
        plotBarCmp(d,lengText);
        ylim([0,4.2]);
        xlim([0,9]);
        set(gcf,'position',[300,300,330,230]);
   end
    %% 绘制三维柱状图
    if 0
        figure
        hold on;
        indexs = [5,9,14,15,20,24];
        node = [1,4,6,7,9,11];
        count = 1;
        h = [];
        bar3Handle = {};
        for index = indexs
            d = joinVector(dv0SimOrderMag(:,index),ov2SimOrderMag(:,index),strSimOrderMag(:,index));
            h = bar3(d,'grouped');
            set(h,'UserData',node(count));
    %         edgeColor = [8,62,88;...
    %             126,85,3;...
    %             52,44,103]./255;
    %         colorFace = [14,158,230;...
    %             248,140,0;...
    %             109,91,240]./255;
            edgeColor = [0,70,57;...
                64,93,52;...
                126,126,52]./255;
            colorFace = [14,158,230;...
                129,192,102;...
                255,255,102]./255;
            set(h(1),'EdgeColor',edgeColor(1,:),'FaceColor',colorFace(1,:));
            set(h(2),'EdgeColor',edgeColor(2,:),'FaceColor',colorFace(2,:));
            set(h(3),'EdgeColor',edgeColor(3,:),'FaceColor',colorFace(3,:));
            bar3Handle{count} = h;
            count = count + 1;
        end
        hbar3=findobj(gca,'Tag','bar3'); % 所有的 bar3
        set(hbar3,{'XData'}, ...
            cellfun(@(X,dx)X+dx-1, ... % 更新x坐标: X+dx
            get(hbar3,{'XData'}), ...% 原x坐标: X
            get(hbar3,{'UserData'}),'Uni',0)) % x坐标增量: dx
        %dasp=daspect; daspect(dasp([1,1,3]))
        set(gca,'XTick',node);
        set(gca,'YTick',1:8);
        %set(gca,'XTickLabel',{'1','4','6','7','9','11'});
        %set(gca,'XTick',1:25);
        axis tight;
        grid on;% grid minor;
        %colormap hsv; caxis([1,inf])
        set(gcf,'color','w');
        set(gcf,'position',[300,300,430,230]);
        zlim([0,4]);
        ylim([0,9]);
        xlim ([0,12]);
        xlabel('node');
        ylabel('order');
        view(-142,48);
        legend()
        export_fig(gcf,fullfile(currentPath,'3dbar',sprintf('i%d',index)),'-m5');
    end
    %% 绘制主要测点的3种情况倍频连线
%     figure
%     indexs = [5,9,14,15,20,24];
%     meaPoint = [1,4,6,7,9,11];
%     count = 1;
%     for index = indexs
%         subplot(2,3,count);
%         h(1) = plot(1:8,dv0SimOrderMag(:,index));
%         hold on;
%         h(2) = plot(1:8,ov2SimOrderMag(:,index));
%         h(3) = plot(1:8,strSimOrderMag(:,index));
%         xlabel(sprintf('measuring point %d',meaPoint(count)));
%         count = count + 1;
%     end
%     set(gcf,'color','w');
%     set(gcf,'position',[300,300,330,230]);
%     zlim([0,4]);
%     ylim([0,9]);

        %export_fig(gcf,fullfile(currentPath,'3dbar',sprintf('i%d',index)),'-m5');
%     d = joinVector(dv0SimOrderMag(:,index),ov2SimOrderMag(:,index),strSimOrderMag(:,index));
%     bar3(d,'grouped');
    %% 绘制某一倍频的各个测点的对比图
%     order = 1;
%     d = joinVector(dv09SimOrderMag(order,:),dv09D2SimOrderMag(order,:),strSimOrderMag(order,:));
%     lengText = {'double vessel 1m distance result','double vessel 1m distance f=10Hz result','straight result'};
%     figure
%     plotLineCmp(d,lengText);  

    %% 0m，0.5m，0.9m，双容拼接，直管对比
%     lengText = {'DV with L2=0m','DV with L2=0.5m','DV with L2=1m'...
%         ,'pipe without DV'};
% %     lengText = {'double vessel spacing 0m result','double vessel spacing 0.5m result','double vessel spacing 1m result','straight result'};
% %       lengText = {'double vessels with L2=1m','double vessels wirh L2=1m V2=1.66m3','pipe without vessels'};
% %     lengText = {'double vessel 0m distance result','straight result'};
%     row = 3;
%     col = 3;
%     
%     for ii=1:8%绘制1到8阶
%          dcell{ii} = joinVector(dv0SimOrderMag(ii,:),dv05SimOrderMag(ii,:),dv09SimOrderMag(ii,:),strSimOrderMag(ii,:));
% %          dcell{ii} = joinVector(dv09SimOrderMag(ii,:),dv09D2SimOrderMag(ii,:),strSimOrderMag(ii,:));
% %          dcell{ii} = joinVector(dvSimOrderMag(ii,5:end),strSimOrderMag(ii,5:end));
%           X{ii} = startDrawPoint:(length(dcell{ii})+startDrawPoint-1);
%     end
%     h = plotLineCmpMult2(X,dcell,row,col,@plot);
%     subplotfun(@(xx) title(sprintf('order %d',xx)),row,col,8);
%     color_style = [245,18,103;...
%     18,175,134;...
%     255,99,56;...
%     131,54,229;...
%     36,100,196;...
%     237,144,10;...
%     
%     ]./255;
%     line_style = {'-';'--';':';'-.'};
%     for ii = 1:length(h)
%         LH = h{ii};
%         for jj = 1:length(LH)
%             set(LH(jj),'color',color_style(jj,:));
%             set(LH(jj),'LineStyle',line_style{jj});
%         end
%     end
%     legend(h{8},lengText,'Location','NorthEastOutside');
%     h = subplot(row,col,2);
%     set(h,'ytick',[0,2.5,5]);
%     subplot(row,col,4)
%     ylabel('amplitude(kPa)');
%     subplot(row,col,8)
%     xlabel('distance(m)');
    %% 绘制3阶，4阶的相对直管脉动抑制率
%     str3 = strSimOrderMag(3,:);
%     sim0str3 = (str3 - dv0SimOrderMag(3,:)) ./ str3 * 100;
%     sim05str3 = (str3 - dv05SimOrderMag(3,:)) ./ str3 * 100;
%     sim09str3 = (str3 - dv09SimOrderMag(3,:)) ./ str3 * 100;
%     figure
%     plotLineCmpWithX(X{3},[sim0str3',sim05str3',sim09str3']);
%     
%     baseCmp  = dv0SimOrderMag(3,:);
%     sim05str3 = (dv05SimOrderMag(3,:) - baseCmp)./baseCmp*100;
%     sim09str3 = (dv09SimOrderMag(3,:) - baseCmp)./baseCmp*100
%     figure
%     plotLineCmpWithX(X{3},[sim05str3',sim09str3']);
    %% 脉动抑制率
%     lengText = {'double vessel 0m distance result','double vessel 0.5m distance result','double vessel 1m distance result'...
%         ,'straight result'};
%     data = joinVector(dvSimRedRate,dv05SimRedRate,dv09SimRedRate,strSimRedRate);
%     figure
%     plotLineCmpWithX(baseFre,data,lengText);
%     xlabel('frequency(Hz)');
%     ylabel('');
%     
    
end

function [mags,index] = findRangMax(Fre,Mag,baseFre,allowDis)
    for ii = 1:length(baseFre)
        k = Fre>(baseFre(ii)-allowDis) & Fre<(baseFre(ii)+allowDis);
        if isempty(find(k == 1, 1))
            mags(ii) = nan;
            index(ii) = nan;
        else
            [mags(ii),index(ii)] = max(Mag(k));%找到范围内最大的值
        end
    end
end

function [orderMag] = getOrderMag(Fres,Mags,baseFre,allowDis)
    for ii = 1:size(Fres,2)
        orderMag(:,ii) = findRangMax(Fres(:,ii),Mags(:,ii),baseFre,allowDis)';
    end
end

function ra = calcReduceRate(indexSperate,data)
    before = max(data(1:indexSperate(1),:));
    after = max(data(indexSperate(2):end,:));
    ra = after./before;
end

function h = plotBarCmp(d,lt)
% 绘制某测点的前10阶幅值
    h = bar(d,'hist');
    colormap summer
    legend(h,lt);
    xlabel('order');
    ylabel('amplitude (kPa)');
    set(gcf,'color','w');
end

function h = plotLineCmp(d,lt)
% 绘制某测点的前10阶幅值
    if nargin == 1
        isLegend = 0;
    else
        isLegend = 1;
    end
    hold on;
    r = colormap(lines);
    for ii = 1:size(d,2)
        h(ii) = plot(d(:,ii));
        set(h(ii),'color',r(ii,:));
    end
    if isLegend
        legend(h,lt);
    end
    grid on;

    xlabel('node');
    ylabel('amplitude (kPa)');
    set(gcf,'color','w');
end

function h = plotLineCmpWithX(x,d,lt)
% 对比图，带x
    if nargin == 2
        isLegend = 0;
    else
        isLegend = 1;
    end
    hold on;
    r = colormap(lines);
    for ii = 1:size(d,2)
        h(ii) = plot(x,d(:,ii));
        set(h(ii),'color',r(ii,:));
    end
    if isLegend
        legend(h,lt);
    end
    grid on;
    set(gcf,'color','w');
end

function h = plotLineCmpMult(dcell,r,c,funPlotHandle)
% 绘制某测点的前10阶幅值
    count = r*c;
    count = min(count,length(dcell));
    hold on;
    for ii = 1:count
        subplot(r,c,ii)
        d = dcell{ii};
        h{ii} = funPlotHandle(d);
    end
end

function h = plotLineCmpMult2(x,y,r,c,funPlotHandle)
% 绘制某测点的前10阶幅值
    count = r*c;
    count = min(count,length(y));
    clr = colormap(lines);
    for ii = 1:count
        subplot(r,c,ii)
        d = y{ii};
        xh = x{ii};
        hold on;
        for jj = 1:size(d,2)
            temp(jj) = funPlotHandle(xh,d(:,jj));
            set(temp(jj),'color',clr(jj,:));
            set(temp(jj),'LineWidth',1.5);
        end
        grid on;
        h{ii} = temp;
    end
    set(gcf,'color','w');
end
