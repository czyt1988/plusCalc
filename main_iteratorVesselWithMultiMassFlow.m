function recorderCell = main_iteratorVesselWithMultiMassFlow
%迭代单一缓冲罐在不同转速下的情况
rcv = 0.11;% 相对余隙容积
pressureRadio=1.5;%压力比（排气压力/吸气压力）
k=1.4;%绝热指数
DCylinder=250/1000;%缸径m
dPipe=98/1000;%管内径m
crank=140/1000;%曲柄长
connectingRod=1.075;%连杆长度
%50deg,0.15MPa(A),Density=1.617
outDensity = 1.617;%1.9167;%环境25度绝热压缩到0.2的温度对应密度
fs = 600;
acousticVelocity = 345;%声速
L1 = 3;
L2 = 6;
l = 0.01;
Dpipe = 0.098;%管道直径（m）
Dv = 0.372;%缓冲罐的直径（m）
Lv = 1.1;%缓冲罐总长
rpm = 250:500;
useCalcTopFreIndex = [1:20];%参与计算的最高频率数，如[1:20]代表前20个最大的频率进行计算，nan为全部参与
accuracy = 1;%计算的精度默认为1，就是每隔1m取一个点
isDamping = 1;%是否有阻尼
isOpening = 0;%是否开口
coeffFriction = 0.04;%管道摩察系数
dcpss = getDefaultCalcPulsSetStruct();
dcpss.calcSection = [0.3,0.7];
dcpss.fs = fs;
dcpss.isHp = 0;
dcpss.f_pass = 7;%通过频率5Hz
dcpss.f_stop = 5;%截止频率3Hz
dcpss.rp = 0.1;%边带区衰减DB数设置
dcpss.rs = 30;%截止区衰减DB数设置
totalSecond = 5;
%%
totalPoints = fs*totalSecond+1;

index = 1;index_rpm = index;
recorderCell{1,index_rpm} = '转速';

index = index + 1;index_mass = index;
recorderCell{1,index_mass} = '质量流量';

index = index + 1;index_time = index;
recorderCell{1,index_time} = '时间';

index = index + 1;index_avg_mass = index;
recorderCell{1,index_avg_mass} = '等效质量流量';

index = index + 1;index_des = index;
recorderCell{1,index_des} = '描述';
index = index + 1;index_x = index;
recorderCell{1,index_x} = 'x值';
index = index + 1;index_plus = index;
recorderCell{1,index_plus} = '压力脉动';
index = index + 1;index_1fre = index;
recorderCell{1,index_1fre} = '1倍频';
index = index + 1;index_2fre = index;
recorderCell{1,index_2fre} = '2倍频';
index = index + 1;index_3fre = index;
recorderCell{1,index_3fre} = '3倍频';
index = index + 1;index_beforemaxPlus = index;
recorderCell{1,index_beforemaxPlus} = '罐前压力脉动最大值';
index = index + 1;index_aftermaxPlus = index;
recorderCell{1,index_aftermaxPlus} = '罐后压力脉动最大值';

recorderCell{1,index_des+8} = '描述';
recorderCell{1,index_x+8} = 'x值';
recorderCell{1,index_plus+8} = '压力脉动';
recorderCell{1,index_1fre+8} = '1倍频';
recorderCell{1,index_2fre+8} = '2倍频';
recorderCell{1,index_3fre+8} = '3倍频';
recorderCell{1,index_beforemaxPlus+8} = '罐前压力脉动最大值';
recorderCell{1,index_aftermaxPlus+8} = '罐后压力脉动最大值';
for i=1:length(rpm)
    recorderCell{i+1,index_rpm} = rpm(i);
    [massFlow,~,recorderCell{i+1,index_avg_mass}] = ...
        massFlowMaker(DCylinder,dPipe,rpm(i)...
        ,crank,connectingRod,outDensity,'rcv',rcv,'k',k,'pr',pressureRadio,'fs',fs);
    while length(massFlow)<totalPoints %生成指定长度的数据
        massFlow = [massFlow,massFlow];
    end
    massFlow = massFlow(1:totalPoints);
    time = linspace(0,totalSecond,totalPoints);%0:1/fs:totalSecond;
    recorderCell{i+1,index_mass} = massFlow;
    recorderCell{i+1,index_time} = time;
    
    calcDatas = funIteratorVesselPipeLinePlusCalc(time,massFlow...
        ,fs,rpm(i)/60*2,acousticVelocity,L1,L2,l,Dpipe,Dv,Lv,'Dv'...
        ,'multfre',[rpm(i)/60*2,rpm(i)/60*2*2,rpm(i)/60*2*3]...%计算的倍频
        ,'isOpening',isOpening...
        ,'useCalcTopFreIndex',useCalcTopFreIndex...
        ,'isDamping',isDamping...
        ,'coeffFriction',coeffFriction...
        ,'accuracy',accuracy...
        ,'dcpss',dcpss...
    );
    recorderCell{i+1,index_des} = calcDatas{3,1};
    recorderCell{i+1,index_x} = calcDatas{3,2};
    recorderCell{i+1,index_plus} = calcDatas{3,3};
    recorderCell{i+1,index_1fre} = calcDatas{3,4};
    recorderCell{i+1,index_2fre} = calcDatas{3,5};
    recorderCell{i+1,index_3fre} = calcDatas{3,6};
    recorderCell{i+1,index_beforemaxPlus} = calcDatas{3,7};
    recorderCell{i+1,index_aftermaxPlus} = calcDatas{3,8};
    
    recorderCell{i+1,index_des+8} = calcDatas{2,1};
    recorderCell{i+1,index_x+8} = calcDatas{2,2};
    recorderCell{i+1,index_plus+8} = calcDatas{2,3};
    recorderCell{i+1,index_1fre+8} = calcDatas{2,4};
    recorderCell{i+1,index_2fre+8} = calcDatas{2,5};
    recorderCell{i+1,index_3fre+8} = calcDatas{2,6};
    recorderCell{i+1,index_beforemaxPlus+8} = calcDatas{2,7};
    recorderCell{i+1,index_aftermaxPlus+8} = calcDatas{2,8};
    
end

currentPath = fileparts(mfilename('fullpath'));
savePath = fullfile(currentPath,'iteratorVesselWithMultiMassFlowRes.mat');
save(savePath,'recorderCell');
%system('shutdown -s');
end