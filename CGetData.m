classdef CGetData
    %获取数据的类
    %   Detailed explanation goes here
    
    properties
        data = [];
        currentPath = [];
    end
    
    methods
        function td = CGetData()
            currentPath = fileparts(mfilename('fullpath'));
        end % 构造函数
        
        function setData(obj,dataType)
            switch(lower(dataType))
                case 'thrfre'
                    load(fullfile(currentPath,'理论-0m间距-f小于80.mat'));
                    [obj.data.fre,obj.data.mag] = fft_byColumn(pppp,3600);
                    obj.data.mag = obj.data.mag./1000;
            end
        end
    end
    
end

