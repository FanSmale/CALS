%%%%%%%%%%   汪敏  2019.07.18       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cost-sensitive active learning for incomplete data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%======================================================================
clc;
clear all;
warning off;
tic;
%%======================================================================

%% 读取原始数据
dataOri = 'E:\研究生\主动学习\论文\实验部分\重新测试\excel\irisDecision.xlsx';
% dataOri = 'E:\CALS_Thesis\Well data\Well05.xlsx';
[dataOri,~] = xlsread(dataOri);
[n, d] = size(dataOri);  %数据集n行，d维

%% 设置缺失值预填补方案
fillingMethod = 1;     %% 1-加权平均值预填补；2-其他样条插值预填补(包含多种插补)；3-BPCA填补；
A = 0.1 ; B = 1 - A;   %% 设置缺失率；

switch fillingMethod
    case 1
        %% 产生随机缺失指示数据
        unlabeledData = dataOri(:,1:d-1);  %%无标签数据
        unlabeledData = unlabeledData';
        [un,uc] = size(unlabeledData);
        randNumber = [0 1];               %%随机产生一定比例的0、1个数；
%         A = 0.5 ; B = 1 - A;
        prob = [A B];                      %%设置0、1出现的概率；
        randomX = randsrc(un,uc,[randNumber; prob]);
        
        %% 产生随机缺失数据
        data = unlabeledData.*randomX;
        data = data';
        data = [data,dataOri(:,end)];
        randomX = [ones(1,uc); randomX];        %每个样本都要增加一个x0=1;
    case 2
        %% 产生随机缺失指示数据
        unlabeledData = dataOri(:,1:d-1);  %%无标签数据
        unlabeledData = unlabeledData';
        [un,uc] = size(unlabeledData);
        randNumber = [NaN 1];               %%随机产生一定比例的NaN、1个数；
%         A = 0.5 ; B = 1 - A;
        prob = [A B];                      %%设置0、1出现的概率；
        randomX = randsrc(un,uc,[randNumber; prob]);
        
        %% 产生随机缺失数据
        data = unlabeledData.*randomX;
        data = data';
        data = [data,dataOri(:,end)];
        randomX = [ones(1,uc); randomX];        %每个样本都要增加一个x0=1;
        
        data_F = fillmissing(data,'spline');
        % spline：分段三次样条插值；pchip：保形分段三次样条插值；makima：修正Akima三次Hermite插值；        
        data_filling = [ones(1,n);(data_F(:,1:end-1))'];
        data(find(isnan(data) == 1)) = 0;
        randomX(find(isnan(randomX) == 1)) = 0;
    case  3
         %% 产生随机缺失指示数据
        unlabeledData = dataOri(:,1:d-1);  %%无标签数据
        unlabeledData = unlabeledData';
        [un,uc] = size(unlabeledData);
        randNumber = [NaN 1];               %%随机产生一定比例的NaN、1个数；
%         A = 0.5 ; B = 1 - A;
        prob = [A B];                      %%设置0、1出现的概率；
        randomX = randsrc(un,uc,[randNumber; prob]);
        
        %% 产生随机缺失数据
        data = unlabeledData.*randomX;
        data = data';
        data = [data,dataOri(:,end)];
        randomX = [ones(1,uc); randomX];        %每个样本都要增加一个x0=1;
        dataM = data(:,1:end-1);
        
        data_F = BPCA_filling(dataOri,dataM);
        % spling：分段三次样条插值；pchip：保形分段三次样条插值；makima：修正Akima三次Hermite插值；        
        data_filling = [ones(1,n);(data_F(:,1:end-1))'];
        data(find(isnan(data) == 1)) = 0;
        randomX(find(isnan(randomX) == 1)) = 0;
end

%%%%%%%%调整dc测试
%%%%%%%%%%%归一化%%%%%%%%%%%%%%%%%%%%%%
% ratio = floor(0.1*n);      % 10%的训练集，并向下取整
% ratio = floor(sqrt(n));     % 以样本数的平方根作为训练集的个数
ratio = 10;
distanceSelection = 1;     
normalization = 1;     % 0-无归一化；1-有归一化；
tempNumber = 1;

switch fillingMethod
    case 1
        for dc = 0.01:0.1:1
            [missRate,accuracy,totalCost,averCost,buyNumbers,buyLable]= fun_CALS(data,dataOri,randomX, dc,ratio,distanceSelection,normalization);
            result(tempNumber,:) = [dc,missRate,accuracy,totalCost,averCost,buyNumbers];
            buy(tempNumber) = {buyLable};
            toc
            tempNumber = tempNumber + 1
        end
    case 2
        for dc = 0.01:0.1:1
            [missRate,accuracy,totalCost,averCost,buyNumbers,buyLable]= fun_CALS_filling(data,dataOri,randomX, dc,ratio,distanceSelection,normalization,data_filling);
            result(tempNumber,:) = [dc,missRate,accuracy,totalCost,averCost,buyNumbers];
            buy(tempNumber) = {buyLable};
            toc
            tempNumber = tempNumber + 1
        end
    case 3
        for dc = 0.01:0.1:1
            [missRate,accuracy,totalCost,averCost,buyNumbers,buyLable]= fun_CALS_filling(data,dataOri,randomX, dc,ratio,distanceSelection,normalization,data_filling);
            result(tempNumber,:) = [dc,missRate,accuracy,totalCost,averCost,buyNumbers];
            buy(tempNumber) = {buyLable};
            toc
            tempNumber = tempNumber + 1
        end
end
result
[Value, Index] = min(result(:,5));     % 找出最小的平均代价
resultMin = result(Index,:)
buyMin = buy{Index};     %购买的标签

% xlswrite('E:\CALS_Thesis\Well data\incompleteData\Well05_50.xlsx',data);
% xlswrite('E:\CALS_Thesis\Well data\incompleteData\randomX\Well05_50.xlsx',randomX);

toc
