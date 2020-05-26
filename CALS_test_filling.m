%%%%%%%%%%   ����  2019.07.18       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cost-sensitive active learning for incomplete data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%======================================================================
clc;
clear all;
warning off;
tic;
%%======================================================================

%% ��ȡԭʼ����
dataOri = 'E:\�о���\����ѧϰ\����\ʵ�鲿��\���²���\excel\irisDecision.xlsx';
% dataOri = 'E:\CALS_Thesis\Well data\Well05.xlsx';
[dataOri,~] = xlsread(dataOri);
[n, d] = size(dataOri);  %���ݼ�n�У�dά

%% ����ȱʧֵԤ�����
fillingMethod = 1;     %% 1-��Ȩƽ��ֵԤ���2-����������ֵԤ�(�������ֲ岹)��3-BPCA���
A = 0.1 ; B = 1 - A;   %% ����ȱʧ�ʣ�

switch fillingMethod
    case 1
        %% �������ȱʧָʾ����
        unlabeledData = dataOri(:,1:d-1);  %%�ޱ�ǩ����
        unlabeledData = unlabeledData';
        [un,uc] = size(unlabeledData);
        randNumber = [0 1];               %%�������һ��������0��1������
%         A = 0.5 ; B = 1 - A;
        prob = [A B];                      %%����0��1���ֵĸ��ʣ�
        randomX = randsrc(un,uc,[randNumber; prob]);
        
        %% �������ȱʧ����
        data = unlabeledData.*randomX;
        data = data';
        data = [data,dataOri(:,end)];
        randomX = [ones(1,uc); randomX];        %ÿ��������Ҫ����һ��x0=1;
    case 2
        %% �������ȱʧָʾ����
        unlabeledData = dataOri(:,1:d-1);  %%�ޱ�ǩ����
        unlabeledData = unlabeledData';
        [un,uc] = size(unlabeledData);
        randNumber = [NaN 1];               %%�������һ��������NaN��1������
%         A = 0.5 ; B = 1 - A;
        prob = [A B];                      %%����0��1���ֵĸ��ʣ�
        randomX = randsrc(un,uc,[randNumber; prob]);
        
        %% �������ȱʧ����
        data = unlabeledData.*randomX;
        data = data';
        data = [data,dataOri(:,end)];
        randomX = [ones(1,uc); randomX];        %ÿ��������Ҫ����һ��x0=1;
        
        data_F = fillmissing(data,'spline');
        % spline���ֶ�����������ֵ��pchip�����ηֶ�����������ֵ��makima������Akima����Hermite��ֵ��        
        data_filling = [ones(1,n);(data_F(:,1:end-1))'];
        data(find(isnan(data) == 1)) = 0;
        randomX(find(isnan(randomX) == 1)) = 0;
    case  3
         %% �������ȱʧָʾ����
        unlabeledData = dataOri(:,1:d-1);  %%�ޱ�ǩ����
        unlabeledData = unlabeledData';
        [un,uc] = size(unlabeledData);
        randNumber = [NaN 1];               %%�������һ��������NaN��1������
%         A = 0.5 ; B = 1 - A;
        prob = [A B];                      %%����0��1���ֵĸ��ʣ�
        randomX = randsrc(un,uc,[randNumber; prob]);
        
        %% �������ȱʧ����
        data = unlabeledData.*randomX;
        data = data';
        data = [data,dataOri(:,end)];
        randomX = [ones(1,uc); randomX];        %ÿ��������Ҫ����һ��x0=1;
        dataM = data(:,1:end-1);
        
        data_F = BPCA_filling(dataOri,dataM);
        % spling���ֶ�����������ֵ��pchip�����ηֶ�����������ֵ��makima������Akima����Hermite��ֵ��        
        data_filling = [ones(1,n);(data_F(:,1:end-1))'];
        data(find(isnan(data) == 1)) = 0;
        randomX(find(isnan(randomX) == 1)) = 0;
end

%%%%%%%%����dc����
%%%%%%%%%%%��һ��%%%%%%%%%%%%%%%%%%%%%%
% ratio = floor(0.1*n);      % 10%��ѵ������������ȡ��
% ratio = floor(sqrt(n));     % ����������ƽ������Ϊѵ�����ĸ���
ratio = 10;
distanceSelection = 1;     
normalization = 1;     % 0-�޹�һ����1-�й�һ����
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
[Value, Index] = min(result(:,5));     % �ҳ���С��ƽ������
resultMin = result(Index,:)
buyMin = buy{Index};     %����ı�ǩ

% xlswrite('E:\CALS_Thesis\Well data\incompleteData\Well05_50.xlsx',data);
% xlswrite('E:\CALS_Thesis\Well data\incompleteData\randomX\Well05_50.xlsx',randomX);

toc
