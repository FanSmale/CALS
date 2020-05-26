function [missRate,accuracy,totalCost,averCost,queryNumbers,buyLable]=fun_CALS(data,dataOri,randomX,dc,ratio,distanceSelection,normalization)

%% ��ȡȱʧ����
[n, d] = size(data);  %���ݼ�n�У�dά
X = data(:,1:d-1);  %%ȱʧ�ޱ�ǩ����
Y = data(:,d);      %%��ʵ��ǩ

%% ��ȡԭʼ����
OriX = dataOri(:,1:d-1);  %%�����ޱ�ǩ����
oneArray = ones(1,n);
OriX = [oneArray; OriX']; %ÿ��������Ҫ����һ��x0=1;

%% ��ʼ������
H = cell(1,n);      %%��ʼ��һ����cell����
labels = [];
queryInstances = [];

%% ģ�Ͳ���
k = numel(unique(Y));           % ͳ�Ʊ�ǩ������
misCostValue = 2;  % ������������
teachCostValue = 1;     % ���ù����ǩ�Ľ�ʦ����
attributeCostValue = 0.2;     % ���ù������ԵĴ���
lambda = 1e-4;         % Weight decay parameter

upgradeTheta = 0;

%%======================================================================
%% ��ʼ��
for i = 1:1:n
    isLabeled(i) = 0;           %%���Ѵ����ʵ�����б��
    instanceCost(i) = 0;        %%�����ǩ�Ĵ���
    Labels(i) = 0;              %%��¼ʵ���ı�ǩ��Ϣ
    minallCost(i) = 0;          %%��ʼ����¼��С�ܴ��۵ľ���
    attributeCost(i) = 200;     %%��ʼ�����Դ���Ϊ�Ӵ�ֵ
end

%% �������
switch distanceSelection
    case 1
        Dists = pdist(X,'cityblock');  %% �����پ���
    case 2
        Dists = pdist(X,'cosine');  %% �н����Ҿ���
    case 3
        Dists = pdist(X,'euclidean');  %%ŷ����¾���
    case 4
        Dists = pdist(X,'minkowski');  %% �ɿƷ�˹������
    case 5
        Dists = pdist(X,'hamming');  %% ��������
    case 6
        Dists = pdist(X,'seuclidean');  %%��׼ŷ����¾���
    case 7
        Dists = pdist(X,'mahalanobis');  %%�����ŵ��˹����Mahalanobis distance ���Ͼ���
    case 8
        Dists = pdist(X,'chebychev');  %% �б�ѩ�����
    otherwise
        Dists = pdist(X,'jaccard');  %% �ܿ��¾���
end
Dists = squareform(Dists);
maxDist = max(max(Dists));

%% ��������ԣ�Ҳ���Ǽ����ܶȺ;���
dc = dc * maxDist;
rho = zeros(n, 1);
for i=1:n-1
    for j=i+1:n
        rho(i)=rho(i)+exp(-(Dists(i,j)/dc)*(Dists(i,j)/dc));
        rho(j)=rho(j)+exp(-(Dists(i,j)/dc)*(Dists(i,j)/dc));
    end
end

delta = zeros(n, 1);
master = -ones(n, 1);
[~, ordrho] = sort(rho, 'descend');
delta(ordrho(1)) = maxDist;
for i = 2:n
    delta(ordrho(i)) = maxDist;
    for j = 1:i-1
        if Dists(ordrho(i), ordrho(j)) < delta(ordrho(i))
            delta(ordrho(i)) = Dists(ordrho(i), ordrho(j));
            master(ordrho(i)) = ordrho(j);
        end
    end
end
gamma = rho .* delta;
[~, desInd] = sort(gamma, 'descend');               %%%%%%%desIndΪ���������е�����
maxGamma = max(gamma);

%% �������ݼ�,��һ������
X = [oneArray; X']; %ÿ��������Ҫ����һ��x0=1;
if normalization == 1
    X = mapminmax(X);
    X = X.*randomX;
    OriX = mapminmax(OriX);
end

%% �������Լ�Ȩƽ��ֵ
for i = 2:1:d
    nonAttributeNumbers(i) = length(nonzeros(X(i,:)));  %%%%%%%%%ÿ�����Է�ȱֵԪ�ظ���
end

for i = 2:1:d
    averValue(i-1) = sum(X(i,:))./nonAttributeNumbers(i);
end
averValue = averValue';
averValue = [1;averValue];

%%======================================================================
%%%%%%%%%%% ����1��������ʼѵ������ѵ��theta
%%%%%%%%%%% ������ʼѵ���� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%ͨ��ѡ������ķ�ʽ������ʼѵ���� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trainingSet = desInd(1:1:ratio);
buyLable = trainingSet;
trainingInput = X(:,trainingSet);     %%%%%%%ѵ����
labels = Y(trainingSet);              %%%%%%%ѵ������ǩ
Labels(trainingSet) = Y(trainingSet);
Labels = Labels';
isLabeled(trainingSet) = 1;

%%%%%%%%%%% ѵ��theta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = 0.0005 * randn(d, k);     % ���������ʼtheta
[cost, grad] = softmax_regression_vec(theta,trainingInput,labels,lambda);
options.maxIter = 100;
softmaxModel = softmaxTrain(d, k, lambda, trainingInput,labels, options);
theta = softmaxModel.optTheta;  %���»��theta����

%%======================================================================
%%%%%%%%%%% ����2����ѭ�� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while sum(isLabeled) ~= n
    
    %%%%%% ����2.1 ѵ�������£�����ѵ��ģ�� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if upgradeTheta == 1
        inputData = X(:,trainingSet);
        [cost, grad] = softmax_regression_vec(theta,inputData,labels,lambda);
        options.maxIter = 3;
        softmaxModel = softmaxTrain(d,k,lambda,inputData,labels,options);
        theta = softmaxModel.optTheta; %theta����
        upgradeTheta = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% ����2.2 �������Դ��� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:1:n
        if isLabeled(i) == 1            %%ʵ������Ѵ�����������ʵ��
            fa(i) = 100;
            f(i) = 100;
            continue;
        end
        
        tempIncomData = X(:,i);        % �����걸���ݼ��ĵ���ʵ����ֵ��tempIncomData���Ӷ����ı�ԭ���걸���ݼ�
        tempRandomX = randomX(:,i);
        
        predProbability(:,i) = theta'*tempIncomData;                           %%����z1��z2��z3
        predProbability(:,i) = exp(predProbability(:,i));                      %%����e^z1��e^z2��e^z3
        sumProb(i) = sum(predProbability(:,i));                                %%����e^z1+e^z2+e^z3
        %%������ֵĸ��ʺ���������
        completeP(i) = max(predProbability(:,i))./sumProb(i);                  %%�������������Ը���ֵ�����ȱʧ����ֵλ��
        attributeCost(i) = (1-completeP(i)).*misCostValue;
        fa(i) = attributeCost(i);                                              %%��¼
        
        ifMiss = ismember(tempRandomX,0);                                        %%���ʵ��û��ȱֵ������Ҫ�������Թ���,ֱ�Ӽ�����������
        if sum(ifMiss) == 0
            continue;
        end
        
        indexMiss = find(tempRandomX == 0);                                        %% ��¼tempIncomData������ȱʧ���Եĵ�ַ����
        numMiss = length(find(tempRandomX == 0));                                  %% ��¼tempIncomData������ȱʧ���Եĸ���
        countNum = 1;
        for j = 1:1:numMiss
            scheme = combntns(indexMiss,j);                                          %%ѭ�������ͬjʱindexMiss��ȱʧֵ��ַ���������
            numScheme = size(scheme,1);
            for jj = 1:1:numScheme
                tempIncomData(scheme(jj,:)) = averValue(scheme(jj,:));               %%�������������ʱ��averValue��ȱֵλ�õ�ֵ����tempIncomData����Ӧλ��
                %%����ÿһ��ĸ���
                predProbability(:,i) = theta'*tempIncomData;                         %%����z1��z2��z3
                predProbability(:,i) = exp(predProbability(:,i));                    %%����e^z1��e^z2��e^z3
                sumProb(i) = sum(predProbability(:,i));                              %%����e^z1+e^z2+e^z3
                
                %%������ֵĸ��ʺ���������
                updateP(i,countNum) = max(predProbability(:,i))./sumProb(i);         %%�������������Ը���ֵ�����ȱʧ����ֵλ��
                misclassification = (1-updateP(i,countNum)).*misCostValue;
                updateCost(i,countNum) = misclassification + j*attributeCostValue;   %%����������Թ��򷽰��µ��ܴ���
                
                if updateCost(i,countNum) < attributeCost(i)                         %%��¼�������Թ��򷽰��µ���С�ܴ���
                    attributeCost(i) = updateCost(i,countNum);
                    H(i) = {scheme(jj,:)};                                           %%����cell���飬���Ը���
                    fa(i) = attributeCost(i);                                        %%��¼
                end
                countNum = countNum + 1;
                tempIncomData = X(:,i);                                              %%��ԭ���걸ʵ��tempIncomData
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% ����2.3 �����ǩ���ۺ�ʵ������ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:1:n
        if isLabeled(i) == 1                                                               %% ʵ������Ѵ�����������ʵ��
            fl(i) = 100;
            f(i) = 100;
            continue;
        end
        %%%%%%%%% ��ǩ���ۼ��㷽ʽһ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    fl(i) = (1 - gamma(i)/maxGamma) * teachCostValue;
        
        %%%%%%%%% ��ǩ���ۼ��㷽ʽ�� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        c = d - 1;
        Wi(i) = length(find(randomX(:,i) == 0));                                              %% ÿһ��������ȱʧ���ݸ���
        b(i,:) = ((c - Wi(i))/c) * (gamma(i)/maxGamma);                                 %% ������
        fl(i) = (1 - b(i)) * teachCostValue;                                            %% ��ʦ����
        
        %%%%%%%%% ÿ��ʵ����ʵ�����ۣ�������С��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f(i) = min(fa(i), fl(i));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% ����2.4 ���Դ��ۺͻ�ȡ��ǩ���۱Ƚϣ�����ѵ���� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [valueMinF,indexMinF] = min(f);                                                 %% �õ�ʵ�������е���Сֵ������������ַ
    tempIncomData = X(:,indexMinF);
    tempOriData = OriX(:,indexMinF);
    if fa(indexMinF) < fl(indexMinF)
        if isempty(H{indexMinF}) == 1                                            %% ���cell����H�ж�Ӧλ��Ϊ�գ�����Ҫ��temprandomX��ֵ
        else
            tempIncomData(H{indexMinF},:) = tempOriData(H{indexMinF},:);         %% ��������
        end
        
        [~,pred]= max(theta'*tempIncomData);                                         %% theta'*data���������ǩ�ĸ��ʣ�������������������pred
        Labels(indexMinF) = pred;                                                %% ��Ԥ��ı�ǩpred�ŵ���ǩԤ��ľ�����
        isLabeled(indexMinF) = 1;
        
    else
        %%�������Դ��۴��ڽ�ʦ���ۣ������ǩ
        Labels(indexMinF) = Y(indexMinF);
        queryInstances = [queryInstances,indexMinF];
        H{indexMinF} = [];
        isLabeled(indexMinF) = 1;
        index = ismember(tempIncomData,0);
        if sum(index) == 0
            trainingSet = [trainingSet; indexMinF];
            labels = [labels;Y(indexMinF)];                                      %%�ڴ˴����±�ǩ�б�
            upgradeTheta = 1;
        end
    end
end

%%======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����3 ����׼ȷ�ʺ�ƽ������ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempInCorrect = 0;
totalCost = 0;
buyAttribute = 0;
for i = 1:1:n
    buyAttribute = buyAttribute + attributeCostValue*length(H{i});              %% �������Ե��ܴ���
    if (Labels(i) ~= Y(i))
        tempInCorrect = tempInCorrect + 1;
    end
end

buyLable = [buyLable',queryInstances];     % ��¼���й��������
queryNumbers = length(queryInstances)+ ratio;
totalCost = (length(queryInstances)+ ratio)*teachCostValue + buyAttribute + tempInCorrect * misCostValue;      %% �����ܴ���
count0 = length(find(randomX == 0));
missRate = count0/(n*(d-1));
accuracy = (n - tempInCorrect - queryNumbers)/(n - queryNumbers);
averCost = totalCost/n;
% result = [missRate,accuracy,totalCost,averCost,length(queryInstances)]

end