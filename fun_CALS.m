function [missRate,accuracy,totalCost,averCost,queryNumbers,buyLable]=fun_CALS(data,dataOri,randomX,dc,ratio,distanceSelection,normalization)

%% 读取缺失数据
[n, d] = size(data);  %数据集n行，d维
X = data(:,1:d-1);  %%缺失无标签数据
Y = data(:,d);      %%真实标签

%% 读取原始数据
OriX = dataOri(:,1:d-1);  %%完整无标签数据
oneArray = ones(1,n);
OriX = [oneArray; OriX']; %每个样本都要增加一个x0=1;

%% 初始化参数
H = cell(1,n);      %%初始化一个空cell数组
labels = [];
queryInstances = [];

%% 模型参数
k = numel(unique(Y));           % 统计标签种类数
misCostValue = 2;  % 设置误分类代价
teachCostValue = 1;     % 设置购买标签的教师代价
attributeCostValue = 0.2;     % 设置购买属性的代价
lambda = 1e-4;         % Weight decay parameter

upgradeTheta = 0;

%%======================================================================
%% 初始化
for i = 1:1:n
    isLabeled(i) = 0;           %%对已处理的实例进行标记
    instanceCost(i) = 0;        %%购买标签的代价
    Labels(i) = 0;              %%记录实例的标签信息
    minallCost(i) = 0;          %%初始化记录最小总代价的矩阵
    attributeCost(i) = 200;     %%初始化属性代价为加大值
end

%% 距离计算
switch distanceSelection
    case 1
        Dists = pdist(X,'cityblock');  %% 曼哈顿距离
    case 2
        Dists = pdist(X,'cosine');  %% 夹角余弦距离
    case 3
        Dists = pdist(X,'euclidean');  %%欧几里德距离
    case 4
        Dists = pdist(X,'minkowski');  %% 闵科夫斯基距离
    case 5
        Dists = pdist(X,'hamming');  %% 汉明距离
    case 6
        Dists = pdist(X,'seuclidean');  %%标准欧几里德距离
    case 7
        Dists = pdist(X,'mahalanobis');  %%马哈拉诺比斯距离Mahalanobis distance 马氏距离
    case 8
        Dists = pdist(X,'chebychev');  %% 切比雪夫距离
    otherwise
        Dists = pdist(X,'jaccard');  %% 杰卡德距离
end
Dists = squareform(Dists);
maxDist = max(max(Dists));

%% 计算代表性，也就是计算密度和距离
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
[~, desInd] = sort(gamma, 'descend');               %%%%%%%desInd为按降序排列的索引
maxGamma = max(gamma);

%% 修正数据集,归一化设置
X = [oneArray; X']; %每个样本都要增加一个x0=1;
if normalization == 1
    X = mapminmax(X);
    X = X.*randomX;
    OriX = mapminmax(OriX);
end

%% 计算属性加权平均值
for i = 2:1:d
    nonAttributeNumbers(i) = length(nonzeros(X(i,:)));  %%%%%%%%%每个属性非缺值元素个数
end

for i = 2:1:d
    averValue(i-1) = sum(X(i,:))./nonAttributeNumbers(i);
end
averValue = averValue';
averValue = [1;averValue];

%%======================================================================
%%%%%%%%%%% 步骤1：构建初始训练集，训练theta
%%%%%%%%%%% 构建初始训练集 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%通过选择代表点的方式产生初始训练集 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trainingSet = desInd(1:1:ratio);
buyLable = trainingSet;
trainingInput = X(:,trainingSet);     %%%%%%%训练集
labels = Y(trainingSet);              %%%%%%%训练集标签
Labels(trainingSet) = Y(trainingSet);
Labels = Labels';
isLabeled(trainingSet) = 1;

%%%%%%%%%%% 训练theta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = 0.0005 * randn(d, k);     % 随机产生初始theta
[cost, grad] = softmax_regression_vec(theta,trainingInput,labels,lambda);
options.maxIter = 100;
softmaxModel = softmaxTrain(d, k, lambda, trainingInput,labels, options);
theta = softmaxModel.optTheta;  %重新获得theta矩阵

%%======================================================================
%%%%%%%%%%% 步骤2：主循环 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while sum(isLabeled) ~= n
    
    %%%%%% 步骤2.1 训练集更新，重新训练模型 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if upgradeTheta == 1
        inputData = X(:,trainingSet);
        [cost, grad] = softmax_regression_vec(theta,inputData,labels,lambda);
        options.maxIter = 3;
        softmaxModel = softmaxTrain(d,k,lambda,inputData,labels,options);
        theta = softmaxModel.optTheta; %theta矩阵
        upgradeTheta = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% 步骤2.2 计算属性代价 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:1:n
        if isLabeled(i) == 1            %%实例如果已处理，则跳过该实例
            fa(i) = 100;
            f(i) = 100;
            continue;
        end
        
        tempIncomData = X(:,i);        % 将不完备数据集的单个实例赋值给tempIncomData，从而不改变原不完备数据集
        tempRandomX = randomX(:,i);
        
        predProbability(:,i) = theta'*tempIncomData;                           %%计算z1，z2，z3
        predProbability(:,i) = exp(predProbability(:,i));                      %%计算e^z1，e^z2，e^z3
        sumProb(i) = sum(predProbability(:,i));                                %%计算e^z1+e^z2+e^z3
        %%计算出现的概率和误分类代价
        completeP(i) = max(predProbability(:,i))./sumProb(i);                  %%将计算所得属性概率值存放在缺失属性值位置
        attributeCost(i) = (1-completeP(i)).*misCostValue;
        fa(i) = attributeCost(i);                                              %%记录
        
        ifMiss = ismember(tempRandomX,0);                                        %%如果实例没有缺值，则不需要进行属性购买,直接计算误分类代价
        if sum(ifMiss) == 0
            continue;
        end
        
        indexMiss = find(tempRandomX == 0);                                        %% 记录tempIncomData矩阵中缺失属性的地址索引
        numMiss = length(find(tempRandomX == 0));                                  %% 记录tempIncomData矩阵中缺失属性的个数
        countNum = 1;
        for j = 1:1:numMiss
            scheme = combntns(indexMiss,j);                                          %%循环输出不同j时indexMiss中缺失值地址的排列组合
            numScheme = size(scheme,1);
            for jj = 1:1:numScheme
                tempIncomData(scheme(jj,:)) = averValue(scheme(jj,:));               %%将各种排列组合时，averValue中缺值位置的值赋给tempIncomData中相应位置
                %%计算每一类的概率
                predProbability(:,i) = theta'*tempIncomData;                         %%计算z1，z2，z3
                predProbability(:,i) = exp(predProbability(:,i));                    %%计算e^z1，e^z2，e^z3
                sumProb(i) = sum(predProbability(:,i));                              %%计算e^z1+e^z2+e^z3
                
                %%计算出现的概率和误分类代价
                updateP(i,countNum) = max(predProbability(:,i))./sumProb(i);         %%将计算所得属性概率值存放在缺失属性值位置
                misclassification = (1-updateP(i,countNum)).*misCostValue;
                updateCost(i,countNum) = misclassification + j*attributeCostValue;   %%计算各种属性购买方案下的总代价
                
                if updateCost(i,countNum) < attributeCost(i)                         %%记录各种属性购买方案下的最小总代价
                    attributeCost(i) = updateCost(i,countNum);
                    H(i) = {scheme(jj,:)};                                           %%采用cell数组，可以更改
                    fa(i) = attributeCost(i);                                        %%记录
                end
                countNum = countNum + 1;
                tempIncomData = X(:,i);                                              %%还原不完备实例tempIncomData
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% 步骤2.3 计算标签代价和实例代价 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:1:n
        if isLabeled(i) == 1                                                               %% 实例如果已处理，则跳过该实例
            fl(i) = 100;
            f(i) = 100;
            continue;
        end
        %%%%%%%%% 标签代价计算方式一 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    fl(i) = (1 - gamma(i)/maxGamma) * teachCostValue;
        
        %%%%%%%%% 标签代价计算方式二 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        c = d - 1;
        Wi(i) = length(find(randomX(:,i) == 0));                                              %% 每一个样本中缺失数据个数
        b(i,:) = ((c - Wi(i))/c) * (gamma(i)/maxGamma);                                 %% 收益率
        fl(i) = (1 - b(i)) * teachCostValue;                                            %% 教师代价
        
        %%%%%%%%% 每个实例的实例代价（两者最小）%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f(i) = min(fa(i), fl(i));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% 步骤2.4 属性代价和获取标签代价比较，更新训练集 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [valueMinF,indexMinF] = min(f);                                                 %% 得到实例代价中的最小值，及其索引地址
    tempIncomData = X(:,indexMinF);
    tempOriData = OriX(:,indexMinF);
    if fa(indexMinF) < fl(indexMinF)
        if isempty(H{indexMinF}) == 1                                            %% 如果cell矩阵H中对应位置为空，则不需要给temprandomX赋值
        else
            tempIncomData(H{indexMinF},:) = tempOriData(H{indexMinF},:);         %% 购买属性
        end
        
        [~,pred]= max(theta'*tempIncomData);                                         %% theta'*data求得三个标签的概率，将概率最大的索引赋给pred
        Labels(indexMinF) = pred;                                                %% 将预测的标签pred放到标签预测的矩阵中
        isLabeled(indexMinF) = 1;
        
    else
        %%购买属性代价大于教师代价，购买标签
        Labels(indexMinF) = Y(indexMinF);
        queryInstances = [queryInstances,indexMinF];
        H{indexMinF} = [];
        isLabeled(indexMinF) = 1;
        index = ismember(tempIncomData,0);
        if sum(index) == 0
            trainingSet = [trainingSet; indexMinF];
            labels = [labels;Y(indexMinF)];                                      %%在此处更新标签列表
            upgradeTheta = 1;
        end
    end
end

%%======================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 步骤3 计算准确率和平均代价 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempInCorrect = 0;
totalCost = 0;
buyAttribute = 0;
for i = 1:1:n
    buyAttribute = buyAttribute + attributeCostValue*length(H{i});              %% 购买属性的总代价
    if (Labels(i) ~= Y(i))
        tempInCorrect = tempInCorrect + 1;
    end
end

buyLable = [buyLable',queryInstances];     % 记录所有购买的样本
queryNumbers = length(queryInstances)+ ratio;
totalCost = (length(queryInstances)+ ratio)*teachCostValue + buyAttribute + tempInCorrect * misCostValue;      %% 计算总代价
count0 = length(find(randomX == 0));
missRate = count0/(n*(d-1));
accuracy = (n - tempInCorrect - queryNumbers)/(n - queryNumbers);
averCost = totalCost/n;
% result = [missRate,accuracy,totalCost,averCost,length(queryInstances)]

end