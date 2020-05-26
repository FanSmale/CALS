%T_estimate
%kestimate
function T = T_estimate(M_missing,k)
global M_temp;

[m,n] = size(M_missing);
[miss_gene,~] = find(isnan(M_missing));

Ts = [0,1e-6,1e-5,1e-4,1e-3,1e-2,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1];

% determine pretending missing positions p.
% estimate the number of missing positions per line.
rho = length(find(isnan(M_missing)))/m;
% rho=5 --> e=1, rho=10 --> e=2, ...
e = floor(rho/5);
% at least 1 position should be pretended as a missing value.
e = max([1 e]);
% the number of pretended missing values should be less than
% 20% of non-missing positions.....
e = min([e n*0.2*(1-0.2)]);
% pretended missing positions among the non-missing positions.
if e<0.5
    e = 0.5;
end
p = [1:2:2*e];

answer = [];
%k = k_min; kcount=0;
% while(k < k_max)
%     kcount=kcount+1; guess{kcount}=[];  k = ceil(k*k_jump);
% end
guess{length(Ts)} = [];

%fprintf(fid, 'Estimating missing values to determine k-value....\n');
for i = 1:m
    
    missidxj = find(isnan(M_missing(i,:))); len_miss = length(missidxj);
    nomissidxj = find(~isnan(M_missing(i,:)));  len_nomiss = length(nomissidxj);
    M_others = M_temp; %M_temp is not full
    if size(M_temp,1) == size(M_missing,1) %M_temp is full
        M_others(i,:) = [];
    end
    % determine the number of artifical missing values among p
    missnum = min([len_miss length(p)]);
    % too small nomissing entries to make artificial missing entries
    if (len_nomiss < 10)
        missnum = 1;
    end
    p1 = p(1:missnum);
    
    % generate the artificial missing values
    missidxj = nomissidxj(p1); nomissidxj(p1) = [];len_miss = length(missidxj);len_nomiss = length(nomissidxj);
    [A,B,w] = similargene_T(i,missidxj,nomissidxj);
    answer = [answer; M_missing(i,missidxj)'];
    R = B'*A;
    
    dises = zeros(size(M_others,1),len_miss);
    for miss_j = 1:len_miss %for the jth missing value in the target gene
        for neighbor_s = 1:size(M_others,1)
            dises(neighbor_s,miss_j) = distance4j_T(M_others,len_miss,R,w,nomissidxj,miss_j,neighbor_s);
        end %for neighbor_s
        [~,sort_idx] = sort(dises(:,miss_j));
        neighbors4missj = M_others(sort_idx,:);
        %for T survey
        Tcount = 0;
        for tt = Ts
            fprintf('surveying..., i=%f, T=%f\n',i,tt);
            Tcount = Tcount+1;
            correlated_js = findCorrelated_T(R,miss_j,tt);
            w_j = w(:,correlated_js);
            Aj = neighbors4missj(:,nomissidxj(correlated_js));Aj=Aj(1:k,:);
            Bj = neighbors4missj(:,missidxj(miss_j));Bj=Bj(1:k,:);
            Sub = [nan(1,size(Bj,2)),w_j;
                Bj,Aj];
            transFlag  =  0;
            if(size(Sub,1) < size(Sub,2)) % for BPCA, it is preferred that size(Sub,1) > size(Sub,2)
                Sub  =  Sub';
                transFlag  =  1;
            end
            temp = BPCAfill(Sub);
            if (transFlag == 1)
                temp  =  temp';
            end
            guess{Tcount} = [guess{Tcount}; temp(1,1:size(Bj,2))'];
        end %for tt
    end % for miss_j
end%i

Tcount = 0;
for tt = Ts
    Tcount = Tcount+1;
    nrmse(Tcount) = sqrt( mean( (guess{Tcount}-answer).^2 ) ) / std( answer );
end %for tt

% try to find global minimum
[minnrmse,idx] = min(nrmse);
T = Ts(idx);

function [A,B,w] = similargene_T(i,missidxj,nomissidxj)
%global M_temp;
global gene_include;
global m_gene_include;
% L2-norm distance calculation
mm1 = 1;
mm2 = m_gene_include;

AA1 = M_temp(i,nomissidxj);
BB1 = M_temp(gene_include,nomissidxj);
AA2 = sum(AA1.^2,2); BB2 = sum(BB1.^2,2);
distance = repmat(AA2,1,mm2) + repmat(BB2',mm1,1)-2*AA1*BB1';
[~,sortidx] = sort(distance);
% gene number
gene=gene_include(sortidx);
% remove the same gene
gene(gene==i) = [];

A = M_temp(gene,nomissidxj);
B = M_temp(gene,missidxj);
w = M_temp(i,nomissidxj);
end

function dis = distance4j_T(M_others,len_miss,R,w,nomissidxj,miss_j,neighbor_s)
fraction = 0;denominator = 0;
for v = 1:length(nomissidxj)
    fraction = fraction+R(miss_j,v)^2*(w(1,v) - M_others(neighbor_s,nomissidxj(v)))^2;
    denominator = denominator+R(miss_j,v)^2;
end %for v
dis = sqrt(fraction/denominator);
end

function js = findCorrelated_T(R,miss_j,T0)
js = [];
for j = 1:size(R,2)
    if abs(R(miss_j,j)) >= T0*max(abs(R(miss_j,:)));
        js = [js,j];
    end %if
end %for j
end

end% T_estimate