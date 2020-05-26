%kestimate
function [mink]  =  kestimate(M_missing)
global m_gene_include;
[m,n] = size(M_missing);
[miss_gene,~] = find(isnan(M_missing));

%m_gene_include = length(gene_include);
k_max = m_gene_include;
%k_min = ceil(1/4*size(M_missing,1));
k_min = 10;
k_jump = 1.5;


% determine pretending missing positions p.
% estimate the number of missing positions per line.
rho = length(find(isnan(M_missing)))/m;
% rho = 5 --> e = 1, rho = 10 --> e = 2, ...
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
k  =  k_min; kcount = 0;
while(k < k_max)
    kcount = kcount+1; guess{kcount} = [];  %k  =  ceil(k*k_jump);
    k = k+ceil(m/20);
end

% fprintf(fid, 'Estimating missing values to determine k-value....\n');
for i = 1:m
    
    missidxj = find(isnan(M_missing(i,:))); len_miss = length(missidxj);
    nomissidxj = find(~isnan(M_missing(i,:)));  len_nomiss = length(nomissidxj);
    %if (len_nomiss < minexp || len_nomiss < 2)
    %            fprintf('%dth gene: skip due to nomiss_exp(%d)<%.2f or < 2\n', ...
    %                i,len_nomiss,minexp);
    %      elseif (len_miss > 0)
    
    % determine the number of artifical missing values among p
    missnum = min([len_miss length(p)]);
    % too small nomissing entries to make artificial missing entries
    if (len_nomiss < 10)
        missnum = 1;
    end
    p1 = p(1:missnum);
    %          if fig == 1
    %fprintf('%dth gene: apply llsq_l2 --- %d --> %d missings\n',i, len_miss, length(p1));
    %end
    % generate the artificial missing values
    missidxj = nomissidxj(p1); nomissidxj(p1) = [];
    [A,B,w]  =  similargene_k(i,missidxj,nomissidxj);
    answer = [answer; M_missing(i,missidxj)'];
    
    % for k survey
    k  =  k_min; kcount = 0;
    while(k < k_max)
        fprintf('estimating k = %d\n',k);
        kcount = kcount+1;
        Apart = A(1:k,:); Bpart = B(1:k,:);
        Sub = [nan(1,size(Bpart,2)),w;
            Bpart,Apart];
        transFlag  =  0;
        if(size(Sub,1) < size(Sub,2))
            Sub  =  Sub';
            transFlag  =  1; % for BPCA, it is preferred that size(Sub,1) > size(Sub,2)
        end
        temp = BPCAfill(Sub);
        if (transFlag == 1) 
            temp  =  temp';
        end
        guess{kcount} = [guess{kcount}; temp(1,1:size(Bpart,2))'];
        %k  =  ceil(k*k_jump);
        k = k+ceil(m/20);
    end %while k
    %end%if
    
end%i

k  =  k_min; kcount = 0; xk = []; 
mincount = 0;
while(k < k_max)
    kcount = kcount+1;
    nrmse(kcount)  =  sqrt( mean( (guess{kcount}-answer).^2 ) ) / std( answer );
    %fprintf(fid, 'nrmsw(k = %d): %f\n', k, nrmse(kcount));
    if(k > n) && (mincount == 0)
        mincount = kcount;
    end
    xk = [xk k];  %k  =  ceil(k*k_jump);
    k = k+ceil(m/20);
end

% try to find global minimum
[minnrmse,minkidx] = min(nrmse); mink = xk(minkidx);
%end %function


function [A,B,w]  =  similargene_k(i,missidxj,nomissidxj)
% L2-norm distance calculation
global M_temp;
global gene_include;
%global m_gene_include;
mm1 = 1;  mm2 = m_gene_include;

AA1 = M_temp(i,nomissidxj);
BB1 = M_temp(gene_include,nomissidxj);
AA2 = sum(AA1.^2,2); BB2 = sum(BB1.^2,2);
distance = repmat(AA2,1,mm2)+repmat(BB2',mm1,1)-2*AA1*BB1';
[~,sortidx] = sort(distance);
% gene number
gene = gene_include(sortidx);
% remove the same gene
gene(gene == i) = [];

A = M_temp(gene,nomissidxj);
B = M_temp(gene,missidxj);
w = M_temp(i,nomissidxj);
end

end %kestimate