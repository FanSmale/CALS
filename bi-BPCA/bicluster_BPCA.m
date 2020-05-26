% The missing values in matrices are represented by 'NaN' in the Matlab .mat
% file. Before running this method, all missing values in your data should
% be converted as the form of 'NaN'.

% Folders attachment, bi-BPCA and BPCA should be added into current working
% directory

% M. Fanchi, C. Cheng, and Y. Hong, "A Bicluster-Based Bayesian Principal
% Component Analysis Method for Microarray Missing Value Estimation," 
% Biomedical and Health Informatics, IEEE Journal of, vol. 18, pp. 863-871, 2014.

%iterative BPCA with bicluster
function [E]  =  bicluster_BPCA(M_missing,k,T0,itrs)
[m,~] = size(M_missing);
global M_temp;
old_E = M_missing;

%begin iteration
iteration = itrs;
while true
    for i = 1:m % for every row in the original gene
        if ~iscomplete(M_missing(i,:))
            fprintf('dealing with gene %0.0f\n',i);
            missidxj = find(isnan(M_missing(i,:))); len_miss = length(missidxj);
            nomissidxj = find(~isnan(M_missing(i,:)));  len_nomiss = length(nomissidxj);
            M_others = M_temp; %M_temp is not full
            if size(M_temp,1) == size(M_missing,1) %M_temp is full
                M_others(i,:) = [];
            end
            
            [A,B,w]  =  similargene(M_missing,i,missidxj,nomissidxj);
            A = A(1:k,:);B = B(1:k,:);
            R = B'*A;
            %calculate the weighted Enclidean distance for the jth missing value in
            %the target gene
            dises = zeros(size(M_others,1),len_miss);
            for miss_j = 1:len_miss %for the jth missing value in the target gene
                for neighbor_s = 1:size(M_others,1)
                    dises(neighbor_s,miss_j) = distance4j(M_others,len_miss,R,w,nomissidxj,miss_j,neighbor_s);
                end %for neighbor_s
                [~,sort_idx] = sort(dises(:,miss_j));
                neighbors4missj = M_others(sort_idx,:);
                correlated_js = findCorrelated(R,miss_j,T0);
                w_j = w(:,correlated_js);
                Aj = neighbors4missj(:,nomissidxj(correlated_js));Aj = Aj(1:k,:);
                Bj = neighbors4missj(:,missidxj(miss_j));Bj = Bj(1:k,:);
                alpha = M_missing(i,missidxj(miss_j));
                bicluster = [alpha,w_j;Bj,Aj];
                complete_bicluster = BPCAfill(bicluster);
                old_E(i,missidxj(miss_j)) = complete_bicluster(1,1);
            end % for miss_j
        end %if ~iscomplete
    end%i
    
    if iteration == itrs
        E = old_E;
        return;
    else
        M_temp = old_E;
        iteration = iteration + 1;
    end %if
end %while
end %function

function [A,B,w]  =  similargene(M_missing,i,missidxj,nomissidxj)
global M_temp;
global gene_include;
global m_gene_include;
% L2-norm distance calculation
mm1 = 1;  mm2 = m_gene_include;
AA1 = M_temp(i,nomissidxj);
BB1 = M_temp(gene_include,nomissidxj);
AA2 = sum(AA1.^2,2);
BB2 = sum(BB1.^2,2);
distance = repmat(AA2,1,mm2)+repmat(BB2',mm1,1)-2*AA1*BB1';

[~,sortidx] = sort(distance);
% gene number
gene = gene_include(sortidx);

% remove the same gene
gene(gene == i) = [];

A = M_temp(gene,nomissidxj);
B = M_temp(gene,missidxj);
w = M_missing(i,nomissidxj);

end

function dis = distance4j(M_others,len_miss,R,w,nomissidxj,miss_j,neighbor_s)
fraction = 0;denominator = 0;
for v = 1:size(M_others,2)-len_miss
    fraction = fraction+R(miss_j,v)^2*(w(1,v)-M_others(neighbor_s,nomissidxj(v)))^2;
    denominator = denominator+R(miss_j,v)^2;
end %for v
dis = sqrt(fraction/denominator);
end

function js = findCorrelated(R,miss_j,T0)
js = [];
%T0 = 1e-4;
for j = 1:size(R,2)
    if abs(R(miss_j,j))>= T0*max(abs(R(miss_j,:)));
        js = [js,j];
    end %if
end %for j
end