% The missing values in matrices are represented by 'NaN' in the Matlab .mat
% file. Before running this method, all missing values in your data should
% be converted as the form of 'NaN'.

% Folders attachment, bi-BPCA and BPCA should be added into current working
% directory

% M. Fanchi, C. Cheng, and Y. Hong, "A Bicluster-Based Bayesian Principal
% Component Analysis Method for Microarray Missing Value Estimation," 
% Biomedical and Health Informatics, IEEE Journal of, vol. 18, pp. 863-871, 2014.

function E = f_bicluster_BPCA(M_missing,itrs)
fprintf('Estimating missing values by bi-BPCA. The survey of parameters may take a long time.\n');
global M_temp;
global gene_include;
global m_gene_include;
[m,~] = size(M_missing);
[miss_gene,~] = find(isnan(M_missing));

%Find neighbors in BPCA imputed matrix or in the complete part of the
%original matrix
if m-length(miss_gene)<400
    %too many incomplete genes(rows),impute  the matrix by row average initially for LLS
    fprintf('Impute the matrix by BPCA initially\n');
    M_temp = BPCAfill(M_missing);
    gene_include  =  1:m;
    m_gene_include = length(gene_include);
    %fprintf('size of M_temp %0.0f*%0.0f\n',size(M_temp,1),size(M_temp,2));
    k = kestimate(M_missing);
    fprintf('the best k is %0.0f\n',k);
    
else
    fprintf('When finding neighbors,ignore genes that have missing values\n');
    M_temp = M_missing;
    gene_include = setdiff(1:m, miss_gene);
    m_gene_include = length(gene_include);
    k = kestimate(M_missing);
    fprintf('the best k is %0.0f\n',k);
    
end %if

T0 = T_estimate(M_missing,k);
% T0 = 0;
E = bicluster_BPCA(M_missing,k,T0,itrs);
end