function M = BPCA_initmodel(y,q)
% model initialization for 
% Bayesian PCA with missing value estimation
% released version Nov. 26 2002

%parameter
%    y: the matrix with missing values(no empty rows)
%    q: the number of axis axes
%    M.mu : estimated mean row vector
%    M.W  : estimated principal axes matrix
%           for example M.W(:,1) is the 1st. principal axis vector.
%    M.tau : estimated precision (inverse variance) of
%            the residual error.
[N,d] = size(y);

M.N = N;
M.q = q;
M.d = d;

M.yest = y;
M.missidx = cell(N,1);
M.nomissidx = cell(N,1);
M.gnomiss = [];
M.gmiss = [];
for i=1:N
%   M.missidx{i} = find(y(i,:)>900);
%   M.nomissidx{i} = find(y(i,:)<900);
  M.missidx{i} = find(isnan(y(i,:)));
  M.nomissidx{i} = find(~isnan(y(i,:)));
  if length(M.missidx{i}) == 0 % a row that has no missing values
    M.gnomiss = [M.gnomiss i]; %Complete rows
  else
    M.gmiss = [M.gmiss i]; %incomplete rows
    M.yest(i,M.missidx{i}) = 0; %set missing values to zero
  end
end

ynomiss = y(M.gnomiss,:);

covy = cov(M.yest);

[U,S,V] = svds(covy,q);

M.mu = zeros(1,d);
for j=1:d
  %idx = find(y(:,j)<900);
  idx= find(~isnan(y(:,j)));
  M.mu(j) = mean(y(idx,j)); %mean of all no missing values in this column
end

M.W  = U * sqrt(S);
M.tau = 1/( sum(diag(covy)) -sum(diag(S)) );
taumax = 1e10;
taumin = 1e-10;
M.tau = max( min( M.tau, taumax), taumin );

M.galpha0 = 1e-10;
M.balpha0 = 1;
M.alpha = (2*M.galpha0 + d)./(M.tau*diag(M.W'*M.W)+2*M.galpha0/M.balpha0);

M.gmu0  = 0.001;

M.btau0 = 1;
M.gtau0 = 1e-10;
M.SigW = eye(q);
