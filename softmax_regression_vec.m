function [f,g] = softmax_regression_vec(theta, X,y,lambda )  
%这里的n指数据有多维（包括新加的x0=1这一维，也就是θj向量的维度）就是inputSize
%这里y是1,2....k，从1开始的 
  % Arguments:  
  %   theta - A vector containing the parameter values to optimize.  
  %       In minFunc, theta is reshaped to a long vector.  So we need to  
  %       resize it to an n-by-(num_classes-1) matrix.  
  %       Recall that we assume theta(:,num_classes) = 0.  
  m=size(X,2);%X每一列是一个样本，m是指有m个样本  
  n=size(X,1);  %n指代的前面说了
  theta=reshape(theta, n, []); %也就是把theta设置成这样矩阵：inputSize行也就是n行，每一列是一个θj，有k列。
  % initialize objective value and gradient.  
  f = 0;  
  g = zeros(size(theta));  
  h = theta'*X;%h是k行m列的矩阵，见图1.
  a = exp(h);  
  p = bsxfun(@rdivide,a,sum(a)); % sum(a)是一个行向量，每个元素是a矩阵的每一列的和。然后运用bsxfun(@rdivide，，）
  %是a矩阵的第i列的每个元素除以 sum(a)向量的第i个元素。得到的p矩阵大小和图1一样，每个元素如图2.
  c = log(p); %然后我们取log2的对数，c矩阵大小和图1一样，每个元素如图3
  i = sub2ind(size(c), y',1:size(c,2)); %y',1:size(c,2)这两个向量必须同时是行向量或列向量
  %因为我们接下来每一个样本xi对应的yi是几，就去找到p的每一列中，所对应的第几个元素就是要找的，如图4.首先使用sub2ind
  %sub2ind: 在matlab中矩阵是按一列一列的存储的，比如A=[1 2 3;4 5 6]
%那么A（2）=4，A(3)=2...而这个函数作用就是比如 sub2ind（size（A）,2,1）就是返回A的第2行第一列的元素存储的下标，因为
%A（2）=4，所以存储的下标是2，所以这里返回2.这里sub2ind（size（A）,2,1）的2,1也可以换成向量[a1,a2..],[b1,b2..]但是注意
%这两个向量必须同时是行向量或列向量，而不能一个是行向量一个是列向量。所以返回的
%第一个元素是A的第a1行第b1列的元素存储的下标，返回的第,二个元素是A的第a2行第b2列的元素存储的下标...i是一个向量，c（i）得到的
%向量的每一个元素就是p中每一列你前面要找的的元素。
  values = c(i);  
  f = -(1/m)*sum(values)+ lambda/2 * sum(theta(:) .^ 2);  %这个就是cost function 
  d = full(sparse(1:m,y,1)); %d为一个稀疏矩阵，有m行k列（k是类别的个数），这个矩阵的（1，y（1））、（2，y（2））
  %....(m,y(m))位置都是1。
  % l = size(theta,2) - max(y);
  if size(d,2)<size(theta,2)
      d(:,(size(d,2)+1):size(theta,2)) = 0;
  end
  g = (1/m)*X*(p'-d)+ lambda * theta; %这个g和theta矩阵的结构一样。 
  g=g(:); % 再还原成向量的形式，这里（：）和reshape都是按列进行的，所以里面位置并没有改变。