function [f,g] = softmax_regression_vec(theta, X,y,lambda )  
%�����nָ�����ж�ά�������¼ӵ�x0=1��һά��Ҳ���Ǧ�j������ά�ȣ�����inputSize
%����y��1,2....k����1��ʼ�� 
  % Arguments:  
  %   theta - A vector containing the parameter values to optimize.  
  %       In minFunc, theta is reshaped to a long vector.  So we need to  
  %       resize it to an n-by-(num_classes-1) matrix.  
  %       Recall that we assume theta(:,num_classes) = 0.  
  m=size(X,2);%Xÿһ����һ��������m��ָ��m������  
  n=size(X,1);  %nָ����ǰ��˵��
  theta=reshape(theta, n, []); %Ҳ���ǰ�theta���ó���������inputSize��Ҳ����n�У�ÿһ����һ����j����k�С�
  % initialize objective value and gradient.  
  f = 0;  
  g = zeros(size(theta));  
  h = theta'*X;%h��k��m�еľ��󣬼�ͼ1.
  a = exp(h);  
  p = bsxfun(@rdivide,a,sum(a)); % sum(a)��һ����������ÿ��Ԫ����a�����ÿһ�еĺ͡�Ȼ������bsxfun(@rdivide������
  %��a����ĵ�i�е�ÿ��Ԫ�س��� sum(a)�����ĵ�i��Ԫ�ء��õ���p�����С��ͼ1һ����ÿ��Ԫ����ͼ2.
  c = log(p); %Ȼ������ȡlog2�Ķ�����c�����С��ͼ1һ����ÿ��Ԫ����ͼ3
  i = sub2ind(size(c), y',1:size(c,2)); %y',1:size(c,2)��������������ͬʱ����������������
  %��Ϊ���ǽ�����ÿһ������xi��Ӧ��yi�Ǽ�����ȥ�ҵ�p��ÿһ���У�����Ӧ�ĵڼ���Ԫ�ؾ���Ҫ�ҵģ���ͼ4.����ʹ��sub2ind
  %sub2ind: ��matlab�о����ǰ�һ��һ�еĴ洢�ģ�����A=[1 2 3;4 5 6]
%��ôA��2��=4��A(3)=2...������������þ��Ǳ��� sub2ind��size��A��,2,1�����Ƿ���A�ĵ�2�е�һ�е�Ԫ�ش洢���±꣬��Ϊ
%A��2��=4�����Դ洢���±���2���������ﷵ��2.����sub2ind��size��A��,2,1����2,1Ҳ���Ի�������[a1,a2..],[b1,b2..]����ע��
%��������������ͬʱ������������������������һ����������һ���������������Է��ص�
%��һ��Ԫ����A�ĵ�a1�е�b1�е�Ԫ�ش洢���±꣬���صĵ�,����Ԫ����A�ĵ�a2�е�b2�е�Ԫ�ش洢���±�...i��һ��������c��i���õ���
%������ÿһ��Ԫ�ؾ���p��ÿһ����ǰ��Ҫ�ҵĵ�Ԫ�ء�
  values = c(i);  
  f = -(1/m)*sum(values)+ lambda/2 * sum(theta(:) .^ 2);  %�������cost function 
  d = full(sparse(1:m,y,1)); %dΪһ��ϡ�������m��k�У�k�����ĸ��������������ģ�1��y��1��������2��y��2����
  %....(m,y(m))λ�ö���1��
  % l = size(theta,2) - max(y);
  if size(d,2)<size(theta,2)
      d(:,(size(d,2)+1):size(theta,2)) = 0;
  end
  g = (1/m)*X*(p'-d)+ lambda * theta; %���g��theta����Ľṹһ���� 
  g=g(:); % �ٻ�ԭ����������ʽ�����������reshape���ǰ��н��еģ���������λ�ò�û�иı䡣