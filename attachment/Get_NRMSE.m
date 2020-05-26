function out_NRMSE=Get_NRMSE(Mcomplete,Mmissing,Mimputed)
% NRMSE(normalised root mean squared error) Calculator
% Mcomplete: the original complete matrix
% Mmissing:  the artificial matrix that has been removed some entries as missing values
% Mimputed:  the imputed matrix

if size(Mcomplete)~=size(Mmissing)|size(Mcomplete)~=size(Mimputed)|size(Mcomplete)~=size(Mimputed)
    fprintf('dimensions do not match\r\n');
    return;
end %if

idx=find(isnan(Mmissing));
Guess=Mimputed(idx);
Answer=Mcomplete(idx);
out_NRMSE=sqrt(mean((Guess-Answer).^2))/std(Answer);
end %function