function M_missing = GetMissing(M_complete,mRate,para)
% mRate is the percentage of the missing value,for example,mRate=10 makes a
% missing matrix with 10% missing values.para is the style of the missing
% values.'P1' means get mRate% missing entries by random; 'P2' means get
% mRate continuous missing entries by random.
switch para
    % P1
    case 'P1'
        targetNums=size(M_complete,1)*size(M_complete,2)*mRate/100;
        totalMissingNums=0;
        
        while totalMissingNums<=targetNums
            
            randline=fix(rand(1)*size(M_complete,1)+1);
            randcolumn=fix(rand(1)*size(M_complete,2)+1);
            M_complete(randline,randcolumn)=NaN;
            totalMissingNums=totalMissingNums+1;
            
        end %while
        M_missing=M_complete;
        
        %P2
    case 'P2'
        targetNums=size(M_complete,1)*size(M_complete,2)*mRate/100;
        totalMissingNums=0;
        while totalMissingNums<=targetNums
            
            randline=fix(rand(1)*size(M_complete,1)+1);
            randcolumn=fix(rand(1)*size(M_complete,2)+1);
            randlen=fix(rand(1)*mRate/100*size(M_complete,2)+1);
            if randcolumn+randlen<=size(M_complete,2)
                M_complete(randline,randcolumn:randcolumn+randlen)=NaN;
            end %if
            totalMissingNums=totalMissingNums+randlen+1;
            
        end% while
        M_missing=M_complete;
end
end