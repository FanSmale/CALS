%function to see whether InputM line is complete
function out=iscomplete(M)
n=0;
for i=1:size(M,2)
    if isnan(M(1,i))
        break;
    end %if
    n=n+1;
end %for i

if size(M,2)==n
    out=true;
else
    out=false;
end
end