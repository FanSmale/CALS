%function missing data setup
function M_complete=GetComplete(filename)
load(filename);

%get the complete matrix from the selected dataset
M_complete=zeros(0,size(GeneData,2));%the matirx name of the loaded data is 'GeneData'(PreSeted)
for i=1:size(GeneData,1)
    if iscomplete(GeneData(i,:))
        M_complete=[M_complete;(GeneData(i,:))];%assume that the data are log2 transformed before analysis
    end % if
end % for i

end %function