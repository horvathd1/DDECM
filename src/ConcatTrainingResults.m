function [res] = ConcatTrainingResults(experiment,result)
sepres=GetTrainingResults(experiment,result);
for i=1:length(sepres)
    for j=1:length(sepres{1})
        res{j}=cellfun(@(x)x{j},sepres,'UniformOutput',false);
        longdim=size(res{j}{1})==length(res{j}{1});
        if longdim(2)
            res{j}=[res{j}{:}];
        else
            res{j}=vertcat(res{j}{:});
        end
    end
end
end

