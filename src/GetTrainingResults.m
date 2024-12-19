function [res] = GetTrainingResults(experiment,result)
dirname="ExperimentProject1/Results/";
dirname=dirname+dir(dirname+experiment+"_Result"+result+"*").name+"/";
res={};
list=dir(dirname+"T*");
for i=1:length(list)
    res{end+1}=load(dirname+list(i).name+"/output.mat").outputs;
end
end

