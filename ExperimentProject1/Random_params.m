function [res] = Random_params(hyparams,monitor)
% Hyparams:
%   - n: number of param combinations
%   - seed: seed
%   - param_range
%   - tau_range
rng=RandStream("mt19937ar","Seed",hyparams.seed);
n=hyparams.n;
params=rng.rand(5,n)-0.5*[ones(4,1);0];
param_range=hyparams.param_range;
tau_range=hyparams.tau_range;
params=params.*[2*param_range*ones(4,1);tau_range];
anal_pls=AnalyticPL_DL(params(1,:),params(2,:),params(3,:),params(4,:),params(5,:));
models=Model.empty(0,n);
monitor.Metrics=["RelErr","MaxRelErr","AnalPL","ApproxPL"];
monitor.groupSubPlot("Err",["RelErr","MaxRelErr"]);
monitor.groupSubPlot("PL",["AnalPL","ApproxPL"]);
monitor.recordMetrics(1,RelErr=0,MaxRelErr=0,AnalPL=0,ApproxPL=0);
rel_errs=zeros(1,n);
for i=1:n
    param=params(:,i);
    models(i)=Inv_params(struct(a=param(1),b=param(2),c=param(3),d=param(4),tau=param(5),eq="DL"), ...
        [],toplot=false);
    disp("i"+string(models(i).PL))
    anal_pl=anal_pls(end,i);
    relerr=abs((models(i).PL-anal_pl)/anal_pl)*100;
    maxrelerr=monitor.MetricData.MaxRelErr(end,2);
    if relerr>maxrelerr
        maxrelerr=relerr;
    end
    monitor.recordMetrics(i+1,RelErr=relerr,MaxRelErr=maxrelerr,AnalPL=anal_pl,ApproxPL=models(i).PL);
    rel_errs(i)=relerr;
    monitor.Progress=i/n*100;
    if ~isempty(monitor)&&monitor.Stop
        break
    end
end
res={models,anal_pls,rel_errs};
end

