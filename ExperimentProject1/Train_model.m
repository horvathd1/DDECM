function [mdl] = Train_model(hyparams,monitor,varargin)
%MULTIICTRAINING Train the model 
%   params need to contain the following:
%   - zeta,j: Parameters of the DDE
%       or a,b,c,d,tau
%   - ic: value of the ic to train for (in case of DL the divisor)
%   - num_per: number of periods to calculate the time series for
%   - num_steps: number of timesteps in a period
%   - to_cut: The first 100*to_cut% of the time series will be discarded
%   - reltol: The relative tolerance for the numeric integral
%   - eq: Name of the equation (DL or MT)
%   Varargin: toplot
p=inputParser;
p.addParameter("toplot",true);
p.parse(varargin{:});
toplot=p.Results.toplot;
if hyparams.eq=="MT"
    zeta=hyparams.zeta;
    j=hyparams.j;
    params=Param(zeta,j);
else
    a=hyparams.a;
    b=hyparams.b;
    c=hyparams.c;
    d=hyparams.d;
    tau=hyparams.tau;
    params=Param_DL(a,b,c,d,tau);
end
T=hyparams.num_per*2*pi/params.omega;
num_steps=hyparams.num_steps*hyparams.num_per;
to_cut=T*hyparams.to_cut;
maxstep=T/num_steps;
if hyparams.eq=="MT"
    sol=Sol(@(t)hyparams.ic*[sin(params.omega*t);params.omega*cos(params.omega*t)], ...
        T,params,'timestep',maxstep,'int_reltol',hyparams.reltol).CutTransient(to_cut);
else
    a=params.a;
    b=params.b;
    c=params.c;
    d=params.d;
    omega=params.omega;
    ics=[(-1/2).*c.*d.^(-1)+(-1/2).*(d.^(-2).*(c.^2+4.*d)).^(1/2),(-1).* ...
  a.*((-1/2).*c.*d.^(-1)+(-1/2).*(d.^(-2).*(c.^2+4.*d)).^(1/2)).^2+( ...
  -1).*b.*((-1/2).*c.*d.^(-1)+(-1/2).*(d.^(-2).*(c.^2+4.*d)).^(1/2)) ...
  .^3+((-1/2).*c.*d.^(-1)+(-1/2).*(d.^(-2).*(c.^2+4.*d)).^(1/2)).*( ...
  omega.^(-2).*(1+(-1).*omega.^4)).^(1/2);(-1/2).*c.*d.^(-1)+(1/2).* ...
  (d.^(-2).*(c.^2+4.*d)).^(1/2),(-1).*a.*((-1/2).*c.*d.^(-1)+(1/2).* ...
  (d.^(-2).*(c.^2+4.*d)).^(1/2)).^2+(-1).*b.*((-1/2).*c.*d.^(-1)+( ...
  1/2).*(d.^(-2).*(c.^2+4.*d)).^(1/2)).^3+((-1/2).*c.*d.^(-1)+(1/2) ...
  .*(d.^(-2).*(c.^2+4.*d)).^(1/2)).*(omega.^(-2).*(1+(-1).*omega.^4) ...
  ).^(1/2)];
    sol=DelayedLienard(@(t)hyparams.ic*[sin(params.omega*t);params.omega*cos(params.omega*t)], ...
        T,params,'timestep',maxstep,'int_reltol',hyparams.reltol).CutTransient(to_cut);
end
if toplot
    figure(Name="Phase-plane plot of the training time series projected" + ...
        " to the center manifold")
    hold on
    scatter(sol.eta(1,1),sol.eta(2,1),"filled",'MarkerFaceColor','red')
    plot(sol.eta(1,:),sol.eta(2,:),'Color','blue')
end
mdl=Model(params,{sol},{},FitParams(3));
rmse_train=mdl.rmse;
if hyparams.eq=="MT"
    anal_coefs=AnalyticCoefs(params);
else
    anal_coefs=zeros(7,2);
end
if toplot
Plotter.BarPlotCoefs(anal_coefs,mdl,1);
Plotter.BarPlotCoefs(anal_coefs,mdl,2);
[fig1,fig2]=Plotter.PlotResiduals(mdl,sol,length(sol.t)/10);
set(fig1,"Name","Residuals for equation 1")
set(fig2,"Name","Residuals for equation 2")
fixes={[],[2,4,6,7],[1,3,5]};
names=["nonlinear", "quadratic","cubic"];
for i=1:2
    for fix_id=1:3
        fix=fixes{fix_id};
        fig=Plotter.CompareProjection(params,sol,i,"fixcoefs",fix);
        set(fig,"Name","Comparison of the projection for equation "+string(i)+...
            " of the "+names(fix_id)+" terms");
    end
end
end
anal_coefs=[params.omega,-params.omega;AnalyticCoefs(params)];
anal_pl=ComputePL_coefs(anal_coefs);
plerr=abs(mdl.PL-anal_pl);
coefErr=norm(anal_coefs(2:8,:)-mdl.coefs(2:8,:));
if ~isempty(monitor)
    monitor.Metrics=["RMSETraining","PLError","PLErrorPercent","CoefError",...
        "AnalyticPL","ApproximatedPL"];
    recordMetrics(monitor,1,RMSETraining=rmse_train, ...
        PLError=plerr,PLErrorPercent=abs(plerr/anal_pl)*100, CoefError=coefErr, ...
        AnalyticPL=anal_pl,ApproximatedPL=mdl.PL);
end
end