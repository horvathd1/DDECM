classdef Plotter
    properties
        
    end
    
    methods
        function obj = Plotter()
        end
    end
    methods(Static)
        function PlotDerivatives(model,sol)
            eta=sol.eta;
            deta=sol.deta;
            x=sol.t;
            coef1=model.coefs(:,1);
            coef2=model.coefs(:,2);
            figure;
            hold on
            plot(x,deta(1,:))
            plot(x,coef1'*[eta(2,:);calculateMonomials(eta(1,:),eta(2,:),n)])
            figure;
            hold on
            plot(x,deta(2,:))
            plot(x,coef2'*[eta(1,:);calculateMonomials(eta(1,:),eta(2,:),n)])
        end
        function PlotTrajectories(model,sol)
            eta=sol.eta;
            x=sol.t;
            coef1=model.coefs(:,1);
            coef2=model.coefs(:,2);
            ic=eta(:,1);
            [t,y]=ode45(@(t,y) [coef1'*[y(2);calculateMonomials(y(1),y(2),n)]; ...
                coef2'*[y(1);calculateMonomials(y(1),y(2),n)]],[0,1000],ic);
            figure;
            disp(CustomError(t,y',x,eta))
            plot(x,eta(1,:));
            hold on
            plot(t,y(:,1),"r")
            legend({'Test set', 'Prediction'})
            xlabel('t')
            ylabel('$\xi_1$','Interpreter','latex')
            figure;
            plot(x,eta(2,:));
            hold on
            plot(t,y(:,2),"r")
            legend({'Test set', 'Prediction'})
            xlabel('t')
            ylabel('$\xi_1$','Interpreter','latex')
        end
        function [fig1,fig2]= PlotResiduals(model,sol,n_avg)
            eta=sol.eta;
            deta=sol.nonlin;
            coef1=model.coefs(:,1);
            coef2=model.coefs(:,2);
            pred=[coef1(2:end,:)'*[model.calculateMonomials(eta(1,:),eta(2,:),model.fitParams.n)];...
                coef2(2:end,:)'*[model.calculateMonomials(eta(1,:),eta(2,:),model.fitParams.n)]];
            resid=deta-pred;
            [~,idx1]=sort(pred(1,:));
            [~,idx2]=sort(pred(2,:));
            beta=1-1/n_avg;
            resid_avg1=Plotter.ExpWeightedAvg(resid(1,idx1),beta);
            resid_avg2=Plotter.ExpWeightedAvg(resid(2,idx2),beta);
            resid_std1=Plotter.ExpWeightedStd(resid(1,idx1),resid_avg1,beta);
            resid_std2=Plotter.ExpWeightedStd(resid(2,idx2),resid_avg2,beta);
            fig1=figure;
            hold on
            plot(pred(1,idx1),resid(1,idx1),'o')
            plot(pred(1,idx1),resid_avg1,'r-','LineWidth',2)
            plot(pred(1,idx1),resid_avg1+2*resid_std1,'r--', ...
                pred(1,idx1),resid_avg1-2*resid_std1,'r--','LineWidth',2)
            fig2=figure;
            hold on
            plot(pred(2,idx2),resid(2,idx2),'o')
            plot(pred(2,idx2),resid_avg2,'r-','LineWidth',2)
            plot(pred(2,idx2),resid_avg2+2*resid_std2,'r--', ...
                pred(2,idx2),resid_avg2-2*resid_std2,'r--','LineWidth',2)
        end
        function [res]= ExpWeightedAvg(data,beta)
            res=zeros(1,length(data));
            res(1)=data(1)*(1-beta);
            for i=2:length(data)
                res(i)=res(i-1)*beta+(1-beta)*data(i);
            end
            res=res./(1-beta.^(1:length(data)));
        end
        function [res]= ExpWeightedStd(data,data_avg,beta)
            stds=(data-data_avg).^2;
            res=sqrt(Plotter.ExpWeightedAvg(stds,beta));
        end
        function PlotAmplitudes(model,analytic,kmin)
            k=model.sols_train{1}.params.k;
            gamma=model.sols_train{1}.params.gamma;
            karr=k-kmin:0.00005:k;
            ampli=sqrt(gamma/model.PL*(karr-k));
            analytic_ampli=sqrt(gamma/analytic*(karr-k));
            figure;
            hold on
            plot(karr,analytic_ampli,'r-','LineWidth',2)
            plot(karr,ampli,'b--','LineWidth',2)
            legend(["Ground truth","Data-driven"])
            xlabel("p")
            ylabel("Amplitude")
            set(gca,'Fontsize',26)
        end
        function [fig]= PlotErrorComparison(error1,error2)
            error1=sort(error1,'descend');
            error2=sort(error2,'descend');
            fig=figure;
            hold on
            plot(error1);
            plot(error2);
        end
        function [fig]=BarCoefs(anal_coefs,approx_coefs,index)
            fig=figure(Name="Bar plot of the coefficients for equation "+ ...
                string(index));
            plot_coefs=~isnan(anal_coefs(:,index)');
            categories=["y_2^2","y_2^3","y_{1}y_2","y_{1}y_2^2","y_1^2","y_1^{2}y_2","y_1^3"];
            bar(reordercats(categorical(categories),[1,6,4,2,3,7,5]), ...
                [anal_coefs(plot_coefs,index)';approx_coefs(plot_coefs,index)']')
            legend(["Analytic","Approximated"])
            set(gca,'FontSize',24)
            ax=gca;
            ax.YAxis.Exponent=0;
        end
        function [fig]=BarPlotCoefs(anal_coefs,model,index)
            fig=Plotter.BarCoefs(anal_coefs,model.coefs(2:end,:),index);
        end
        function [anal_eta,anal_model]=ComputeAnalyticProjection(params,sol)
            anal_coefs=[params.omega,-params.omega;...
                AnalyticCoefs(params)];
            anal_model=Model(params,{sol},{sol}, ...
                FitParams(3));
            anal_model.coefs=anal_coefs;
            anal_eta=anal_model.CalculateTrajectories(1).y;
        end
        
        function res=SubtractCoefs(orig,coefs,eta,fix,index)
            size_coef=size(coefs);
            to_sub_coefs=zeros(size_coef(1)-1,1);
            to_sub_coefs(fix)=coefs(fix,index);
            monoms=Model.calculateMonomials(eta(1,:),eta(2,:),3);
            res=orig-to_sub_coefs'*monoms;
        end

        function fig=CompareProjection(params,sol,index,varargin)
            p=inputParser;
            p.addParameter("fixcoefs",[])
            p.parse(varargin{:});
            fig=figure;
            [anal_eta,anal_model]=Plotter.ComputeAnalyticProjection(params,sol);
            anal_deta=gradient(anal_eta)./gradient(sol.t);
            anal_nonlin=[anal_deta(1,:)-params.omega*anal_eta(2,:);
                anal_deta(2,:)+params.omega*anal_eta(1,:)];
            sol_nonlin=sol.nonlin(index,:);
            if ~isempty(p.Results.fixcoefs)
                anal_nonlin=Plotter.SubtractCoefs(anal_nonlin,anal_model.coefs, ...
                    anal_eta,p.Results.fixcoefs,index);
                sol_nonlin=Plotter.SubtractCoefs(sol_nonlin,anal_model.coefs, ...
                    anal_eta,p.Results.fixcoefs,index);
            end
            hold on
            plot(sol.t,sol_nonlin)
            plot(sol.t(2:end-1),anal_nonlin(index,2:end-1))
            legend(["Projected","Analytic coefficient"])
            xlabel("t")
            ylabel("eta")
        end
        function fig=ProjectionPhasePlane(params,sol,varargin)
            p=inputParser;
            p.addParameter("nonlin",false);
            p.addParameter("trim",100);
            p.parse(varargin{:});
            fig=figure;
            anal_eta=Plotter.ComputeAnalyticProjection(params,sol);
            hold on
            if ~p.Results.nonlin
                plot(sol.eta(1,:),sol.eta(2,:))
                plot(anal_eta(1,:),anal_eta(2,:))
            else
                plot(sol.nonlin(1,:),sol.nonlin(2,:))
                anal_deta=gradient(Plotter.ComputeAnalyticProjection(params,...
                    sol))./gradient(sol.t);
                trim=p.Results.trim;
                plot(anal_deta(1,trim:end-trim)-params.omega*anal_eta(2,trim:end-trim), ...
                    anal_deta(2,trim:end-trim)+params.omega*anal_eta(1,trim:end-trim))
            end
            legend(["Projected","Analytic coefficient"])
        end
        function fig=Histogram(values,maxval,varargin)
            % varargin: indices: indices of the dimensions other than the 
            % first in values
            % hist_params: parameters to give to the histogram function
            % xlabel: xlabel to use
            % mean(=NaN): whether to draw a line at the mean value
            % linespec(={}): Parameters to pass to xline
            p=inputParser;
            p.addParameter("indices",0);
            p.addParameter("hist_params",{});
            p.addParameter("xlabel","Relative error of the Poincaré-Lyapunov constant [%]")
            p.addParameter("mean",false);
            p.addParameter("linespec",{"Linewidth",2,"Color","black"});
            p.parse(varargin{:});
            fig=figure;
            hold on
            if class(p.Results.indices)=="cell"
                val=values(values(:,p.Results.indices{:})<=maxval, ...
                    p.Results.indices{:});
                m=mean(values(values(:,p.Results.indices{:}), ...
                    p.Results.indices{:}),"all");
                histogram(val,p.Results.hist_params{:})
            else
                val=values(abs(values)<=maxval);
                m=mean(values,"all");
                histogram(val,p.Results.hist_params{:})
            end
            if p.Results.mean
                xline(m,p.Results.linespec{:});
            end
            xlabel(p.Results.xlabel)
            ylabel("Number of runs")
            set(gca,'FontSize',24)
            if min(values)>0
                xlim([0,maxval])
            else
                xlim([-maxval,maxval])
            end
            ax=gca;
            ax.YAxis.Exponent=0;
            ax.XAxis.Exponent=0;
            ytickformat("%.0f")
            %xtickformat("%.0f")
        end
        function fig=PlotErrorSurface(errfun,xval,yval)
            [x_mesh,y_mesh]=meshgrid(xval,yval);
            z=x_mesh;
            for ind_x=1:length(xval)
                for ind_y=1:length(yval)
                    z(ind_y,ind_x)=errfun(xval(ind_x),yval(ind_y));
                end
            end
            fig=figure;
            surf(x_mesh,y_mesh,z)
        end
        function fig=PlotRMSESurface(sol,fit_coefs,index1,index2,xval,yval)
            function rmse=rmsefun(x,y)
                mod_fit_coefs=fit_coefs;
                mod_fit_coefs(index1(1),index1(2))=fit_coefs(index1(1),index1(2))+x;
                mod_fit_coefs(index2(1),index2(2))=fit_coefs(index2(1),index2(2))+y;
                rmse=RandomEqUtil.ComputeRMSE(sol,mod_fit_coefs,"all");
            end
            fig=Plotter.PlotErrorSurface(@rmsefun,xval,yval);
        end
        function fig=addIndexTooltip(fig,indices)
            fig.DataTipTemplate.DataTipRows(end+1)=[dataTipTextRow("index",indices)];
        end
        function tab=TableCoefs(model,index,digits)
            tab=RandomEqUtil.TableCoefsOmega(AnalyticCoefs(model.params), ...
                model.coefs(2:end,:),model.coefs(1,1),model.coefs(1,1), ...
                index,digits);
        end
    end
end

