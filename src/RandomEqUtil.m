classdef RandomEqUtil
    properties
        
    end
    
    methods(Static)
        function terms=nonlinterms(y1,y2)
            terms=[y2.^2;y2.^3;y1.*y2;y1.*y2.^2;y1.^2;y1.^2.*y2;y1.^3];
        end
        function terms=nonlinterms_str(index)
            terms=char(["y"+string(3-index),"y2^2","y2^3","y1y2","y1y2^2","y1^2","y1^2y2","y1^3"]');
        end
        function rmse=ComputeRMSE(sol,fit_coefs,index)
            if string(index)=="all"
                rmse=RandomEqUtil.ComputeRMSE(sol,fit_coefs,1)+...
                    RandomEqUtil.ComputeRMSE(sol,fit_coefs,2);
                return
            end
            term1=fit_coefs(1,1)*[sol.y(2,:);-sol.y(1,:)];
            term1=term1(index,:);
            rmse=norm((term1+ ...
                (RandomEqUtil.nonlinterms(sol.y(1,:),sol.y(2,:))'* ...
                fit_coefs(2:end,index))')- ...
                gradient(sol.y(index,:))./gradient(sol.x));
        end
        function res=Predict(coefs,y1,y2)
            res=RandomEqUtil.nonlinterms(y1,y2)'*coefs;
        end
        function fig=PlotNonlinComparison(sol,omega,coefs,index)
            preds=RandomEqUtil.Predict(coefs,sol.y(1,:),sol.y(2,:))';
            term1=omega*[sol.y(2,:);-sol.y(1,:)];
            term1=term1(index,:);
            dy=gradient(sol.y(index,:))./gradient(sol.x);
            fig=figure;
            hold on
            plot(sol.x(2:end-1),dy(2:end-1)-term1(2:end-1));
            plot(sol.x,preds(index,:));
        end
        function [x_eq,y_eq]=ComputeEquilibria(coefs)
            syms x y;
            coefs_a=coefs(:,1);
            coefs_b=coefs(:,2);
            terms=[y^2;y^3;x*y;x*y^2;x^2;x^2*y;x^3];
            [x_eq,y_eq]=vpasolve([coefs_a(1)*y;coefs_b(1)*x]+ ...
                [coefs_a(2:end)';coefs_b(2:end)']*terms,[x;y]);
        end
        function tab=TableCoefsOmega(coefs,fit_coefs,omega,fit_omega,index,digits)
                tab=RandomEqUtil.TableCoefs(RandomEqUtil.ExtendCoefs(coefs,omega), ...
                    RandomEqUtil.ExtendCoefs(fit_coefs,fit_omega),index,digits);
        end
        function tab=TableCoefs(coefs,fit_coefs,index,digits)
            vars=char(["y1","y2^2","y2^3","y1y2","y1y2^2","y1^2","y1^2y2","y1^3"]');
            orig_index=(1:8)';
            tab=table(orig_index,vars,round(coefs(:,index),digits),round(fit_coefs(:,index),digits), ...
                VariableNames=["Original Index","Variable","Original","Computed"]);
            tab=tab([1,6,4,2,8,7,5,3],:);
        end
        function res=MaxErrAtSmallCoef(coefs,fit_coefs, ...
                small_coef_threshold,max_err_threshold)
            [maxerr,ids]=max(abs((coefs-fit_coefs)./fit_coefs*100),[],"all");
            res=maxerr>max_err_threshold&...
                coefs(ids)/max(abs(coefs),[],"all")<small_coef_threshold;
        end
        function fig=PlotVectorField(coefs,xval,yval,varargin)
            [x,y]=meshgrid(xval,yval);
            nonlin=reshape(RandomEqUtil.nonlinterms(x,y)',size(x,1),size(x,2),7);
            u=coefs(1,1)*y+squeeze(tensorprod(coefs(2:end,1),nonlin,1,3));
            v=coefs(1,2)*x+squeeze(tensorprod(coefs(2:end,2),nonlin,1,3));
            [equil_x,equil_y]=RandomEqUtil.ComputeEquilibria(coefs);
            equil_x=equil_x(real(equil_x)==equil_x);
            equil_y=equil_y(real(equil_y)==equil_y);
            fig=figure;
            hold on
            streamline(x,y,u,v,x,y,varargin{:});
            scatter(equil_x,equil_y);
            xlim([min(xval),max(xval)])
            ylim([min(yval),max(yval)])
        end
        function coefs=ExtendCoefs(coefs,omega)
            coefs=[omega,-omega;coefs];
        end
        function pl_coefs=getPLCoefs(coefs)
            % Coefs is Nx7x2
            pl_coefs=zeros(size(coefs,1),10);
            keep=logical([1,0,1,1,1,0,1,1,1,1,0,1,1,0]');
            for i=1:size(coefs,1)
                coef=reshape(coefs(i,:,:),[14,1]);
                pl_coefs(i,:)=coef(keep);
            end
        end
        function pl_coefs=getNonPLCoefs(coefs)
            % Coefs is Nx7x2
            pl_coefs=zeros(size(coefs,1),4);
            keep=~logical([1,0,1,1,1,0,1,1,1,1,0,1,1,0]');
            for i=1:size(coefs,1)
                coef=reshape(coefs(i,:,:),[14,1]);
                pl_coefs(i,:)=coef(keep);
            end
        end
        function tab=TableCoefSummary(values,digits)
            % Values is Nx7x2
            terms=[RandomEqUtil.nonlinterms_str(1);RandomEqUtil.nonlinterms_str(2)];
            tab=table((1:7)',terms(2:8,:),reshape(mean(values,1),[7,2]),reshape(max(values,[],1),[7,2]), ...
                reshape(std(values,1),[7,2]),'VariableNames',["Index","Term","Mean","Max","Std"]);
            tab(:,3:end).Variables=round(tab(:,3:end).Variables,digits);
            tab=tab([5,3,1,7,6,4,2],:);
        end
        function [coefs,dists,omegas,relerrs,abserrs]=...
                    ComputeHistoryVariables(params,models)
            coefs=cellfun(@(x)[x.coefs_a(2:end),x.coefs_b(2:end)],params,UniformOutput=false);
            dists=cellfun(@(x)norm(x.ic0)*5,params);
            omegas=cellfun(@(x)x.coefs_a(1),params);
            relerrs=zeros(length(coefs),7,2);
            for i=1:length(coefs)
                relerrs(i,:,:)=abs((coefs{i}-models{i})./coefs{i}*100);
            end
            abserrs=zeros(length(coefs),7,2);
            for i=1:length(coefs)
                abserrs(i,:,:)=abs((coefs{i}-models{i}));
            end
        end
        function [ic]=ComputeIC(coefs_a,coefs_b)
            x=sym("x");
            y=sym("y");
            terms=[y^2;y^3;x*y;x*y^2;x^2;x^2*y;x^3];
            [xsol,ysol]=vpasolve([coefs_a(1)*y;coefs_b(1)*x]+ ...
                [coefs_a(2:end)';coefs_b(2:end)']*terms,[x;y]);
            xsol=xsol(real(xsol)==xsol);
            ysol=ysol(real(ysol)==ysol);
            if length(xsol)==1
                ic=[0.01,0];
            else
                ic=double(min(sqrt(xsol(xsol~=0).^2+ysol(ysol~=0).^2)*[1,0]))/5;
                ic=[min(ic(1),0.01),0];
            end
        end
        function coefs=CoefFromModels(models,omega)
            coefs=RandomEqUtil.ExtendCoefs(cell2mat( ...
                cellfun(@(x)x.Coefficients.Estimate,models, ...
                UniformOutput=false)),omega);
        end
        function dist=DistanceFromStart(sol)
            dist=(sol.y(1,end)-sol.y(1,1)).^2+(sol.y(2,end)-sol.y(2,1)).^2;
        end
    end
end

