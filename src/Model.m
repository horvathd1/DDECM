classdef Model
    properties
        sols_train;
        sols_test;
        preds_test;
        err_preds;
        coefs;
        powers;
        fitParams;
        PL;
        params;
        rmse;
    end
    
    methods
        function obj = Model(params,sols_train,sols_test,fitparams)
            obj.params=params;
            obj.sols_train=sols_train;
            obj.sols_test=sols_test;
            obj.fitParams=fitparams;
            obj=obj.fitCMDE();
            obj.preds_test=cell(1,length(obj.sols_test));
            for i=1:length(obj.preds_test)
                obj.preds_test{i}=obj.CalculateTrajectories(i);
                obj.err_preds{i}=norm(obj.preds_test{i}.y-obj.sols_test{i}.eta,'fro');
                obj.err_preds{i}=obj.err_preds{i}(1);
            end
            obj=obj.ComputePL();
        end
        function [obj]=fitCMDE(obj)
            sol_train=obj.sols_train{1};
            for i=2:length(obj.sols_train)
                sol_train=sol_train.Append(obj.sols_train{i});
            end
            eta=sol_train.eta;
            nonlin=sol_train.nonlin;
            [monoms,obj.powers]=obj.calculateMonomials(eta(1,:),eta(2,:), ...
                obj.fitParams.n);
            X=monoms(obj.fitParams.free,:);
            scale=max(abs(X),[],2);
            X=X./scale;
            anal_coefs=AnalyticCoefs(obj.params);
            fix_coefs=anal_coefs(obj.fitParams.fix,1)';
            to_sub=fix_coefs*monoms(obj.fitParams.fix,:);
            Y=(nonlin(1,:)-to_sub)';
            scale2=max(abs(Y));
            Y=Y./scale2;
            coef=lsqlin(X',Y,[],[]);
            coef1=[obj.params.omega;obj.insertFreeCoefs(fix_coefs,scale2.*coef./scale)];
            fix_coefs=anal_coefs(obj.fitParams.fix,2)';
            to_sub=fix_coefs*monoms(obj.fitParams.fix,:);
            Y=(nonlin(2,:)-to_sub)';
            scale2=max(abs(Y));
            Y=Y./scale2;
            coef=lsqlin(X',Y,[],[]);
            coef2=[-obj.params.omega;obj.insertFreeCoefs(fix_coefs,scale2.*coef./scale)];
            obj.coefs=[coef1 coef2];
            obj.rmse=obj.ComputeRMSE(sol_train);
        end

        function coefs=insertFreeCoefs(obj,fix_coefs,free_coefs)
            coefs=fix_coefs;
            j=1;
            for i=obj.fitParams.free
                coefs=[coefs(1:i-1),free_coefs(j),coefs(i:end)];
                j=j+1;
            end
            coefs=coefs';
        end
        
        function [rmse]=ComputeRMSE(obj,sol)
        %COMPUTERMSE: Computes the RMSE given a Sol object
        monoms=obj.calculateMonomials(sol.eta(1,:),sol.eta(2,:), ...
            obj.fitParams.n);
        X={monoms',monoms'};
        preds={X{1}*obj.coefs(2:end,1),X{2}*obj.coefs(2:end,2)};
        rmse=sqrt(1/length(sol.nonlin)*(norm(preds{1}'-sol.nonlin(1,:))^2+...
                   norm(preds{2}'-sol.nonlin(2,:))^2));
        end

        function [str]=ToString(obj)
            str=obj.cmdeToString(obj.coefs(:,1),obj.powers,'y1');
            str=strcat(str,obj.cmdeToString(obj.coefs(:,2),obj.powers,'y2'));
        end
        function obj=ComputePL(obj)
            obj.PL=ComputePL_coefs(obj.coefs);
        end
        function [sol]=CalculateTrajectories(obj,index)
            eta=obj.sols_test{index}.eta;
            ic=eta(:,1);
            coef1=obj.coefs(:,1);
            coef2=obj.coefs(:,2);
            T=obj.sols_test{index}.t(end);
            maxstep=max(diff(obj.sols_test{index}.t));
            [t,y]=ode45(@(t,y) [coef1'*[y(2);obj.calculateMonomials(y(1),y(2),obj.fitParams.n)]; ...
                coef2'*[y(1);obj.calculateMonomials(y(1),y(2),obj.fitParams.n)]], ...
                [0,T],ic,odeset('MaxStep',maxstep));
            sol.t=obj.sols_test{index}.t;
            sol.y=interp1(t,y,sol.t)';
        end
        function [deta]=Predict(obj,eta)
            deta(1,:)=obj.coefs(:,1)'*...
                [eta(2,:);obj.calculateMonomials(eta(1,:),eta(2,:),obj.fitParams.n)];
            deta(2,:)=obj.coefs(:,2)'*...
                [eta(2,:);obj.calculateMonomials(eta(1,:),eta(2,:),obj.fitParams.n)];
        end
    end
    methods(Static)
        function [monoms,powers]=calculateMonomials(y1,y2,n)
            monoms=zeros((n+1)^2,length(y1));
            powers=zeros((n+1)^2,2);
            index=1;
            for i=0:3
                for j=0:3
                    if (i+j<2) || (i+j>3)
                        continue
                    end
                    powers(index,:)=[i,j];
                    monoms(index,:)=y1.^i.*y2.^j;
                    index=index+1;
                end
            end
            for i=0:n
                for j=0:n
                    if (i+j<4) || (i+j>n)
                        continue
                    end
                    powers(index,:)=[i,j];
                    monoms(index,:)=y1.^i.*y2.^j;
                    index=index+1;
                end
            end
            powers(index:end,:)=[];
            monoms(index:end,:)=[];
        end
        function [de]=cmdeToString(coef,pow,var)
            if var=="y1"
                var2="y2";
            else
                var2="y1";
            end
            de=strcat(var,"'=",string(coef(1)),"*",var2);
            for i=1:length(coef)-1
                if coef(i+1)>0
                    de=strcat(de,"+",string(coef(i+1)),"*",Model.powToStr(pow(i,1),pow(i,2)));
                elseif coef(i+1)<0
                    de=strcat(de,string(coef(i+1)),"*",Model.powToStr(pow(i,1),pow(i,2)));
                end
            end
        end
        
        function [str]=powToStr(p1,p2)
            if p1==0
                str=strcat("y2^",string(p2));
            elseif p2==0
                str=strcat("y1^",string(p1));
            else
                str=strcat("y1^",string(p1),"*y2^",string(p2));
            end
        end
    end
end

