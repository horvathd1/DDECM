classdef Sol
    properties
        t;
        y;
        eta;
        nonlin;
        params;
        int_reltol;
        R;
        d1;
        d2;
    end
    
    methods
        function obj = Sol(x0,T,params,varargin)
            p=inputParser();
            p.addParameter("timestep",0.1);
            p.addParameter("int_reltol",6);
            p.parse(varargin{:});
            obj.int_reltol=p.Results.int_reltol;
            obj.params=params;
            sol=obj.Equation(x0,T,p.Results.timestep);
            obj=obj.AddLS();
            obj.t=sol.x;
            obj.y=sol.y;
            [y1,y2,nl1,nl2]=obj.ProjectionCM(x0);
            obj.t=obj.t(2:end-4);
            obj.y=obj.y(:,2:end-4);
            obj.eta=[arrayfun(y1,obj.t);arrayfun(y2,obj.t)];
            obj.nonlin=[arrayfun(nl1,obj.t);arrayfun(nl2,obj.t)];
        end
        function obj=AddLS(obj)
            gamma=obj.params.gamma;
            omega=obj.params.omega;
            zeta=obj.params.zeta;
            obj.R=[0,0;obj.params.p,0];
            obj.d1=[2*gamma*(2*zeta^2+2*zeta+1);2*gamma*zeta];
            obj.d2=[2*gamma*omega*zeta;2*gamma*omega];
        end
        function [y1,y2,nonlin1,nonlin2]=ProjectionCM(obj,ic)
            omega=obj.params.omega;
            tau=obj.params.tau;
            n1=@(theta)cos(omega*theta).*obj.d1-sin(omega*theta).*obj.d2;
            n2=@(theta)sin(omega*theta).*obj.d1+cos(omega*theta).*obj.d2;
            % Calculate z_t
            z=@(t)interp1(obj.t',obj.y',t,'linear','extrap')';
            dz=@(t)interp1(obj.t',(gradient(obj.y)./gradient(obj.t))',t, ...
                'linear','extrap')';
            if isnumeric(ic)
                zt=@(t,phi) z(t+phi).*(t+phi>=0)+ic.*(t+phi<0);
                dzt=@(t,phi) dz(t+phi).*(t+phi>=0)+0.*(t+phi<0);
            else
                syms f(x)
                f(x)=ic(x);
                df=diff(f,x);
                dic=matlabFunction(df);
                zt=@(t,phi) z(t+phi).*(t+phi>=0)+ic(t+phi).*(t+phi<0);
                dzt=@(t,phi) dz(t+phi).*(t+phi>=0)+dic(t+phi).*(t+phi<0);
            end
            % Calculate y1=(n1,z_t), y2=(n2,z_t)
            integrand=@(t,d,fun)diag(n1(d+obj.params.tau)'*obj.R*fun(t,d))';
            dot_bilin=@(t,fun)dot(n1(0),fun(t,0));
            y1=@(t)dot_bilin(t,zt)+integral(@(d)integrand(t,d,zt),-tau,0, ...
                "RelTol",obj.int_reltol);
            integrand2=@(t,d,fun)diag(n2(d+obj.params.tau)'*obj.R*fun(t,d))';
            dot_bilin2=@(t,fun)dot(n2(0),fun(t,0));
            y2=@(t)dot_bilin2(t,zt)+integral(@(d)integrand2(t,d,zt),-tau,0, ...
                "RelTol",obj.int_reltol);
            nonlin1=@(t)dot_bilin(t,dzt)-omega*dot_bilin2(t,zt)+...
                integral(@(d)integrand(t,d,dzt)-omega*integrand2(t,d,zt),-tau,0, ...
                "RelTol",obj.int_reltol);
            nonlin2=@(t)dot_bilin2(t,dzt)+omega*dot_bilin(t,zt)+...
                integral(@(d)integrand2(t,d,dzt)+omega*integrand(t,d,zt),-tau,0, ...
                "RelTol",obj.int_reltol);
        end
        function [sol]=Equation(obj,x0,T,maxstep)
            rhs=@(t,y,Z)[y(2);
                -2*obj.params.zeta*y(2)-(1+obj.params.p)*y(1)+...
                obj.params.p*Z(1,1)+3*obj.params.p/10*((y(1)-Z(1,1)).^2-...
                (y(1)-Z(1,1)).^3)];
            sol=dde23(rhs,obj.params.tau, x0,[0 T],ddeset('MaxStep',maxstep));
        end
        function [obj]=AppendSol(obj,sol)
            if obj.params==sol.params
                sol.t=sol.t+obj.t(end)-sol.t(1)+obj.t(end)-obj.t(end-1);
                obj.t=[obj.t sol.t];
                obj.y=[obj.y sol.y];
                obj.eta=[obj.eta sol.eta];
                obj.nonlin=[obj.nonlin sol.nonlin];
            else
                throw(MException('Sol:ParamNotEqual', ...
                    'To append two sol objects, the parameters must match'))
            end
        end
        function [obj]=CutTransient(obj,time)
            obj.y=obj.y(:,obj.t>time);
            obj.eta=obj.eta(:,obj.t>time);
            obj.nonlin=obj.nonlin(:,obj.t>time);
            obj.t=obj.t(:,obj.t>time)-time;
        end
    end
end

