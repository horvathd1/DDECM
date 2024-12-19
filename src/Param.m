classdef Param
    properties
        zeta;
        j;
        omega;
        gamma;
        p;
        tau;
    end
    
    methods
        function obj = Param(zeta,j)
            %PARAM_DL Construct an instance of this class
            %   Detailed explanation goes here
            obj.zeta=zeta;
            obj.j=j;
            [obj.omega,obj.gamma,obj.p,obj.tau]=obj.ComputeLS(obj.zeta,obj.j);
        end

        function res=eq(obj,obj2)
            res =obj.zeta==obj2.zeta && obj.j==obj2.j;
        end
    end
    methods(Static)
        function [omega,gamma,p,tau]=ComputeLS(zeta,j)
            omega=sqrt(1+2*zeta);
            p=2*zeta*(zeta+1);
            tau=2*(j*pi-atan(1/omega))/omega;
            gamma=1/(2*(1+zeta)^2*(1+zeta*tau));
        end

    end
end

