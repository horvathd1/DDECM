classdef DelayedLienard < Sol
    properties

    end

    methods

        function obj = DelayedLienard(x0, T, params, varargin)
            %DELAYEDLIENARD Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@Sol(x0, T, params, varargin{:});
        end

        function [sol] = Equation(obj, x0, T, maxstep)
            rhs = @(t, y, Z)[y(2) - obj.params.k * y(1) + obj.params.a * y(1) ^ 2 + obj.params.b * y(1) ^ 3;
                             -Z(1, 1) + obj.params.c * Z(1, 1) ^ 2 + obj.params.d * Z(1, 1) ^ 3];
            sol = dde23(rhs, obj.params.tau, x0, [0 T], ddeset('MaxStep', maxstep));
        end

        function obj = AddLS(obj)
            omega = obj.params.omega;
            k=obj.params.k;
            tau=obj.params.tau;
            Omega=4*obj.params.gamma/(omega^2*(2+k*tau));
            obj.R = [0, 0; -1, 0];
            obj.d1 = Omega*[omega^2*(1+1/2*k*tau); 1/2*(k-omega^2*tau)];
            obj.d2 = Omega*[-omega/2*(k-omega^2*tau); omega*(1+1/2*k*tau)];
        end

    end

end
