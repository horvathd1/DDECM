classdef Param_DL
    properties
        a; b; c; d; tau; k; omega; gamma
    end

    methods

        function obj = Param_DL(a, b, c, d, tau)
            %PARAM_DL Construct an instance of this class
            %   Detailed explanation goes here
            obj.a = a;
            obj.b = b;
            obj.c = c;
            obj.d = d;
            obj.tau = tau;
            obj.k = ComputeCritParam(tau);
            [obj.omega, obj.gamma] = obj.ComputeLS(obj.k, obj.tau);
        end

        function res = eq(obj, obj2)
            res = obj.a == obj2.a && obj.b == obj2.b && obj.c == obj2.c && ...
                obj.d == obj2.d && obj.tau == obj2.tau && obj.k == obj2.k;
        end

    end

    methods (Static)

        function [omega, gamma] = ComputeLS(k, tau)
            syms om
            omega = double(solve(k ^ 2 * om ^ 2 == 1 - om ^ 4, om));
            omega = omega(omega >= 0 & omega <= 1&imag(omega)==0);
            gamma = (omega ^ 2 * (2 + k * tau)) / ((k - tau * omega ^ 2) ^ 2 + (2 * omega + k * tau * omega) ^ 2);
        end

    end

end
