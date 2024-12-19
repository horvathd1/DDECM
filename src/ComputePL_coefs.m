function PL=ComputePL_coefs(coefs)
%   ComputePL_coefs: Compute the Poincaré-Lyapunov constant from coefficients
    coef1=coefs(:,1);
    coef2=coefs(:,2);
    PL=1/(8*coef1(1))*((coef1(6)+coef1(2))*(-coef1(4)+coef2(6)-coef2(2))+ ...
    (coef2(6)+coef2(2))*(coef2(4)+coef1(6)-coef1(2)))+1/8*(3*coef1(8)+ ...
    coef1(5)+coef2(7)+3*coef2(3));
end