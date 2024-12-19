function [coefs] = ComputeCoefPL(coefs,pl,varargin)
% varargin: row, col if coefs contains no NaNs.
if ~any(isnan(coefs),"all")
    coefs(varargin{1},varargin{2})=NaN;
end
[row,col]=find(isnan(coefs));
coef_sym=sym("coef",[8,2]);
PL=1/(8*coef_sym(1,1))*((coef_sym(6,1)+coef_sym(2,1))*(-coef_sym(4,1)+coef_sym(6,2)-coef_sym(2,2))+ ...
    (coef_sym(6,2)+coef_sym(2,2))*(coef_sym(4,2)+coef_sym(6,1)-coef_sym(2,1)))+1/8*(3*coef_sym(8,1)+ ...
    coef_sym(5,1)+coef_sym(7,2)+3*coef_sym(3,2));
for i=1:8
    for j=1:2
        if i~=row || j~=col
            PL=subs(PL,coef_sym(i,j),coefs(i,j));
        end
    end
end
coef=vpasolve(PL==pl);
coefs(row,col)=coef;
end

