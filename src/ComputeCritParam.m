function [k]=ComputeCritParam(tau)
    syms om
    om=vpasolve(1/om*atan(sqrt((1-om^4)/om^4))==tau,om,[0,1]);
    k=double(sqrt((1-om^4)/om^2));
end