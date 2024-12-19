function [PLs]=AnalyticPL_DL(a,b,c,d,tau)
    PLs=[a;b;c;d;tau;zeros(1,length(a))];
    for i=1:length(a)
        PLs(6,i)=ComputePL(PLs(1,i),PLs(2,i),PLs(3,i),PLs(4,i),PLs(5,i));
    end
end

function [PL]=ComputePL(a,b,c,d,tau)
    param=Param_DL(a,b,c,d,tau);
    zeta=2*param.omega/(5+12*param.omega^4-8*param.omega^6);
    Omega=4*param.gamma/(param.omega^2*(2+param.k*tau));
    d12=Omega*1/2*(param.k-param.omega^2*tau);
    d22=Omega*param.omega*(1+1/2*param.k*tau);
    l1=3/8*d*param.omega^2+a^2/4*param.omega*zeta*(1+4*param.omega^2-2*param.omega^4)-...
    a*c/2*param.k*param.omega*zeta*(1+param.omega^2+param.omega^4)+c^2/4*param.omega*zeta*...
    (11/2+param.k^2+2*param.omega^2+12*param.omega^4-12*param.omega^6);
    l2=3/8*b*param.omega-3/8*d*param.k*param.omega+a^2/2*param.k*param.omega^2*zeta*(1-param.omega^2)+...
        a*c/4*zeta*(7/2+param.omega^2+10*param.omega^4-10*param.omega^6)+c^2/4*param.k*zeta*...
        (-11/2+param.omega^2-12*param.omega^4+12*param.omega^6);
    PL=l1*d12+l2*d22;
end