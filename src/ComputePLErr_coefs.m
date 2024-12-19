function [err] = ComputePLErr_coefs(coefs,coef_errs)
    inv_omega=1/(8*coefs(1,1));
    err=zeros(7,2);
    err(2,2)=3/8;
    err(7,1)=3/8;
    err(4,1)=1/8;
    err(6,2)=1/8;
    err([1,5],1)=inv_omega*(-coefs(4,1)+coefs(6,2)-coefs(2,2));
    err(1,1)=err(1,1)-inv_omega*(coefs(6,2)+coefs(2,2));
    err(3,1)=-inv_omega*(err(1,1)+err(5,1));
    err(5,1)=err(5,1)+inv_omega*(coefs(6,2)+coefs(2,2));
    err([1,5],2)=inv_omega*(coefs(4,2)+coefs(6,1)-coefs(2,1));
    err(1,2)=err(1,2)-inv_omega*(coefs(6,1)+coefs(2,1));
    err(3,2)=inv_omega*(coefs(6,2)+coefs(2,2));
    err(5,2)=err(5,2)+inv_omega*(coefs(6,1)+coefs(2,1));
    err=err.*coef_errs;
end

