function [B,zeta_w] = sketch_index_noproj(W,Psi,sk,sp)
%
%
    n1 = size(W,1);
    n2 = size(Psi,2);
    
    U = rsvd_subGaussian_sketch(W,Psi,floor(n2*sp),1);
    U = U.^2;

    xi = sk/n2*sum(U,2);
    zeta = min([ones(n1,1),xi],[],2);
    B = find((rand(n1,1)-zeta)<0);
    zeta_w = 1./zeta(B);
end

