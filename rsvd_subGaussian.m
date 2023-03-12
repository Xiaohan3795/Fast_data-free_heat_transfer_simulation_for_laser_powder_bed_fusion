function [U] = rsvd_subGaussian(A,r,k2)
%
%
    omega1 = randn(size(A,2),r);
    B = A * omega1;
    [Q,~] = qr(B,0);
    omega2 = randn(k2,size(A,1));
    C = (omega2*Q)\(omega2*A);
    [U,~,~] = svd(C);
    U = Q*U;
end

