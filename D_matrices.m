function [Dr,Dk,Dm] = D_matrices(nsk, nsm, msh, ut1_k, sigma, varepsilon)

    % Dr
    ux=[];
    for q=1:length(msh.Gqsrf{5}.weights)
        ux =[ux; dot(msh.Gqsrf{5}.Gqsrf_temp_compute{q}, ut1_k(msh.srf'))'];
    end
    neS4 = length(ux);
    Dr = sparse(1:neS4,1:neS4,sigma*varepsilon*ux.^3.*msh.Gqsrf{5}.wa);
    
    % Dk
    ux=[];
    for q=1:length(msh.Gq{nsk}.weights)
        ux =[ux; dot(msh.Gq{nsk}.Gq_temp_compute{q}, ut1_k(msh.simp'))'];
    end
    [kappa] = kappa_find(ux);
    kappa = kappa*1e-3;
    bb = kron((kappa.*msh.Gq{nsk}.wa)', [1,1,1]);
    Dk = reshape(bb, 3*size(msh.simp,1), length(msh.Gq{nsk}.weights));
    neK = size(Dk,1);
    Dk =sparse(1:neK,1:neK,sum(Dk,2));
    
    % Dm
    ux=[];
    for q=1:length(msh.Gq{nsm}.weights)
        ux =[ux; dot(msh.Gq{nsm}.Gq_temp_compute{q}, ut1_k(msh.simp'))'];
    end
    [rho,c] = rhoc_find(ux);
    rho = rho *1e-9;
    new = length(rho);
    Dm = sparse(1:new,1:new,rho.*c.*msh.Gq{nsm}.wa);
end

