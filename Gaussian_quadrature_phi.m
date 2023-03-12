function Phi_qc = Gaussian_quadrature_phi(msh,Gp)


%Number of elements 
ne = size(msh.simp,1);


for q=1:length(Gp.weights)
    
    %Form and vectorise p
    p = Gp.coord(:,:,q);
    pt = p.'; 
    pv = pt(:); %3ne x 1 column vector

    %Left multiply D by pv
    P = sparse(1:3*ne,1:3*ne,pv,3*ne,3*ne,3*ne);

    %Scale the rows of D with pv
    S = P*msh.phi_grad;

    %Form the martix Phi: ne x nn
    Phi_q = S(1:3:end,:) + S(2:3:end,:) + S(3:3:end,:);

    %to this we should add the constant term of each basis function so that
    %evaluated on its node has value of 1

    Phi_qc{q} = Phi_q + msh.phi_cons;


end %for q integration points
