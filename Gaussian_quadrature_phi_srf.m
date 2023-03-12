function Phi_qc = Gaussian_quadrature_phi_srf(msh,Gp,Gamma1)


    % Number of srface triangles 
    ne = size(msh.srf,1);
    
    notGamma1 = setdiff(1:length(msh.srf),Gamma1);
    
    for q=1:length(Gp.weights)
        
        %Form and vectorise p
        p = Gp.coord(:,:,q);
        pt = p.'; 
        pv = pt(:); %3ne x 1 column vector
    
        %Left multiply D by pv
        P = sparse(1:3*ne,1:3*ne,pv,3*ne,3*ne,3*ne);
    
        %Scale the rows of D with pv
        S = P*msh.Bgrad;
    
        %Form the martix Phi: ne x nn
        Phi_q = S(1:3:end,:) + S(2:3:end,:) + S(3:3:end,:);
    
        %to this we should add the constant term of each basis function so that
        %evaluated on its node has value of 1
    
        Phi_qc{q} = Phi_q + msh.bcons;
        
        for i=1:length(notGamma1)
            Phi_qc{q}(notGamma1(i),:) = 0;
        end
    
    
    end %for q integration points

end








