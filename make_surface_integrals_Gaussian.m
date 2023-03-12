function [C,a] = make_surface_integrals_Gaussian(msh,h,u_a,sigma,varepsilon,Gam1ind)



    Phi_qc1 = Gaussian_quadrature_phi_srf(msh,msh.Gqsrf{1},Gam1ind);
    Phi_qc2 = Gaussian_quadrature_phi_srf(msh,msh.Gqsrf{2},Gam1ind);
    
    
    % Number of srface triangles 
    ns = size(msh.srf,1);

    %% Code for the vector surface integral Sa

    spdarea = sparse(1:ns,1:ns,msh.srfarea,ns,ns,ns);

    q=1;
    Sa = h*u_a*Phi_qc1{q}.'*msh.srfarea;
    Sa4 = sigma*varepsilon*u_a^4*Phi_qc1{q}.'*msh.srfarea;

    C = Phi_qc2{q}.'*h*msh.Gqsrf{2}.weights(q)*spdarea*Phi_qc2{q};
    for q=2:length(msh.Gqsrf{2}.weights)
        %evaluate S1 at the integration point 
        C = C + Phi_qc2{q}.'*h*msh.Gqsrf{2}.weights(q)*spdarea*Phi_qc2{q};
    end %for q integration points

    a = Sa+Sa4;

end

