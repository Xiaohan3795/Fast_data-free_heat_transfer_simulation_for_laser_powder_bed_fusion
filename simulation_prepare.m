function [bpv, bpv_dif, C, a, Phim,Phir,barPhir,barPhik,Phik] = simulation_prepare(msh,h,u_a,sigma,varepsilon,ns,bpt)


        % the nodes on the bottom surface
        bp = find(msh.csrf(:,3)==0);
        bpv = unique(msh.srf(bp,:));

        % the nodes not on the bottom surface
        bpv_dif = setdiff(1:length(msh.vtx),bpv);

        % surfaces not on the bottom surface
        Gam1ind = find(msh.csrf(:,3)>0);

        % the matrix of convection heat loss: C
        % the vector of ambient temperature in heat loss: a
        [C,a] = make_surface_integrals_Gaussian(msh,h,u_a,sigma,varepsilon,Gam1ind);
  
        % Phim
        Phi_qc = Gaussian_quadrature_phi(msh,msh.Gq{ns});
        Phim = [];
        for q = 1:length(msh.Gq{ns}.weights)
            Phim = [Phim; Phi_qc{q}(:,bpv_dif)];
        end 
        % Phir
        Phi_qc_srf_5 = Gaussian_quadrature_phi_srf(msh,msh.Gqsrf{5},Gam1ind);
        Phir = [];
        for q = 1:length(msh.Gqsrf{5}.weights)
            Phir = [Phir; Phi_qc_srf_5{q}(:,bpv_dif)];
        end 
        % barWr
        barPhir = [];
        for q = 1:length(msh.Gqsrf{5}.weights)
            barPhir = [barPhir; Phi_qc_srf_5{q}(:,bpv)*(ones(length(bpv),1)*bpt)];
        end
        % barPhik
        barPhik = msh.phi_grad(:,bpv)*(ones(length(bpv),1)*bpt);
        % Phik
        Phik = msh.phi_grad(:,bpv_dif);

end