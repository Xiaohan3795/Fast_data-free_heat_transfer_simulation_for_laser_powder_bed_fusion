function [proj_Mt,proj_Mt_k_bar] = proj_sketch_M(msh,ut1_k,B_M,B_M_w,proj_Wm,nsm,dti,WM_k_bar)

                        ux=[];
                        for q=1:length(msh.Gq{nsm}.weights)
                            ux =[ux; dot(msh.Gq{nsm}.Gq_temp_compute{q}, ut1_k(msh.simp'))'];
                        end
                        [rho,c] = rhoc_find(ux(B_M));
                        rho = rho *1e-9;

                        new = length(rho);
                        DD = sparse(1:new,1:new,rho.*c.*msh.Gq{nsm}.wa(B_M).*B_M_w);

                        proj_M = proj_Wm.'*DD;
                        proj_Mt = proj_M*proj_Wm*dti;
                        proj_Mt_k_bar = proj_M*WM_k_bar(B_M)*dti;
end

