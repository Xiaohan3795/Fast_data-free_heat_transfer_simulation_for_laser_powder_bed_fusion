function [proj_K,proj_barK] = proj_sketch_K(msh,ut1_k,B_K,B_K_w,proj_Wk,barWk,nsk)

                        ux=[];
                        BK = [];
                        BK_w = [];
                        nsimp = size(msh.simp,1);
                        for q=1:length(msh.Gq{nsk}.weights)
                            ux =[ux; dot(msh.Gq{nsk}.Gq_temp_compute{q}, ut1_k(msh.simp'))'];
                            BK = [BK;B_K(:,2)+(q-1)*nsimp];
                            BK_w = [BK_w;B_K_w];
                        end
                        [kappa] = kappa_find(ux(BK));
                        kappa = kappa*1e-3;
                        K = reshape(kappa.*msh.Gq{nsk}.wa(BK).*BK_w, length(B_K(:,2)), length(msh.Gq{nsk}.weights)); %ais
                        neK = size(K,1);
                        K=sparse(1:neK,1:neK,sum(K,2));
                        
                        proj_barK = proj_Wk.'*K*barWk(B_K(:,1));
                        proj_K = proj_Wk.'*K*proj_Wk;
end

