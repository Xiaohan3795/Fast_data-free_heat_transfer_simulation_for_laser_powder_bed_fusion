function [proj_C,proj_Wm,proj_Wr,proj_Wk,proj_lf,B_M,B_M_w,B_K,B_K_w,B_S,B_S_w,proj_time,sketch_time] = pre_proj_sketch_new(Psi,a,C,bpv,bpt,Wm,Wr,Wk,skm,skk,sks,bpv_dif,l,sp)


            
            %randomized sketching
            tic;
            [B_M,B_M_w] = sketch_index_noproj(Wm,Psi,skm,sp);
            [B_K_dif, B_K_w] = sketch_index_noproj(Wk,Psi,skk,sp); 
            B_K(:,1) = B_K_dif;
            B_K(:,2) = ceil(B_K_dif/3);
            [B_S, B_S_w] = sketch_index_noproj(Wr,Psi,sks,sp); 
            sketch_time = toc;
            
            
          
            % projection of the selected rows
            tic;
            proj_C = Psi.'*C(bpv_dif,bpv_dif)*Psi;
            proj_Wm = Wm(B_M,:)*Psi;
            proj_Wr =Wr(B_S,:)*Psi;
            proj_Wk = Wk(B_K(:,1),:)*Psi;
            proj_lf = Psi.'*(l(bpv_dif) + a(bpv_dif)-C(bpv_dif,bpv)*(ones(length(bpv),1)*bpt));
            proj_time = toc;
            


 

            
            

end

