function [proj_R,proj_barR] = proj_sketch_R(msh,ut1_k,B_S,sigma,varepsilon,B_S_w,proj_Wr,barWr)

                        ux=[];
                        for q=1:length(msh.Gqsrf{5}.weights)
                            ux =[ux; dot(msh.Gqsrf{5}.Gqsrf_temp_compute{q}, ut1_k(msh.srf'))'];
                        end
                        neS4 = length(B_S);
                        S4_ = sparse(1:neS4,1:neS4,sigma*varepsilon*ux(B_S).^3.*msh.Gqsrf{5}.wa(B_S).*B_S_w);

                        proj_R =  proj_Wr.'* S4_ * proj_Wr;
                        proj_barR = proj_Wr.'* S4_* barWr(B_S);
end

