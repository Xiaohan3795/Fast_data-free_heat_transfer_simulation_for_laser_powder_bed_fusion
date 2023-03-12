function [rts,l] = load_vector(msh, tn,vn,dt,rx,lx,ly,d2,fm)
% the load vector l and Gaussian heat flux function f
   
    
        Gam1ind = find(msh.csrf(:,3)> max(msh.vtx(:,3))-0.00001);

        rts = [];
        for t=0:tn
            sn = vn*dt*t;
            x_d = mod(sn,rx-lx);      
            f = [lx+x_d, ly,  max(msh.vtx(:,3)), d2, fm];  
            rts= [rts; f(1:3)];
            % moving heat source
            l{t+1} = make_surface_source(msh,f,Gam1ind);
        end
        
end

