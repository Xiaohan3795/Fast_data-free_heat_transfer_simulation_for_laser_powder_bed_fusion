clear all;
clc;

load simple_case;

% first run FOM, then run surrogate
model_str{1} = 'FOM';
model_str{2} = 'SURROGATE'; 

% set up parameters
[dim_Psi,eta,fm,h,u_a,bpt,d2,sigma,varepsilon,e,n_mu,n_u,n_t,lx,ly,rx,um,v,dt,dti,tn,nsk,nsm,run_max,sp,skm,skk,sks] = parameter_setup();

% constant matrices and vectors in simulation
[bpv,bpv_dif,C,a,Phim,Phir,barPhir,barPhik,Phik] = simulation_prepare(msh,h,u_a,sigma,varepsilon,nsm,bpt);

% the positions of laser beam ceter and load vector
[rts,l] = load_vector(msh, tn,v,dt,rx,lx,ly,d2,fm);

% record all temperatures
temp_FOM = [];
temp_S = [];
    
% initial condition
u0 = u_a*ones(size(msh.vtx,1),1);  
u_t=u0;
u_t(bpv)=bpt;

for model_i = 1:2
    model = model_str{model_i};
    
    for t=0:tn

        t

        % run FOM for the first n_t time steps
        if strcmp(model, 'SURROGATE') & t>=n_t
            model_t = 'SURROGATE';
            Phim_k_bar = Phim*u_t(bpv_dif);
        else
            model_t = 'FOM';
        end

        switch model_t
            case 'FOM'
                %load vector and constant right hand side vectors
                lf = l{t+1}(bpv_dif)+a(bpv_dif)-C(bpv_dif,bpv)*(ones(length(bpv),1)*bpt);

            case 'SURROGATE'
                % the generation of projection basis
                [Psi, pj_time] = make_gauss_orth(msh,temp_S(bpv_dif,t),rts(t,:),rts(t+1,:),bpv_dif,dim_Psi,eta,temp_S(bpv_dif,t-(n_u-1):t),n_mu,um,d2); 
                proj_basis_time(t+1) = pj_time;

                % sketching and projection
                [proj_C,proj_Wm,proj_Wr,proj_Wk,proj_lf,B_M,B_M_w,B_K,B_K_w,B_S,B_S_w,p_time,s_time] = pre_proj_sketch_new(Psi,a,C,bpv,bpt,Phim,Phir,Phik,skm,skk,sks,bpv_dif,l{t+1},sp);
                sketch_time(t+1) = s_time;
                proj_time(t+1) = p_time;

            otherwise
                error('wrong model selection');
        end

        % initial ut1_0
        ut1_k = u_t;
        error = e+1;
        run = 0;

        % Picard iterations
        while error>e

            switch model_t
                case 'FOM'
                    tic;
                    % diagonal temperature-dependent matrices
                    [Dr,Dk,Dm] = D_matrices(nsk, nsm, msh,ut1_k,sigma,varepsilon);
                    % radiation heat loss
                    R =  Phir.' * Dr * Phir;
                    barR = Phir.'* Dr* barPhir;
                    % stiffness matrix
                    K = Phik.'*Dk*Phik;
                    barK = Phik.'*Dk*barPhik;
                    % mass matrix
                    M = Phim.'*Dm*Phim;
                    Mt = M*dti;
                    % compute the high-dimensional temperature
                    A = Mt+K+C(bpv_dif,bpv_dif)+R;
                    b = lf + Mt * u_t(bpv_dif)- barK-barR;
                    ut1_k1_ = A\b;
                    % time cost of Picard iteration
                    trecord_FOM(run+1, t+1) = toc;

                case 'SURROGATE'
                    tic;
                    % proj_R  proj_barR
                    [proj_R,proj_barR] = proj_sketch_R(msh,ut1_k,B_S,sigma,varepsilon,B_S_w,proj_Wr,barPhir);
                    % proj_K  proj_barK
                    [proj_K,proj_barK] = proj_sketch_K(msh,ut1_k,B_K,B_K_w,proj_Wk,barPhik,nsk);
                    % proj_Mt  proj_Mt_k_bar
                    [proj_Mt,proj_Mt_k_bar] = proj_sketch_M(msh,ut1_k,B_M,B_M_w,proj_Wm,nsm,dti,Phim_k_bar);
                    % compute the low-dimensional temperature
                    A = proj_Mt+proj_K+proj_C+proj_R;
                    b = proj_lf + proj_Mt_k_bar - proj_barK-proj_barR;
                    ut1_k1_ = A\b;
                    % reconstruct temperature
                    ut1_k1_ = Psi*ut1_k1_;
                    % time cost of Picard iteration
                    trecord_S(run+1, t+1) = toc;

                otherwise
                    error('wrong model selection');
            end
            % compute the error
            ut1_k1(bpv_dif,:) = ut1_k1_; ut1_k1(bpv,:)=bpt;
            error = norm(ut1_k1 - ut1_k)/norm(ut1_k);

            % number of iterations
            run = run+1;
            %maximum number of runs
            if run>run_max
                break;
            end
            ut1_k = ut1_k1;

        end

        u_t=ut1_k1; 
        
        % record temperatures
        switch model
            case 'FOM'
                temp_FOM = [temp_FOM,u_t];
            case 'SURROGATE'
                temp_S = [temp_S,u_t];
        end
        

    end
end
       



%% 2 norm relative error of temperatures

for i = n_u+1:tn+1
     re_2norm(i-n_u) = norm(temp_FOM(:,i)-temp_S(:,i))/norm(temp_FOM(:,i))*100;
end

re_2norm


%% time cost reduction

% time cost of FOM
time_FOM = sum(trecord_FOM(:,n_u+1:end));
% time cost of surrogate
time_surrogate = sum(trecord_S(:,n_u+1:end))+sketch_time(n_u+1:end)+proj_time(n_u+1:end)+proj_basis_time(n_u+1:end);
% percentage of time cost reduction
time_reduction = (time_FOM-time_surrogate)./time_FOM*100;

time_reduction

    
    




