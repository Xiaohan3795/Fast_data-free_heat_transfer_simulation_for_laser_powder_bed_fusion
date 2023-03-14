function [Q, pj_time] = make_gauss_orth(msh,utm1,rtm1,rt,ind,r,eta,u_pre,n_mu,um,d2)

    tic;
    x0m1=rtm1(1); y0m1=rtm1(2); z0m1=rtm1(3);
    x0=rt(1); y0=rt(2); z0=rt(3);

    mshvtx = msh.vtx(ind,:);

    inda = find(utm1>um);
    %Normalise utm1 to max 1
    utm1n = utm1/max(utm1);


    %The regression data
    data = mshvtx(inda,:); ydata = -log(utm1n(inda));
    A = [(data(:,1)-x0m1).^2, (data(:,2)-y0m1).^2, (data(:,3)-z0m1).^2]/2;
    rho = (A.'*A)\A.'*(ydata);
    sigmas = ones(3,1)./sqrt(rho);

    %Pick some values per sigma_dir centred at the estimator above
    logs = log(sigmas);

    for i=1:3
        Sigs(i,:)= exp(logs(i)*eta);
    end


    %Get a range of angles in [0,pi)
    theta = 0;

    %Initialise the basis with constant function
    Phi=ones(size(utm1));


    for i=1:length(theta)

        R_theta = [cos(theta(i)), sin(theta(i)), 0; -sin(theta(i)), cos(theta(i)), 0; 0, 0, 1];

        vtx_th_i = mshvtx*R_theta;


        for jx = 1:size(Sigs,2)
            for jy = 1:size(Sigs,2)
                for jz = 1:size(Sigs,2)

                    a=0.5/Sigs(1,jx)^2; 
                    b=0.5/Sigs(2,jy)^2; 
                    c=0.5/Sigs(3,jz)^2;

                    for imu = 0:n_mu-1
                        Phi = [Phi, exp(-(a * (vtx_th_i(:,1) - (x0+imu*d2)).^2 + ...
                                b * (vtx_th_i(:,2) - y0).^2 + ...
                                c * (vtx_th_i(:,3) - z0).^2))];
                    end

                end
            end
        end
    end

    %Add the normalised utm1 & orthogonalise.
    Q = rsvd_subGaussian([u_pre,Phi],r,1); 
    
    pj_time = toc;
    
end


