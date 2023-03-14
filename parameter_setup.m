function [dim_Psi,eta,fm,h,u_a,bpt,d2,sigma,varepsilon,e,n_mu,n_u,n_t,lx,ly,rx,um,v,dt,dti,tn,nsk,nsm,run_max,sp,skm,skk,sks] = parameter_setup()

    dim_Psi = 80; % the dimension of projection basis
    eta = [0.5, 0.8, 1, 1.2, 1.5]; % the scale factor of the variances of Gaussian functions
    
    sp = 0.6; % percentage of columns to approximate leverage scores
    
    %sketching parameters
    skm = 4000;
    skk = 4000;
    sks = 2500;
    
    fm = 45000; % peak value of Gaussian heat source
    h = 10e-6;  % convection coefficient
    u_a = 20;  % ambient temperature(C)
    bpt = 20;  % temperature of build platform(C)
    d2 = 50e-3; %(mm)
    sigma = 5.67e-14; %W/(mm^2K^4) Stefan-Boltzmann constant
    varepsilon = 0.04; %emissity
    
    e = 1e-5;% picard iteration error tolerance
    
    n_mu = 2; % the number of means for selected Gaussian functions
    n_u = 5; % the number of previous temperatures
    n_t = 5; % the number of time steps for FOM
    
    lx = 0.3; % x coordinate of the start position
    ly = 0.3; % y coordinate of the start position
    rx = 1;  % x coordinate of the end position
    
    um = 30; % the temperature threshold to select nodes for linear regression
    v = 400; % scan speed
    
    dt = 5e-5; % discrete time interval dt(s)
    dti=1/dt; 
    tn=14; % total amount of time steps

    nsk = 3; % Gaussian quadrature rules of stiffness matrix
    nsm = 4; % Gaussian quadrature rules of mass matrix

    run_max = 30; % maximum number of runs
end