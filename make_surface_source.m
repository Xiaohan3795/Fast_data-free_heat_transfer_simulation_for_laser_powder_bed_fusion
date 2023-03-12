function [Sq] = make_surface_source(msh,f,Gamma1)

%This function evaluates the load vector \int_\Omega f \phi_i(x) assuming 
%linear basis functions, and f is a Gaussian function with ampliture mag 
%centred at (fx,fy,fz) with std fs.

%Input: (1) mesh geometry including the gradients operator. This assumes
%           msh processed by [msh] = analyse_mesh(vtx,simp)
%       (2) Quadrature points and weights structure. This is 3 matrices 
%           in nn x np dimensions with the x-, y-, and z- coorfinates of
%           the integration points, and a vector in np dimension with the
%           relevant weights.
%       (3) Constants h and u)a
%       (4) Index Gamma1 into msh.srf to indicate which boundary elements
%       are on \Gamma1 for the surface integrals.

Sq = zeros(size(msh.vtx,1),1);
n_q = 3;
Phi_qc = Gaussian_quadrature_phi_srf(msh,msh.Gqsrf{n_q},Gamma1);

Qs = zeros(size(msh.srf,1),1);
for q = 1:length(msh.Gqsrf{n_q}.weights)
    p = msh.Gqsrf{n_q}.coord(Gamma1,:,q);
    Q = f(5)*exp(-2*( ((f(1)-p(:,1)).^2)/f(4)^2 + ((f(2)-p(:,2)).^2)/f(4)^2 + ((f(3)-p(:,3)).^2)/f(4)^2));
    Qs(Gamma1) = Q.*msh.srfarea(Gamma1);
    Sq = Sq+msh.Gqsrf{n_q}.weights(q)*Phi_qc{q}.'*Qs;
end


end

