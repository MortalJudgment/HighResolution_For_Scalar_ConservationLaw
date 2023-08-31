function [ numflux ] = BurgersLW ( u , v , lambda , maxvel )
% Lax Wendroff numerical flux for Burger equation 
fu = u.^2 ; 
fv = v.^2 ;
alpha = lambda*2*(u+v)/2 ;
numflux = ( fu + fv )/2 - alpha/2.*( fv-fu ) ;
end