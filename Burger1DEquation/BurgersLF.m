function [ numflux ] = BurgersLF ( u , v , lambda , maxvel )
%Lax Friedrich numerical flux for Burger equation 
fu = u.^2 ;
fv = v.^2 ;
numflux = ( fu + fv )/2 - maxvel/2*(v-u) ;
end