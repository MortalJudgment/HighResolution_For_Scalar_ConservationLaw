function [ numflux ] = BurgersRoe ( u , v , lambda , maxvel )
% Roe numerical flux for Burger equation .
% No sonic fix .
fu = u.^2 ;
fv = v.^2 ;
alpha = u+v ;
numflux = (alpha >= 0).*fu + (alpha < 0).*fv ;
end