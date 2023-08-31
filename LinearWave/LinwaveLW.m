function[ numflux ] = LinwaveLW(u,v,lambda,maxvel)
% Lax-Wendroff numerical flux for wave equation
% numflux = (u + v) /2 - lambda / 2*(v - u) ;
numflux = u + 1/2*(1-lambda)*(v-u);
end