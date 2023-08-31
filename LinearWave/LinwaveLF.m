function[ numflux ] = LinwaveLF ( u , v , lambda , maxvel )
% Lax-Friedrich numerical flux for wave equation
numflux = (u + v) /2 - maxvel /2*(v - u) ;
end