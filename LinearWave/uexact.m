% Defind Exact solution for linear advection equation
function u = uexact(x,t,maxvel)
% "vel" stand for velocity

u = u0(x-maxvel*t);
end