% Defind Exact solution for linear advection equation
function u = uexact(x,t,uL,uR)
N = length(x);
for ii = 1:N
    if (uL >= uR)
        s = (uL+uR)/2;      % Consider for Burger's equation
        if (x(ii) < s*t)
            u(ii) = uL;
        else u(ii) = uR;
        end
    else
        if (x(ii) < uL*t)
            u(ii) = uL;
        elseif (x(ii) > uR*t)
            u(ii) = uR;
        else u(ii) = x(ii)/t;
        end
    end
end