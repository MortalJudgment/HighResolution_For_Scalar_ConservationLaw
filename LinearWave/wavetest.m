function u = wavetest(x,a,z,delta,alpha,beta)
N = length(x);
u = zeros(N,1);
for ii = 1:N
    if (-0.8 <= x(ii)) && (x(ii) <= -0.6)
        u(ii) = 1/6*(G(x(ii),beta,z-delta) + G(x(ii),beta,z+delta) + 4*G(x(ii),beta,z));
    elseif (-0.4 <= x(ii)) && (x(ii) <= -0.2)
        u(ii) = 1;
    elseif (0 <= x(ii)) && (x(ii) <= 0.2)
        u(ii) = 1 - abs(10*(x(ii)-0.1));
    elseif (0.4 <= x(ii)) && (x(ii) <= 0.6)
        u(ii) = 1/6*(F(x(ii),alpha,a-delta) + F(x(ii),alpha,a+delta) + 4*F(x(ii),alpha,a));
    else
        u(ii) = 0;
    end
end