function u = single_dis(x,a,b)
N = length(x);
u = zeros(N,1);
for ii = 1:N
    if (-1 <= x(ii)) && (x(ii) <= 0)
        u(ii) = a;
    else
        u(ii) = b;
    end
end