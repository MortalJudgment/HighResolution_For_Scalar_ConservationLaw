function u = single_dis(x,k,uL,uR)
% k : position where dis.
% uL : value left.
% uR : value right
N = length(x);
u = zeros(N,1);
for ii = 1:N
    if (k <= x(ii))
        u(ii) = uR;
    else
        u(ii) = uL;
    end
end