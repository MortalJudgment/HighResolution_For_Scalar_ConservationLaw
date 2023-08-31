function [ xe , ue ] = extend( x , u , h ,m, BCl , ul , BCr , ur )

% Boundary conditions (BC).
% Stand for (sf): "D": Dirichlets; "N" : Neumann ; "P " : periodic

%   x: node - defind for dicreted domain
%   u: solution
%   h: space step size
%   m: number cell extend (ie, 1)
%   BCl: boundary condition on left side
%   BCr: boundary condition on right side
%   ul: BC value on left side
%   ur: BC value on right side
%   whereas ul & ur : BC value only for Dirichlet BC

%% --------------------------------------------------------------------- %%
% Defind basic value: 
xl = min(x);            % Left node
xr = max(x);            % Right node
N = length(u);          % Total value
%--------------------------------%
xe = zeros(N+2*m,1);    % node extended
ue = zeros(N+2*m,1);    % value extended
q = [1:m];

% Extend x
xe(m-q+1)= xl - q*h;        % Left ghost point
xe(N+m+q)= xr + q*h;        % Right ghost point
xe((m+1):(N+m)) = x(1:N);   % Save the old value

% Periodic BC
if (BCl == 'P') || (BCr == 'P')
    ue(m-q+1) = u(N-q) ;        % Left value
    ue(N+m+q) = u(q+1) ;        % Right value
    ue((m+1):(N+m)) = u(1:N) ;  % Old value
end
% Dirichlet BC
if (BCl == 'D')
    ue(m-q+1) = -u(q+1) + 2*ul;
else
    ue(m-q+1) = u(q+1) ;
end

if (BCr == 'D')
    ue(N+m+q) = - u(N-q) + 2*ur ;
else
    ue(N+m+q) = u(N-q) ;
end
ue ((m+1):(N+m)) = u(1:N) ;
end
