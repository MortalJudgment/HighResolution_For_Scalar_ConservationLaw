clc
clear all
close all

L = 2 ; 
FinalTime = 0.2 ; 
N = 2^8 ; 
h = L/N; 
CFL = 0.9;

% Domain 
x = [0:h:2]' ;
% IC
% k : position where dis.
% uL : value left.
% uR : value right
k = 1/4;
uL = 1;
uR = 2;

[u1] = u0(x,k,uL,uR) ;

%%-----------------------------------------------------------------------%%
% Solve Problem

time = 0 ; 
tstep = 0 ;
% Set the time step
Nx = length(x);
maxvel = max( 2*abs(u1) ) ;
k = CFL*h/maxvel ;
% Integrate scheme
u_xx = u1;

% Choosing flux-limiter scheme
% 0: non ; 1: CO; 2: Koren ; 3: Sweby ; 4: OSPRE ; 5: Van Leer
type = 3 ; 
beta = 2;

% Start a loop
while(time < FinalTime)
    if( time + k > FinalTime ) 
        k = FinalTime - time ; 
    end
    
    % Boudary condition
    % Periodic boundary conditions
    [ xe , ue ] = extend( x , u_xx , h , 1 , 'P' , 0 , 'P' , 0 ) ;
%-------------------------------------------------------------------------%
 
    %------------------- Normal case ------------------%
    % Compute numerical flux
    % Lax-Friedrich scheme
%     du = - ( BurgersLF(ue(2:Nx+1),ue(3:Nx+2),0,maxvel )...
%         - BurgersLF(ue(1:Nx),ue(2:Nx+1),0,maxvel) )/h ;
    % Lax-Wendroff scheme
%     du = -(BurgersLW(ue(2:Nx+1), ue(3:Nx+2), k/h , maxvel ) - ...
%         BurgersLW(ue(1:Nx), ue(2:Nx+1), k/h , maxvel ) )/h ;
    %-------------- High Resolution case --------------%
    
    % Compute indicator function
    r = ( ue(2:Nx+1) - ue(1:Nx) )./( ue(3:Nx+2) - ue(2:Nx+1) );
    size(r)
    [ xe , re ] = extend ( x , r , h , 1 , 'N' , 0 , 'N' , 0 ) ;
    rm = 1./re ;
    phiLp = FluxLimit( re(1:Nx),type,beta ) ;
    phiRp = FluxLimit( re(2:Nx+1),type,beta ) ;
    phiLm = FluxLimit( rm(2:Nx+1),type,beta ) ;
    phiRm = FluxLimit( rm(3:Nx+2),type,beta ) ;
    ufilt = (u_xx>=0) ;
    phiL = ufilt.*phiLp + ( 1-ufilt ).*phiLm ;
    phiR = ufilt.*phiRp + ( 1-ufilt ).*phiRm ;
    % Compute left flux - Change numerical flux here
    Fluxlow = BurgersLF( ue(1:Nx),ue(2:Nx+1),0,maxvel ) ;
    Fluxhigh = BurgersLW( ue(1:Nx),ue(2:Nx+1),k/h,maxvel ) ;
    FluxL = Fluxlow + phiL.*( Fluxhigh-Fluxlow ) ;
    % Compute right flux - Change numerical flux here
    Fluxlow = BurgersLF( ue(2:Nx+1),ue(3:Nx+2),0,maxvel ) ;
    Fluxhigh = BurgersLW( ue(2:Nx+1),ue(3:Nx+2),k/h,maxvel ) ;
    FluxR = Fluxlow + phiR.*( Fluxhigh-Fluxlow ) ;
    % Compute RHS
    du = -( FluxR-FluxL )/h ;
%-------------------------------------------------------------------------%
    % Update solution
    u_xx = u_xx + k * du ;

    time = time + k ; 
    tstep = tstep + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---------------------------   Plot   ---------------------------------%%
% Exact solution
uex = uexact(x,FinalTime*maxvel/2 + 0.167,uL,uR);

% error = abs(u_xx-uex);
% errorMax = max(error)
% plot(x,u1,'r--')
hold on
plot(x,uex,'r x-')
plot(x,u_xx,'b x-')
% legend('u exact', 'u Superbee flux-limiter')
axis([0 2 0.85 2.15])
% hold off