clc
close all
clear all

L = 2 ; 
FinalTime = 4 ; 
N = 2^8 ; 
h = L/N; 
CFL = 0.9;

x = [-5:h:5]';
[u1] = u0(x) ;

%%-----------------------------------------------------------------------%%
% Solve Problem

time = 0 ; 
tstep = 0 ;
% Set the time step
Nx = length(x);
k = CFL*h ;
maxvel = 1;
% Integrate scheme
u_xx = u1;

% Choosing flux-limiter scheme
% 0: non ; 1: CO; 2: Koren ; 3: Sweby ; 4: OSPRE ; 5: Van Leer
type = 2 ; 
beta = 1 ;

% Start a loop
while(time < FinalTime)
    if( time + k > FinalTime ) 
        k = FinalTime - time ; 
    end
    
    % Boudary condition
    % Periodic boundary conditions
    [ xe , ue ] = extend ( x , u_xx , h , 1 , 'P' , 0 , 'P' , 0 ) ;
%-------------------------------------------------------------------------%

    %------------------- Normal case ------------------%
    % Compute numerical flux
    % Lax-Friedrich scheme
%     du = -(LinwaveLF(ue(2:Nx+1), ue(3:Nx+2), k/h , maxvel ) - ...
%         LinwaveLF(ue(1:Nx), ue(2:Nx+1), k/h , maxvel ) )/h ;
    
    % Lax-Wendroff scheme
%     du = -(LinwaveLW(ue(2:Nx+1), ue(3:Nx+2), k/h , maxvel ) - ...
%         LinwaveLW(ue(1:Nx), ue(2:Nx+1), k/h , maxvel ) )/h ;
    %-------------- High Resolution case --------------%
    
    % Compute indicator function
  
    r = ( ue(2:Nx+1) - ue(1:Nx) )./( ue(3:Nx+2) - ue(2:Nx+1) ) ;
    [ xe , re ] = extend( x,r,h,1,'N',0,'N',0) ;
    phiL = FluxLimit ( re(1:Nx), type, beta ) ;
    phiR = FluxLimit ( re(2:Nx+1), type, beta ) ;
    
    % Compute left flux in cell
    Fluxlow = LinwaveLF( ue(1:Nx), ue(2:Nx+1), k/h , maxvel ) ;
    Fluxhigh = LinwaveLW( ue(1:Nx), ue(2:Nx+1) ,k/h , maxvel ) ;
    FluxL = Fluxlow + phiL.*( Fluxhigh - Fluxlow ) ;
    
    % Compute right flux in cell
    Fluxlow = LinwaveLF( ue(2:Nx+1), ue(3:Nx+2), k/h, maxvel ) ;
    Fluxhigh = LinwaveLW( ue(2:Nx+1), ue(3:Nx+2), k/h, maxvel ) ;
    FluxR = Fluxlow + phiR.*( Fluxhigh - Fluxlow ) ;
    
    % Compute RHS
    du = -(FluxR - FluxL )/h ;
%-------------------------------------------------------------------------%
    % Update solution
    u_xx = u_xx + k * du ;

    time = time + k ; 
    tstep = tstep + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%---------------------------   Plot   ---------------------------------%%
% Exact solution
uex = uexact(x,FinalTime,maxvel);

error = abs(u_xx-uex);
errorMax = max(error)

hold on
plot(x,uex,'r x-')
plot(x,u_xx,'b x-')
legend('u exact', 'u MinMod flux-limiter')
title(['Comparison between exact solution and numerical solution after ', num2str(FinalTime), ' seconds'])
% axis ([2.5 4.5 -0.5 1.5])
axis ([2 5 -0.5 1.5])

hold off