function HeatDiffusionParabolicPDE_ExplicitMethod(L,T,maxk,n,cond)
% HeatDiffusionParabolicPDE_ExplicitMethod: Heat diffusion parabolic PDE 
% solved by the explicit method. Heat Diffusion in one dimensional wire.
%
% SINTAXIS:
%   HeatDiffusionParabolicPDE_ExplicitMethod(L,T,maxk,n,cond)
%
%     L : Length of the wire.
%     T : Final time.
%  maxk : Number of time steps.
%     n : Number of space steps.
%  cond : Conductivity.
%
% EXAMPLE:
%  HeatDiffusionParabolicPDE_ExplicitMethod(1,1,2500,50,0.25)
%

% Parameters needed to solve the equation within the explicit method.
dt = T/maxk;            % Time steps delta.
dx = L/n;               % Space steps delta.
b = 2.*cond*dt/(dx*dx); % Stability parameter (b <= 1).

% Initial temperature of the wire.
for i = 1:n+1
    x(i) = (i-1)*dx;
    u(i,1) = sin(pi*x(i));
end

% Temperature at the boundary.
for k = 1:maxk+1
    u(1,k) = 0.0;
    u(n+1,k) = 0.0;
    time(k) = (k-1)*dt;
end

% Implementation of the explicit method.
for k = 1:maxk
    for i = 2:n
        u(i,k+1) = u(i,k) + 0.5*b*(u(i-1,k)+u(i+1,k)-2.*u(i,k));
    end
end

% Graphical representation of the temperature.
figure(1)
plot(x,u(:,1),'-',x,u(:,100),'-',x,u(:,300),'-',x,u(:,600),'-')
title('Temperature within the explicit method')
xlabel('X')
ylabel('Y')

figure(2)
mesh(x,time,u')
title('Temperature within the explicit method')
xlabel('X')
ylabel('Temperature')

