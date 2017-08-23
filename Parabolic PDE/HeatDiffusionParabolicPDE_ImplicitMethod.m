function HeatDiffusionParabolicPDE_ImplicitMethod(L,T,maxk,n,cond)
% HeatDiffusionParabolicPDE_ImplicitMethod: Heat diffusion parabolic PDE 
% solved by the implicit method. Heat Diffusion in one dimensional wire.
%
% SINTAXIS:
%   HeatDiffusionParabolicPDE_ImplicitMethod(L,T,maxk,n,cond)
%
%     L : Length of the wire.
%     T : Final time.
%  maxk : Number of time steps.
%     n : Number of space steps.
%  cond : Conductivity.
%
% EXAMPLE:
%  HeatDiffusionParabolicPDE_ImplicitMethod(1,1,2500,50,0.25)
%

% Parameters needed to solve the equation within the implicit method.
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

% Auxiliar matrix.
aa(1:n-2) = -b;
bb(1:n-1) = 1 + 2.*b;
cc(1:n-2) = -b;
MM = inv(diag(bb,0) + diag(aa,-1) + diag(cc,1));

% Implementation of the implicit method.
for k = 2:maxk
    uu = u(2:n,k-1);
    u(2:n,k) = MM*uu;
end

% Graphical representation of the temperature.
figure(1)
plot(x,u(:,1),'-',x,u(:,100),'-',x,u(:,300),'-',x,u(:,600),'-')
title('Temperature within the implicit method')
xlabel('X')
ylabel('Y')

figure(2)
mesh(x,time,u')
title('Temperature within the implicit method')
xlabel('X')
ylabel('Temperature')

