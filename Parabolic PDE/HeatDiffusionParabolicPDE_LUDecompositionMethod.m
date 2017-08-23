function HeatDiffusionParabolicPDE_LUDecompositionMethod(L,T,maxk,n,cond)
% HeatDiffusionParabolicPDE_LUDecompositionMethod: Heat diffusion parabolic PDE 
% solved by the Crank-Nicholson method. Heat Diffusion in one dimensional wire.
%
% SINTAXIS:
%   HeatDiffusionParabolicPDE_LUDecompositionMethod(L,T,maxk,n,cond)
%
%     L : Length of the wire.
%     T : Final time.
%  maxk : Number of time steps.
%     n : Number of space steps.
%  cond : Conductivity.
%
% EXAMPLE:
%  HeatDiffusionParabolicPDE_LUDecompositionMethod(1,1,2500,50,0.25)
%

% Parameters needed to solve the equation within the LU decomposition method.
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

% Auxiliar left and right matrixes.
aal(1:n-2) = -b;
bbl(1:n-1) = 2 + 2.*b;
ccl(1:n-2) = -b;
MMl = diag(bbl,0) + diag(aal,-1) + diag(ccl,1);
[L,U] = lu(MMl);

aar(1:n-2) = b;
bbr(1:n-1) = 2 - 2.*b;
ccr(1:n-2) = b;
MMr = diag(bbr,0) + diag(aar,-1) + diag(ccr,1);

% Implementation of the LU decomposition method.
for k = 2:maxk
    uu = u(2:n,k-1);
    vv = MMr*uu;
    w(1) = vv(1);
    for j = 2:n-1
        w(j) = vv(j) - L(j,j-1)*w(j-1);
    end
    u(n,k) = w(n-1)/U(n-1,n-1);
    for i = n-1:-1:2
        u(i,k) = (w(i-1)-U(i-1,i)*u(i+1,k))/U(i-1,i-1);
    end
end

% Graphical representation of the temperature.
figure(1)
plot(x,u(:,1),'-',x,u(:,100),'-',x,u(:,300),'-',x,u(:,600),'-')
title('Temperature within the LU decomposition method')
xlabel('X')
ylabel('Y')

figure(2)
mesh(x,time,u')
title('Temperature within the LU decomposition method')
xlabel('X')
ylabel('Temperature')

