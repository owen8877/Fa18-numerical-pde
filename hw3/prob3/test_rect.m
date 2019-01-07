clear; %clc

N = 128; h = 1/N;

xs = (0:N)/N; ys = xs';
[X, Y] = meshgrid(xs, ys);

%
% u_truth = ((1-X).*(1-Y)).^2.*Y;
% u0 = ys.*(1-ys).^2;
% f = -(6*Y-4).*((1-X).^2)-2.*(Y.*(1-Y).^2);
% beta = 0;
% g_x0 = - (1-xs).^2 + beta * u_truth(1, :);
% g_1y = 0*ys + beta * u_truth(:, end);
% g_x1 = 0*xs + beta * u_truth(end, :);

%
% u_truth = cos(X*pi/2).*(1-Y).^2.*Y;
% u0 = ys.*(1-ys).^2;
% f = -(6*Y-4).*cos(X*pi/2)+pi^2/4.*cos(X*pi/2).*(1-Y).^2.*Y;
% beta = 0;
% g_x0 = - cos(xs*pi/2) + beta * u_truth(1, :);
% g_1y = (-pi/2).*(1-ys).^2.*ys + beta * u_truth(:, end);
% g_x1 = 0*xs + beta * u_truth(end, :);

%
u_truth = sin(X*2*pi).*sin(Y*2*pi) + (1-X).*Y.*(1-Y);
u0 = ys.*(1-ys);
f = (8*pi^2)*sin(X*2*pi).*sin(Y*2*pi) + 2*(1-X);
beta = 0;
g_x0 = -(sin(xs*2*pi)*2*pi + (1-xs)) + beta * u_truth(1, :);
g_1y = 2*pi*sin(ys*2*pi)-ys.*(1-ys) + beta * u_truth(:, end);
g_x1 = sin(xs*2*pi)*2*pi - (1-xs) + beta * u_truth(end, :);

fv = ffunc(f, u0, g_x0, g_1y, g_x1, const, h, N);

options.maxItr = 1e3;
options.outputInt = 1e1;
options.torlerance = 1e-12;

Au_truth = Afunc(u_truth(:, 2:end), const, beta, h, N);

u = CG(zeros(N+1, N), @(u) Afunc(u, const, beta, h, N), fv, options);
mesh(u)

