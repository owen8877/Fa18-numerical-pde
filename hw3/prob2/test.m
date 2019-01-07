clear; %clc

% name = 'Dirichlet - C^\infty';
% a = 1;
% u0 = @(xs) 0*xs;
% v0 = @(xs) sin(pi*xs) * pi;
% alpha0 = @(x, t) 1; alpha1 = @(x, t) 1;
% beta0 = @(x, t) 0; beta1 = @(x, t) 0;
% g0 = @(t) 0; g1 = @(t) 0;
% history_truth = @(X, T) sin(pi*X) .* sin(pi*T);

% name = 'Neumann - C^\infty';
% a = 1;
% u0 = @(xs) 0*xs;
% v0 = @(xs) cos(pi*xs) * pi;
% alpha0 = @(x, t) 0; alpha1 = @(x, t) 0;
% beta0 = @(x, t) 1; beta1 = @(x, t) 1;
% g0 = @(t) 0; g1 = @(t) 0;
% history_truth = @(X, T) cos(pi*X) .* sin(pi*T);

name = 'Neumann - C^0';
a = 1;
% u0 = @(xs) circshift(max(0, min(xs, 0.5-xs)), (numel(xs)-1)/4);
% v0 = @(xs) circshift(2*(0.5>xs & xs>eps) .* sign(xs-0.25), (numel(xs)-1)/4);
% u0 = @(xs) circshift(max(0, min(xs, 0.5-xs)), 0);
% v0 = @(xs) circshift(2*(0.5>xs & xs>eps) .* sign(xs-0.25), 0);
u0 = @(xs) max(0, min(xs-0.25, 0.5));
v0 = @(xs) (xs>0.25) & (xs<0.75);
alpha0 = @(x, t) 0; alpha1 = @(x, t) 0;
beta0 = @(x, t) 1; beta1 = @(x, t) 1;
g0 = @(t) 0; g1 = @(t) 0;

hs = 1./[8 16 32 64 128];
CFL = 1;

errors = hs*0;
for i = 1:numel(hs)
    h = hs(i);
    dt = h*CFL/a;
    xs = 0:h:1;
    ts = 0:dt:1;

    tic
    [u, history] = explicit_wave(ts, xs, u0, v0, a, alpha0, beta0, alpha1, beta1, g0, g1);
    toc
    [X, T] = meshgrid(xs, ts);
    if exist('history_truth', 'var')
        errors(i) = matNorm(history_truth(X, T)-history);
    end
    if abs(1/h-64) < eps
        figure(2); clf
        mesh(xs, ts, history);
        title([name '; h= ' num2str(h)])
        xlabel('x')
        ylabel('time')
    end
end
figure(1); clf; hold on;
plot(hs, errors, 's-')
title(name)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
grid on
xlabel('h')
ylabel('error')