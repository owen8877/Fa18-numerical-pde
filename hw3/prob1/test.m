clear; %clc

f = @(u) u.^2 / 2;
aFunc = @(u) u;

%%
CFL = 1;
h = 1e-2;
maxT = 0.2;

figure(1)
[u, history, X, T, xs, dump] = helper(@non_cons_upwind, h, maxT, CFL, f, aFunc);
subplot(2, 3, 1); mesh(X, T, history); title('Non conservative upwind'); xlabel('x'); ylabel('t')
subplot(2, 3, 4); contour(X, T, history); title('Non conservative upwind'); xlabel('x'); ylabel('t'); grid on
[u, history, X, T, xs, dump] = helper(@cons_upwind, h, maxT, CFL, f, aFunc);
subplot(2, 3, 2); mesh(X, T, history); title('Conservative upwind'); xlabel('x'); ylabel('t')
subplot(2, 3, 5); contour(X, T, history); title('Conservative upwind'); xlabel('x'); ylabel('t'); grid on
[u, history, X, T, xs, dump] = helper(@lax_wendroff, h, maxT, CFL, f, aFunc);
subplot(2, 3, 3); mesh(X, T, history); title('Lax Wendroff'); xlabel('x'); ylabel('t')
subplot(2, 3, 6); contour(X, T, history); title('Lax Wendroff'); xlabel('x'); ylabel('t'); grid on

%%
CFL = 1;
hs = [1e-1, 1e-3, 1e-4];
maxT = 1;

figure(2)
for j = 1:numel(hs)
    h = hs(j);
    [u, history, X, T, xs, dump1] = helper(@cons_upwind, h, maxT, CFL, f, aFunc);
    [u, history, X, T, xs, dump2] = helper(@lax_wendroff, h, maxT, CFL, f, aFunc);
    
    s = subplot(3, 2, 2*j-1); hold(s, 'on'); title(['Conservative upwind; Comparison at h=' num2str(h)])
    for l = 1:5
        plot(xs, dump1(l, :), '-', 'DisplayName', num2str(l/5))
    end
    title(legend('Location', 'best'), 't'); xlim([-0.1 0.7])
    s = subplot(3, 2, 2*j); hold(s, 'on'); title(['Lax Wendroff; Comparison at h=' num2str(h)])
    for l = 1:5
        plot(xs, dump2(l, :), '-', 'DisplayName', num2str(l/5))
    end
    title(legend('Location', 'best'), 't'); xlim([-0.1 0.7])
end

%%
h = 1e-2;
maxT = 0.2;
figure(3); hold on
[u, ~, ~, ~, xs, ~] = helper(@cons_upwind, h, maxT, 1, f, aFunc);
plot(xs, u, 'DisplayName', '1')
[u, ~, ~, ~, xs, ~] = helper(@cons_upwind, h, maxT, 0.5, f, aFunc);
plot(xs, u, 'DisplayName', '0.5')
[u, ~, ~, ~, xs, ~] = helper(@cons_upwind, h, maxT, 0.05, f, aFunc);
plot(xs, u, 'DisplayName', '0.05')

grid on; xlim([0.05 0.15]); title(legend('Location', 'best'), 'CFL');
xlabel('x'); ylabel('u'); title('Under different CFL settings at t = 0.2')

%%
function [u, history, X, T, xs, dump] = helper(method, h, maxT, CFL, f, aFunc)
    xs = -0.1:h:0.7;
    tau = h*CFL;
    tSpan = 0:tau:maxT;
    u0 = floor((1-sign(xs))/2);

    tic
    [u, history, dump] = method(tSpan, xs, u0, f, aFunc);
    toc
    [X, T] = meshgrid(xs, tSpan);
end