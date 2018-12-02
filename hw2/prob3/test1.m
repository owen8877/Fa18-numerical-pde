clear; %clc
addpath ../lib
strictMode = false;

%% Settings
Nxs = [8 16 24 32];
maxT = 0.1;

theta = -1;
scheme = 'explicit';
Nts = [64 128 256 512]' * [1 2 4 8 16];

errs = zeros(numel(Nxs), size(Nts, 2));

for i = 1:numel(Nxs)
    Nx = Nxs(i);
    Ny = Nx;
    xs = linspace(0, 1, Nx+1); ys = xs;
    [XM, YM] = meshgrid(xs, ys);
    u0 = cos(pi*XM) .* cos(2*pi*YM);
    for j = 1:size(Nts, 2)
        Nt = Nts(i, j);
        uTruth = cos(pi*XM) .* cos(2*pi*YM) * exp(-5*pi^2*maxT);
        [u] = main(Nx, Nt, maxT, u0, theta, scheme, strictMode);

        errs(i, j) = max(max(abs(u-uTruth)));

        figure(1); clf
        subplot(2, 2, 1); mesh(XM, YM, u); xlabel('x'); ylabel('y'); title('numerical')
        subplot(2, 2, 3); mesh(XM, YM, uTruth); xlabel('x'); ylabel('y'); title('truth')
        subplot(2, 2, [2 4]); mesh(XM, YM, u-uTruth); xlabel('x'); ylabel('y');
        title(['error\newlineh=' num2str(1/Nx) ', \tau=' num2str(maxT/Nt) ', scheme=' scheme])
        % return
    end
end

figure(1); clf

s1 = subplot(1, 2, 1); hold(s1, 'on')
for i = 1:numel(Nxs)
    plot(maxT./Nts(i, :), errs(i, :), 's-', 'DisplayName', num2str(1/Nxs(i)));
end
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('\tau')
ylabel('L_\infty error')
title('Convergence order of \tau')
grid on; title(legend(), 'h')

s2 = subplot(1, 2, 2); hold(s2, 'on')
plot(1./Nxs(1:3), errs(sub2ind(size(errs), 1:3, 3:-1:1)), 's-', 'DisplayName', num2str(maxT/64/4))
plot(1./Nxs(1:4), errs(sub2ind(size(errs), 1:4, 4:-1:1)), 's-', 'DisplayName', num2str(maxT/64/8))
plot(1./Nxs(1:4), errs(sub2ind(size(errs), 1:4, 5:-1:2)), 's-', 'DisplayName', num2str(maxT/64/16))
plot(1./Nxs(2:4), errs(sub2ind(size(errs), 2:4, 5:-1:3)), 's-', 'DisplayName', num2str(maxT/64/32))
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('h')
ylabel('L_\infty error')
title('Convergence order of h')
grid on; title(legend(), '\tau')
