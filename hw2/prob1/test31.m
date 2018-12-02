clear; %clc
addpath ../lib
warning off

Nx = 32; maxTs = [0.1, 1, 5, 10, 50]; theta = 0.5;

Nts = maxTs ./ 0.1;
Nerrs = zeros(numel(Nts), 1);
Herrs = zeros(numel(Nts), 1);

for i = 1:numel(Nts)
    Nt = Nts(i);
    err = maxErrHelper(Nx, Nt, maxTs(i), theta, false, true);
    Nerrs(i) = err;
    err = maxErrHelper(Nx, Nt, maxTs(i), theta, false, false);
    Herrs(i) = err;
end

figure(9); clf; hold on
plot(maxTs, Nerrs, 's-', 'DisplayName', '0 order');
plot(maxTs, Herrs, 's-', 'DisplayName', '1 order');

title(['Comparison between two ways addressing boudary conditions, h=' num2str(1/Nx) ', \theta=' num2str(theta)])
legend()
xlabel('T_{max}')
ylabel('L_\infty Error')
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
grid on

function [tinfty, t2] = criticalTau(aMax, h, theta)
    if theta < 0.5
        t2 = 0.5 * h^2 / (1-2*theta) / aMax;
    else
        t2 = Inf;
    end
    tinfty = 0.5 * h^2 / (1-theta) / aMax;
end

function err = maxErrHelper(Nx, Nt, maxT, theta, show, naivety)
    [~, uHistory] = main(Nx, Nt, maxT, theta, struct('boundary', naivety));
    xs = linspace(0, 1, Nx+1);
    [X, T] = meshgrid(xs, linspace(0, maxT, Nt+1));
    uTruth = sin(pi*X).*(1+T).*exp(-T);
    err = max(max(abs(uHistory'-uTruth)));
    
    if show
        figure()
        xs = linspace(0, 1, Nx+1);
        [X, T] = meshgrid(xs, linspace(0, maxT, Nt+1));
        uTruth = sin(pi*X).*(1+T).*exp(-T);
        subplot(2, 2, 1); mesh(X, T, uHistory'); xlabel('x'); ylabel('t'); title('numerical')
        subplot(2, 2, 3); mesh(X, T, uTruth); xlabel('x'); ylabel('t'); title('truth')
        subplot(2, 2, [2 4]); mesh(X, T, uHistory'-uTruth); xlabel('x'); ylabel('t');
        title(['error\newline\tau=' num2str(maxT/Nt), ', h=' num2str(1/Nx), ', \theta=' num2str(theta)])
    end
end