clear; %clc
addpath ../lib

thetaHelper(0, 1, 150:10:250, 1)
thetaHelper(0.25, 1, 60:20:220, 2)
thetaHelper(0.4, 1, [30 40:20:200], 3)

function [tinfty, t2] = criticalTau(aMax, h, theta)
    if theta < 0.5
        t2 = 0.5 * h^2 / (1-2*theta) / aMax;
    else
        t2 = Inf;
    end
    tinfty = 0.5 * h^2 / (1-theta) / aMax;
end

function err = maxErrHelper(Nx, Nt, maxT, theta)
    [~, uHistory] = main(Nx, Nt, maxT, theta, struct());
    xs = linspace(0, 1, Nx+1);
    [X, T] = meshgrid(xs, linspace(0, maxT, Nt+1));
    uTruth = sin(pi*X).*(1+T).*exp(-T);
    err = max(max(abs(uHistory'-uTruth)));
    
%     xs = linspace(0, 1, Nx+1);
%     [X, T] = meshgrid(xs, linspace(0, maxT, Nt+1));
%     uTruth = sin(pi*X).*(1+T).*exp(-T);
%     subplot(2, 2, 1); mesh(X, T, uHistory'); xlabel('x'); ylabel('t'); title('numerical')
%     subplot(2, 2, 3); mesh(X, T, uTruth); xlabel('x'); ylabel('t'); title('truth')
%     subplot(2, 2, [2 4]); mesh(X, T, uHistory'-uTruth); xlabel('x'); ylabel('t'); title('error')
end

function thetaHelper(theta, maxT, Nts, figNum)
    aMax = pi^(-2);
    Nx = 32; h = 1/Nx;
    [tinfty, t2] = criticalTau(aMax, h, theta);
    fprintf('l2: %.3e, linfty: %.3e\n', t2, tinfty);

    errs = 0*Nts;
    for i = 1:numel(Nts)
        Nt = Nts(i);
        errs(i) = maxErrHelper(Nx, Nt, maxT, theta);
    end
    figure(figNum); clf; hold on
    plot(maxT./Nts, errs, 'DisplayName', 'numerical')
    set(gca, 'YScale', 'log')
    plot(t2 * [1 1], ylim, 'r:', 'DisplayName', 'l_2', 'LineWidth', 2)
    plot(tinfty * [1 1], ylim, 'b:', 'DisplayName', 'l_\infty', 'LineWidth', 2)
    title(['\theta=' num2str(theta) ' case, h=1/32'])
    legend()
    xlabel('\tau')
    ylabel('L_\infty Error')
    grid on
end