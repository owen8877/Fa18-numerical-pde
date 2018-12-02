clear; %clc
addpath ../lib
warning off

queries = { ...
    {0.5, 1, floor(5 * (2.^(0:0.5:5)))}, ...
    {0.75, 1, floor(5 * (2.^(0:0.5:5)))}, ...
    {1, 1, floor(5 * (2.^(0:0.5:5)))}, ...
};
queryN = numel(queries);

results = cell(queryN, 1);
for i = 1:queryN
    query = queries{i};
    [err, tinfty] = thetaHelper(query{1}, query{2}, query{3});
    s.err = err;
    s.tinfty = tinfty;
    results{i} = s;
end

figure(10); clf; hold on
for i = 1:queryN
    query = queries{i};
    s = results{i};
    plot(query{2}./query{3}, s.err, 's-', 'DisplayName', num2str(query{1}));
    hold on
    
    yy = interp1(query{2}./query{3}, s.err, s.tinfty, 'linear', 'extrap');
    plot(s.tinfty, yy, 'd')
    
    a = query{2}./query{3};
    b = s.err;
    p = polyfit(log(a), log(b), 1);
    plot(a, exp(log(a)*p(1)+p(2)), '--')
    fprintf('Converge order for theta=%.2f : %.2f\n.', query{1}, p(1));
end
title('\theta>=0.5 case, h=1/32')
f = get(gca, 'Children');
title(legend(f(3:3:end)), '\theta')
xlabel('\tau')
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

function err = maxErrHelper(Nx, Nt, maxT, theta, show)
    [~, uHistory] = main(Nx, Nt, maxT, theta, struct('boundary', true));
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

function [errs, tinfty] = thetaHelper(theta, maxT, Nts)
    aMax = pi^(-2);
    Nx = 32; h = 1/Nx;
    [tinfty, ~] = criticalTau(aMax, h, theta);
    fprintf('linfty: %.3e\n', tinfty);

    errs = 0*Nts;
    for i = 1:numel(Nts)
        Nt = Nts(i);
        errs(i) = maxErrHelper(Nx, Nt, maxT, theta, i == floor(numel(Nts)/2));
    end
end