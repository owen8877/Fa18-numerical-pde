clear; %clc
addpath ../lib
strictMode = false;

%% Settings
% C-0 functions
Nxs = [8 16 32 64 128];
Nts = [4 8 16 32 64 128 256 512];
Nt0s = [100 200 400 800 1600 3200 6400 12800];
thetas = {'best-match', 0.5, 1};

maxT = 0.2;

errs = zeros(numel(Nts), numel(Nxs), numel(thetas)+1);

for i = 1:numel(Nts)
    Nt = Nts(i);

    for j = 1:numel(Nxs)
        Nx = Nxs(j);
    	xs = linspace(0, 1, Nx+1);
        u0 = max(xs, 1-xs)-1; % The test initial condition
        
        [X, T] = meshgrid(xs, linspace(0, maxT, Nt+1));
        Nsample = Nx/2;
        as = coeffGen(Nsample, u0, xs);
        uTruth = solutionGen(Nsample, as, X, T);

        for l = 1:numel(thetas)
            theta = thetas{l};
            
            [u, uHistory] = main(u0, Nx, Nt, maxT, theta, struct('strict', strictMode));
            
            errs(i, j, l) = max(max(abs(uHistory'-uTruth)));
            
            figure(1); clf
            subplot(2, 2, 1); mesh(X, T, uHistory'); xlabel('x'); ylabel('t'); title('numerical')
            subplot(2, 2, 3); mesh(X, T, uTruth); xlabel('x'); ylabel('t'); title('truth')
            subplot(2, 2, [2 4]); mesh(X, T, uHistory'-uTruth); xlabel('x'); ylabel('t');
            title(['error\newlineh=' num2str(1/Nx) ', \tau=' num2str(maxT/Nt) ', \theta=' num2str(theta)])
            return
        end
    end
end

for i = 1:numel(Nt0s)
    Nt = Nt0s(i);

    for j = 1:numel(Nxs)
        Nx = Nxs(j);
    	xs = linspace(0, 1, Nx+1);
        u0 = xs.^2 .* (1-xs); % The test initial condition
        
        [X, T] = meshgrid(xs, linspace(0, maxT, Nt+1));
        Nsample = Nx;
        as = coeffGen(Nsample, u0, xs);
        uTruth = solutionGen(Nsample, as, X, T);
        
        [u, uHistory] = main(u0, Nx, Nt, maxT, 0, struct('strict', strictMode));
            
        errs(i, j, end) = max(max(abs(uHistory'-uTruth)));
    end
end

plotHelper(Nxs, Nts, maxT, errs, thetas, 5, 1)
plotHelper(Nxs, Nts, maxT, errs, thetas, 6, 2)
plotHelper(Nxs, Nts, maxT, errs, thetas, 7, 3)
nerrs = errs;
nerrs(nerrs>10) = 0;
plotHelper(Nxs, Nt0s, maxT, nerrs, {0, 0, 0, 0}, 8, 4)

function plotHelper(Nxs, Nts, maxT, errs, thetas, figNum, thetaSlice)
    [NX, NT] = meshgrid(Nxs, Nts);
    figure(figNum); clf;
    subplot(2, 2, [2 4])
    mesh(1./NX, maxT./NT, errs(:, :, thetaSlice));
    title(['Error plot in log scale. theta=' num2str(thetas{thetaSlice})])
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'ZScale', 'log')
    xlabel('h'); ylabel('\tau'); zlabel('L_\infty error')
    
    subplot(2, 2, 1)
    mesh(1./NX, maxT./NT, errs(:, :, thetaSlice));
    title('In relation to \tau')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'ZScale', 'log')
    xlabel('h'); ylabel('\tau'); zlabel('L_\infty error')
    view([-90 0])
    
    subplot(2, 2, 3)
    mesh(1./NX, maxT./NT, errs(:, :, thetaSlice));
    title('In relation to h')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'ZScale', 'log')
    xlabel('h'); ylabel('\tau'); zlabel('L_\infty error')
    view([0 0])
    set(gcf, 'Position', [200, 200, 660, 400])
end