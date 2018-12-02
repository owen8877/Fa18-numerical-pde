clear; %clc
addpath ../lib
strictMode = false;

%% Settings
% piecewise C-0 functions
Nxs = [8 16 32 64 128];
Nts = [4 8 16 32 64 128 256 512];
Nt0s = [100 200 400 800 1600 3200 6400 12800];
thetas = {'best-match', 0.5, 1};

maxT = 0.2;

errs = zeros(numel(Nts), numel(Nxs), numel(thetas)+1);
l2errs = zeros(numel(Nts), numel(Nxs), numel(thetas)+1);

for i = 1:numel(Nts)
    Nt = Nts(i);

    for j = 1:numel(Nxs)
        Nx = Nxs(j);
    	xs = linspace(0, 1, Nx+1);
        u0 = sign(xs-0.5) .* sin(pi*xs); % The test initial condition
        
        [X, T] = meshgrid(xs, linspace(0, maxT, Nt+1));
        Nsample = 1024;
%         as = coeffGen(Nsample, u0, xs);
        as = zeros(1, Nsample); ks = 1:Nsample/2;
        as(2:2:end) = (4*ks.*(-1).^ks)./(4*ks.^2-1)./pi * 2;
        uTruth = solutionGen(Nsample, as, X, T);

        for l = 1:numel(thetas)
            theta = thetas{l};
            
            [u, uHistory] = main(u0, Nx, Nt, maxT, theta, struct('strict', strictMode));
            
            errs(i, j, l) = max(max(abs(uHistory'-uTruth)));
            l2errs(i, j, l) = mean(mean((uHistory'-uTruth).^2))^0.5;
            
%             figure(1); clf
%             subplot(2, 2, 1); mesh(X, T, uHistory'); xlabel('x'); ylabel('t'); title('numerical')
%             subplot(2, 2, 3); mesh(X, T, uTruth); xlabel('x'); ylabel('t'); title('truth')
%             subplot(2, 2, [2 4]); mesh(X, T, uHistory'-uTruth); xlabel('x'); ylabel('t');
%             title(['error\newlineh=' num2str(1/Nx) ', \tau=' num2str(maxT/Nt) ', \theta=' num2str(theta)])
%             return
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
        l2errs(i, j, end) = mean(mean((uHistory'-uTruth).^2))^0.5;
    end
end

plotHelper(Nxs, Nts, maxT, l2errs, errs, thetas, 5, 1)
plotHelper(Nxs, Nts, maxT, l2errs, errs, thetas, 6, 2)
plotHelper(Nxs, Nts, maxT, l2errs, errs, thetas, 7, 3)
nerrs = errs;
nerrs(nerrs>10) = 0;
nl2errs = l2errs;
nl2errs(nl2errs>1) = 0;
plotHelper(Nxs, Nt0s, maxT, nl2errs, nerrs, {0, 0, 0, 0}, 8, 4)

function plotHelper(Nxs, Nts, maxT, l2errs, errs, thetas, figNum, thetaSlice)
    [NX, NT] = meshgrid(Nxs, Nts);
    figure(figNum); clf;
    subplot(2, 3, 1)
    mesh(1./NX, maxT./NT, errs(:, :, thetaSlice));
    title(['Error plot in log scale. theta=' num2str(thetas{thetaSlice})])
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'ZScale', 'log')
    xlabel('h'); ylabel('\tau'); zlabel('L_\infty error')
    
    subplot(2, 3, 2)
    mesh(1./NX, maxT./NT, errs(:, :, thetaSlice));
    title('In relation to \tau')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'ZScale', 'log')
    xlabel('h'); ylabel('\tau'); zlabel('L_\infty error')
    view([-90 0])
    
    subplot(2, 3, 3)
    mesh(1./NX, maxT./NT, errs(:, :, thetaSlice));
    title('In relation to h')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'ZScale', 'log')
    xlabel('h'); ylabel('\tau'); zlabel('L_\infty error')
    view([0 0])
    
    
    subplot(2, 3, 4)
    mesh(1./NX, maxT./NT, l2errs(:, :, thetaSlice));
    title(['L_2 error plot in log scale. theta=' num2str(thetas{thetaSlice})])
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'ZScale', 'log')
    xlabel('h'); ylabel('\tau'); zlabel('L_2 error')
    
    subplot(2, 3, 5)
    mesh(1./NX, maxT./NT, l2errs(:, :, thetaSlice));
    title('In relation to \tau')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'ZScale', 'log')
    xlabel('h'); ylabel('\tau'); zlabel('L_2 error')
    view([-90 0])
    
    subplot(2, 3, 6)
    mesh(1./NX, maxT./NT, l2errs(:, :, thetaSlice));
    title('In relation to h')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'ZScale', 'log')
    xlabel('h'); ylabel('\tau'); zlabel('L_2 error')
    view([0 0])
    
    
    set(gcf, 'Position', [200, 200, 800, 400])
end