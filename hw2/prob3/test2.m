clear; %clc
addpath ../lib
strictMode = false;

%% Settings
Nxs = [8 16 32 64 128];
maxT = 0.1;

theta = -1;
schemes = {'implicit', 'adi'};
Ntss = {20 * ones(1, 5)' * [1 4 16 32], ...
        10 * ones(1, 5)' * [1 4 16 32]};

errs = zeros(numel(Nxs), size(Ntss{1}, 2), numel(schemes));

for i = 1:numel(Nxs)
    Nx = Nxs(i);
    Ny = Nx;
    xs = linspace(0, 1, Nx+1); ys = xs;
    [XM, YM] = meshgrid(xs, ys);
    u0 = cos(pi*XM) .* cos(2*pi*YM);
    for l = 1:numel(schemes)
        scheme = schemes{l};
        Nts = Ntss{l};
        for j = 1:size(Nts, 2)
            Nt = Nts(i, j);
            uTruth = cos(pi*XM) .* cos(2*pi*YM) * exp(-5*pi^2*maxT);
            [u] = main(Nx, Nt, maxT, u0, theta, scheme, strictMode);
            
            errs(i, j, l) = max(max(abs(u-uTruth)));
            
            figure(1); clf
            subplot(2, 2, 1); mesh(XM, YM, u); xlabel('x'); ylabel('y'); title('numerical')
            subplot(2, 2, 3); mesh(XM, YM, uTruth); xlabel('x'); ylabel('y'); title('truth')
            subplot(2, 2, [2 4]); mesh(XM, YM, u-uTruth); xlabel('x'); ylabel('y');
            title(['error\newlineh=' num2str(1/Nx) ', \tau=' num2str(maxT/Nt) ', scheme=' scheme])
            % return
        end
    end
end

plotHelper(Nxs, Ntss{1}(1, :), maxT, errs, schemes, 5, 1);
plotHelper(Nxs, Ntss{2}(1, :), maxT, errs, schemes, 6, 2);

function plotHelper(Nxs, Nts, maxT, errs, schemes, figNum, thetaSlice)
    [NX, NT] = meshgrid(Nxs, Nts);
    figure(figNum); clf;
    subplot(2, 2, [2 4])
    mesh(1./NX, maxT./NT, errs(:, :, thetaSlice)');
    title(['Error plot in log scale. scheme=' schemes{thetaSlice}])
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'ZScale', 'log')
    xlabel('h'); ylabel('\tau'); zlabel('L_\infty error')
    
    subplot(2, 2, 1)
    mesh(1./NX, maxT./NT, errs(:, :, thetaSlice)');
    title('In relation to \tau')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'ZScale', 'log')
    xlabel('h'); ylabel('\tau'); zlabel('L_\infty error')
    view([-90 0])
    
    subplot(2, 2, 3)
    mesh(1./NX, maxT./NT, errs(:, :, thetaSlice)');
    title('In relation to h')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    set(gca, 'ZScale', 'log')
    xlabel('h'); ylabel('\tau'); zlabel('L_\infty error')
    view([0 0])
    
    set(gcf, 'Position', [200, 200, 800, 400])
end
