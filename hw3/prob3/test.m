clear; %clc
rect = Kf_rect();
tri = Kf_tri();
set_1 = system_set_1();
set_2 = system_set_2();

Ns = [4 8 16 32 64 128 256];
options.maxItr = 1e3;
options.outputInt = 1e2;
options.torlerance = 1e-12;

%%
prob_sets = {set_1, set_2};
for l = 1:2
    prob_set = prob_sets{l};
    figure(l); clf; hold on
    mesh_Afs = {{'rect', rect, 'b'}, {'tri', tri, 'r'}};
    for j = 1:2
        mesh_Af = mesh_Afs{j}{2};
        l2errs = Ns * 0; linftyerrs = Ns * 0;
        for i = 1:numel(Ns)
            N = Ns(i);
            [u, u_truth] = helper(prob_set, mesh_Af, N, options);
            err = u - u_truth(:, 2:end);
            l2errs(i) = matNorm(err)/N;
            linftyerrs(i) = max(abs(err), [], 'all');
        end
        plot(1./Ns, l2errs, [mesh_Afs{j}{3} 's-'], 'DisplayName', [mesh_Afs{j}{1} ' l_2 norm'])
        plot(1./Ns, linftyerrs, [mesh_Afs{j}{3} 'd-'], 'DisplayName', [mesh_Afs{j}{1} ' l_\infty norm'])
    end
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    grid on
    xlim([1/Ns(end) 1/Ns(1)])
    xlabel('h')
    ylabel('error')
    legend('location', 'best')
    title(['System ' num2str(l)])
end
    
%%
function [u, u_truth] = helper(set, mesh_Af, N, options)
    h = 1/N;
    xs = (0:N)/N; ys = xs';
    [X, Y] = meshgrid(xs, ys);

    u_truth = set.u_truth(X, Y);
    u0 = set.u0(xs, ys);
    f = set.f(X, Y);
    beta = set.beta;
    g_x0 = set.g_x0(xs, ys, u_truth);
    g_1y = set.g_1y(xs, ys, u_truth);
    g_x1 = set.g_x1(xs, ys, u_truth);

    fv = mesh_Af.ffunc(f, u0, g_x0, g_1y, g_x1, mesh_Af.const, h, N);
    % Au_truth = mesh_Af.Afunc(u_truth(:, 2:end), mesh_Af.const, beta, h, N);

    u = CG(zeros(N+1, N), @(u) mesh_Af.Afunc(u, mesh_Af.const, beta, h, N), fv, options);
end