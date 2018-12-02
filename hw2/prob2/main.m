function [u, uHistory] = main(u0, Nx, Nt, maxT, theta, options)
    addpath ../lib
    
    %% Solver
    tau = maxT / Nt;
    h = 1 / Nx;
    strictMode = default(options, 'strict', false);

    Mu = tau/h^2;
    if strcmp(theta, 'best-match')
        theta = 0.5 - 1/(12*Mu);
    end
    
    fprintf('h=%.3f\ttau=%.3f\ttheta=%.3f\n', h, tau, theta)
    fprintf('Mu = %.2f.\n', Mu);
    condNum = 2*Mu*(1-2*theta);
    if theta < 0.5
        if condNum > 1
            fprintf('Condnum %.2f > 1, violated!\n', condNum)
            if strictMode
                return
            end
        else
            fprintf('Condnum %.2f <= 1.\n', condNum)
        end
    else
        fprintf('Condnum check skipped (theta >= 0.5).\n')
    end

    maxCond = 2*Mu*(1-theta);
    if maxCond > 1
        fprintf('%.2f; MV thm cond is NOT satisfied!\n', maxCond)
    else
        fprintf('%.2f; MV thm cond is satisfied!\n', maxCond)
    end

    u = u0;
    uHistory = zeros(Nx+1, Nt+1);
    uHistory(:, 1) = u;

    for i = 1:Nt
        t = i * tau;
        % fprintf('Time advanced to %.4f.\n', t);

        lapu = laplacian(u);
        mu = tau / h^2;
        duProp = mu * lapu;

        duPropii = ii(duProp);
        if theta == 0
            duii = duPropii;
        else
            [duii, ~] = CG(duPropii, @(v) v - mu*(theta*laplacian(v)), ...
                duPropii, struct());
        end
        u(2:end-1) = u(2:end-1) + duii;
        uHistory(:, i+1) = u;
    end
end

function out = ii(x)
    out = x(2:end-1);
end