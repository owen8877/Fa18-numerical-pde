function [u] = main(Nx, Nt, maxT, u0, theta, scheme, strictMode)
    %% Solver
    tau = maxT / Nt;
    h = 1 / Nx;

    Mu = 2*tau/h^2;
    if strcmp(theta, 'best-match')
        theta = 0.5 - 1/(12*Mu);
    end
    switch scheme
        case 'adi'
            theta = 0.5;
        case 'explicit'
            theta = 0;
        case 'implicit'
            theta = 1;
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

    for i = 1:Nt
        t = i * tau;
    %     fprintf('Time advanced to %.4f.\n', t);

        lapu = laplacianr_2d(u, 0);
        mu = tau / h^2;
        duProp = mu * lapu;

        if strcmp(scheme, 'adi')
            [duhalf, ~] = CG(duProp/2, @(v) v - mu*(theta*laplacianr_2d(v, 1)), ...
                duProp/2, struct());
            u = u + duhalf;

            lapu = laplacianr_2d(u, 0);
            duProp = mu * lapu;
            [du, ~] = CG(duProp/2, @(v) v - mu*(theta*laplacianr_2d(v, -1)), ...
                duProp/2, struct());
            u = u + du;
        else
            if theta ~= 0
                [du, ~] = CG(duProp, @(v) v - mu*(theta*laplacianr_2d(v, 0)), ...
                    duProp, struct());
                u = u + du;
            else
                u = u + duProp;
            end
        end
    end
end