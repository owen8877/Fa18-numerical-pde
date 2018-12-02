function [u, uHistory] = main(Nx, Nt, maxT, theta, options)
    addpath ../lib
    strictMode = default(options, 'strict', false);
    naivety = default(options, 'boundary', false);
    fprintf('Naivety: %d\n', naivety);

    aFunc = @(xs, t) (1/((1+t)*pi^2)) .* ones(size(xs));
    aMax = pi^(-2);
    alphaFunc = @(t) 1+t^2;
    fFunc = @(xs, t) (-(t-1)*exp(-t)) .* sin(pi*xs);
    gFunc = @(t) pi * (1+t) * exp(-t);

%     Nx = 32;
%     Nt = 20*5;

%     maxT = 1*5;
%     theta = 0.95;
%     strictMode = false;

    tau = maxT / Nt;
    h = 1 / Nx;

    condNum = 2*aMax*tau/h^2*(1-2*theta);
    if theta < 0.5
        if condNum > 1
            fprintf('L2 Condnum %.2f > 1, violated!\n', condNum)
            if strictMode
                return
            end
        else
            fprintf('L2 Condnum %.2f <= 1.\n', condNum)
        end
    else
        fprintf('L2 Condnum check skipped (theta >= 0.5).\n')
    end

    maxCond = 2*aMax*tau/h^2*(1-theta);
    if maxCond > 1
        fprintf('%.2f; MV thm cond (L infty) is NOT satisfied!\n', maxCond)
    else
        fprintf('%.2f; MV thm cond (L infty) is satisfied!\n', maxCond)
    end

    xs = linspace(0, 1, Nx+1);
    if theta == 0
        a_handle = @(xs, t1, t2) aFunc(xs, t1);
    elseif theta == 1
        a_handle = @(xs, t1, t2) aFunc(xs, t2);
    else
    %     a_handle = @(xs, t1, t2) (aFunc(xs, t1) + aFunc(xs, t2)) / 2;
        a_handle = @(xs, t1, t2) aFunc(xs, (t1+t2)/2);
    end

    u0 = sin(pi*xs);
    u = u0;
    uHistory = zeros(Nx+1, Nt+1);
    uHistory(:, 1) = u;

    for i = 1:Nt
        t = i * tau;
%         fprintf('Time advanced to %.4f.\n', t);

        lapu = laplacian(u);
        as = a_handle(xs, t-tau, t);
        fs = fFunc(xs, t)*tau;
        mu = as * tau / h^2;
        duProp = mu .* lapu + fs;

        alpha1 = alphaFunc(t);
        g1 = gFunc(t);
        if naivety
            duProp(1) = g1 - (u(2)-u(1))/h + alpha1*u(1);
        else
            duProp(1) = g1 - (u(2)-u(1))/h + alpha1*(u(1)+u(2))/2;
        end

        musi = si(mu);
    %     [dusi, ~] = CG(si(duProp), @(v) helper(v, musi, theta, h, alpha1), ...
    %         si(duProp), struct());
        siduPropt = si(duProp)';
        [dusi, ~] = minres(@(v) helper(v', musi, theta, h, alpha1, naivety)', ...
            siduPropt, 1e-4, 1e4, [], [], siduPropt);
        u(1:end-1) = u(1:end-1) + dusi';

        uHistory(:, i+1) = u;

    %     muii = ii(mu);
    %     [duii, ~] = CG(ii(duProp), @(v) v - muii.*(theta*laplacian(v)), ...
    %         ii(duProp), struct());
    %     u(2:end-1) = u(2:end-1) + duii;
    %     
    %     uHistory(:, i+1) = u;
    end
end

function out = si(x)
    out = x(1:end-1);
end

function b = helper(v, musi, theta, h, alpha, naivety)
    b = v-musi.*(theta*laplacian(v));
    if naivety
        b(1) = (v(2)-v(1))/h - alpha*v(1);
    else
        b(1) = (v(2)-v(1))/h - alpha*(v(1)+v(2))/2;
    end  
end