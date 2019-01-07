function set = system_set_2()
    set.u_truth = @(X, Y) exp(X.^3+Y.^2);
    set.u0 = @(xs, ys) exp(ys.^2);
    set.f = @(X, Y) -exp(X.^3+Y.^2) .* (6*X + 9*X.^4 + 2 + 4*Y.^2);
    set.beta = 1;
    set.g_x0 = @(xs, ys, u_truth) set.beta * u_truth(1, :);
    set.g_1y = @(xs, ys, u_truth) 3*exp(1+ys.^2) + set.beta * u_truth(:, end);
    set.g_x1 = @(xs, ys, u_truth) 2*exp(xs.^3+1) + set.beta * u_truth(end, :);
end

