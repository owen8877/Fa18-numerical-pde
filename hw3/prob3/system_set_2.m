function set = system_set_2()
    set.u_truth = @(X, Y) ((1-X).*(1-Y)).^2.*Y;
    set.u0 = @(xs, ys) ys.*(1-ys).^2;
    set.f = @(X, Y) -(6*Y-4).*((1-X).^2)-2.*(Y.*(1-Y).^2);
    set.beta = 10;
    set.g_x0 = @(xs, ys, u_truth) - (1-xs).^2 + set.beta * u_truth(1, :);
    set.g_1y = @(xs, ys, u_truth) 0*ys + set.beta * u_truth(:, end);
    set.g_x1 = @(xs, ys, u_truth) 0*xs + set.beta * u_truth(end, :);
end

