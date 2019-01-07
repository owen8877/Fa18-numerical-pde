function set = system_set_1()
    set.u_truth = @(X, Y) sin(X*2*pi).*sin(Y*2*pi) + (1-X).*Y.*(1-Y);
    set.u0 = @(xs, ys) ys.*(1-ys);
    set.f = @(X, Y) (8*pi^2)*sin(X*2*pi).*sin(Y*2*pi) + 2*(1-X);
    set.beta = 0;
    set.g_x0 = @(xs, ys, u_truth) -(sin(xs*2*pi)*2*pi + (1-xs)) + set.beta * u_truth(1, :);
    set.g_1y = @(xs, ys, u_truth) 2*pi*sin(ys*2*pi)-ys.*(1-ys) + set.beta * u_truth(:, end);
    set.g_x1 = @(xs, ys, u_truth) sin(xs*2*pi)*2*pi - (1-xs) + set.beta * u_truth(end, :);
end

