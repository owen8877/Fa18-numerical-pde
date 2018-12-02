function l = laplacian_2d(x)
    left = [x(:, 2:end) x(:, end)];
    right = [x(:, 1) x(:, 1:end-1)];
    up = [x(2:end, :); x(end, :)];
    down = [x(1, :); x(1:end-1, :)];
    l = left+down+up+right - 4*x;
end

