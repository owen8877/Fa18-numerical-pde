function l = laplacianr_2d(x, dir)
    if dir >= 0
        left = [x(:, 2:end) x(:, end-1)];
        right = [x(:, 2) x(:, 1:end-1)];
    end
    if dir <= 0
        up = [x(2:end, :); x(end-1, :)];
        down = [x(2, :); x(1:end-1, :)];
    end
    switch dir
        case 0
            l = left+down+up+right - 4*x;
        case 1
            l = left+right - 2*x;
        case -1
            l = up+down - 2*x;
    end
end

