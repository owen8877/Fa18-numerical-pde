function l = laplacian(x)
    l = [x(2:end) 0] + [0 x(1:end-1)] - 2*x;
end

