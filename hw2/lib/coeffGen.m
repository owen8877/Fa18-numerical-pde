function as = coeffGen(N, u0, xs)
    as = N*0;
    for i = 1:N
        as(i) = 2 * trapz(xs, u0.*sin((i*pi).*xs));
    end
end

