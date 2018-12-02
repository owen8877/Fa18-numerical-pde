function sol = solutionGen(N, as, X, T)
    sol = 0*X;
    for i = 1:N
        sol = sol + as(i) * exp((-i^2*pi^2)*T) .* sin((i*pi).*X);
    end
end

