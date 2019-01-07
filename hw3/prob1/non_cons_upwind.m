function [u, history, dump] = non_cons_upwind(tSpan, xs, u0, f, afunc)
    u = u0;
    h = xs(2) - xs(1);
    history = zeros(numel(tSpan), numel(xs));
    
    history(1, :) = u;
    dump = [];
    for i = 1:numel(tSpan)-1
        dt = tSpan(i+1) - tSpan(i);
        
        a = afunc(u);
        nu = a * dt / h;
        nu = nu(2:end-1);
        
        dup = u(3:end) - u(2:end-1);
        dun = u(2:end-1) - u(1:end-2);
        
        s = ceil((1+sign(a(2:end-1)))/2);
        u(2:end-1) = u(2:end-1) - nu .* ((1-s) .* dup + s .* dun);
        history(i+1, :) = u;
        
        if mod(5*tSpan(i), 1) < 1e-6 && i > 1
            dump = [dump; u];
        end
    end
    dump = [dump; u];
end

