function [u, history, dump] = lax_wendroff(tSpan, xs, u0, f, aFunc)
    u = u0;
    h = xs(2) - xs(1);
    history = zeros(numel(tSpan), numel(xs));
    
    history(1, :) = u;
    dump = [];
    for i = 1:numel(tSpan)-1
        dt = tSpan(i+1) - tSpan(i);
        
        fval = f(u);
        a = aFunc((u(1:end-1)+u(2:end))/2);
        
        u(2:end-1) = u(2:end-1) - dt/2/h * (fval(3:end)-fval(1:end-2)) ...
            + dt^2/2/h^2 * (a(2:end).*(fval(3:end)-fval(2:end-1)) - a(1:end-1).*(fval(2:end-1)-fval(1:end-2)));
        history(i+1, :) = u;
        
        if mod(5*tSpan(i), 1) < 1e-6 && i > 1
            dump = [dump; u];
        end
    end
    dump = [dump; u];
end

