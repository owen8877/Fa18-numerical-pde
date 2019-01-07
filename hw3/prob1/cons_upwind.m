function [u, history, dump] = cons_upwind(tSpan, xs, u0, f, aFunc)
    u = u0;
    h = xs(2) - xs(1);
    history = zeros(numel(tSpan), numel(xs));
    
    history(1, :) = u;
    dump = [];
    for i = 1:numel(tSpan)-1
        dt = tSpan(i+1) - tSpan(i);
        
        fval = f(u);
        a = (fval(1:end-1)-fval(2:end))./(u(1:end-1)-u(2:end));
        mask = (u(1:end-1)-u(2:end)) < 1e-16;
        a(mask) = 0;
        
        s2 = sign(a(2:end)); s1 = sign(a(1:end-1));
%         u(2:end-1) = u(2:end-1) - dt/2/h * ( ...
%             (1+s2) .* fval(2:end-1) ...
%           + (1-s2) .* fval(3:end) ...
%           - (1+s1) .* fval(1:end-2) ...
%           - (1-s1) .* fval(2:end-1) ...
%         );
        u(2:end-1) = u(2:end-1) - dt/2/h * ( ...
            (s1+s2) .* fval(2:end-1) ...
          + (1-s2) .* fval(3:end) ...
          - (1+s1) .* fval(1:end-2) ...
        );
        history(i+1, :) = u;
        
        if mod(5*tSpan(i), 1) < 1e-6 && i > 1
            dump = [dump; u];
        end
    end
    dump = [dump; u];
end

