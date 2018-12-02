function [x, history] = CG(x0, A, b, options)
    if ~isfield(options, 'torlerance')
        options.torlerance = 1e-6;
    end
    x = x0;
    r0 = b - A*x;
    r = r0;
    z = r;
    p = z;
    
    itr = 0;
    history = [];
    while true
        rTz = r' * z;
        Ap = A * p;
        alpha = rTz / (p'*Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        
        history = [history norm(r)/norm(r0)];
        if norm(r) < options.torlerance
            break
        end
        
        z = r;
        beta = r'*z / rTz;
        p = z + beta * p;
        
        itr = itr + 1;
        if mod(itr, 100) == 0
            fprintf('Iteration %d / %.3e.\n', itr, norm(r)/norm(r0));
        end
    end
end

