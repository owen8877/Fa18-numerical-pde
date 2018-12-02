function [x, history] = CG(x0, A, b, options)
    torlerance = default(options, 'torlerance', 1e-5);
    maxItr = default(options, 'maxItr', 1e4);
    mute = default(options, 'mute', false);
    
    x = x0;
    r0 = b - A(x);
    r = r0;
    z = r;
    p = z;
    
    itr = 0;
    history = [];
    r0Norm = matNorm(r0);
    while true
        rTz = matProd(r, z);
        Ap = A(p);
        alpha = rTz / matProd(p, Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        
%         history = [history matNorm(r)/r0Norm];
        if matNorm(r)/r0Norm < torlerance || itr > maxItr
            break
        end
        
        z = r;
        beta = matProd(r, z) / rTz;
        p = z + beta * p;
        
        itr = itr + 1;
        if ~mute && mod(itr, 10000) == 0
            fprintf('Iteration %d / %.3e.\n', itr, norm(r)/norm(r0));
        end
    end
end

