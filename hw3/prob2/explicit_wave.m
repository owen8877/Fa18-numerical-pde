function [u, history] = explicit_wave(ts, xs, u0, v0, a, alpha0F, beta0F, alpha1F, beta1F, g0F, g1F)
    u = u0(xs);
    h = xs(2) - xs(1);
    history = zeros(numel(ts), numel(xs));
    
    history(1, :) = u;
    for i = 1:numel(ts)-1
        dt = ts(i+1) - ts(i);
        t = ts(i);
        nu = a*dt/h;
        
        alpha0 = alpha0F(t); alpha1 = alpha1F(t); g0 = g0F(t);
        beta0 = beta0F(t); beta1 = beta1F(t); g1 = g1F(t);
        
        ubak = u;
        if i == 1
            v0v = v0(xs);
            u(2:end-1) = (nu^2/2) * (u(3:end)+u(1:end-2)) + (1-nu^2)*u(2:end-1) + dt*v0v(2:end-1);
            if abs(beta0) < eps
                u(1) = g0 / alpha0;
            else
                u(1) = nu^2*(ubak(2) + h/beta0*(g0-alpha0*ubak(1))) + (1-nu^2)*u(1) + dt*v0v(1);
            end
            if abs(beta1) < eps
                u(end) = g1 / alpha1;
            else
                u(end) = nu^2*(ubak(end-1) + h/beta1*(g1-alpha1*ubak(end))) + (1-nu^2)*u(end) + dt*v0v(end);
            end
        else
            rhs = (u(3:end)+u(1:end-2)-2*u(2:end-1))*(a*dt/h)^2;
            u(2:end-1) = 2*u(2:end-1) - uold(2:end-1) + rhs;

            if abs(beta0) < eps
                u(1) = g0 / alpha0;
            else
                rhs0 = (ubak(2)-(1+alpha0*h/beta0)*ubak(1)+g0*h/beta0)*(2*(a*dt/h)^2);
                u(1) = 2*u(1) - uold(1) + rhs0;
            end

            if abs(beta1) < eps
                u(end) = g1 / alpha1;
            else
                rhs1 = (ubak(end-1)-(1+alpha1*h/beta1)*ubak(end)+g1*h/beta1)*(2*(a*dt/h)^2);
                u(end) = 2*u(end) - uold(end) + rhs1;
            end
        end
        uold = ubak;
        history(i+1, :) = u;
    end
end

