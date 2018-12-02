function [s, u, u3x, u3y, u4x, u4y] = prob1()
    qQuant = [2; -2];
    qPosition = [2, 1; 2, 2];
    alpha = 1;
    s = @(N, initialGuess) solver(qQuant, qPosition, alpha, N, initialGuess);
    u = @(x, y) uF(x, y, qQuant, qPosition);
    u3x = @(x, y) u3xF(x, y, qQuant, qPosition);
    u3y = @(x, y) u3yF(x, y, qQuant, qPosition);
    u4x = @(x, y) u4xF(x, y, qQuant, qPosition);
    u4y = @(x, y) u4yF(x, y, qQuant, qPosition);
end

function u = uF(x, y, qQuant, qPosition)
	u = ...
        qQuant(1).*log(radius(x-qPosition(1, 1), y-qPosition(1, 2))) + ...
        qQuant(2).*log(radius(x-qPosition(2, 1), y-qPosition(2, 2)));
end

function u3x = u3xF(x, y, qQuant, qPosition)
    dx1 = x-qPosition(1, 1); dx2 = x-qPosition(2, 1);
    dy1 = y-qPosition(1, 2); dy2 = y-qPosition(2, 2);
	u3x = ...
        qQuant(1).*(2*dx1.^3-6*dx1.*dy1.^2)./radius(dx1, dy1).^6 + ...
        qQuant(2).*(2*dx2.^3-6*dx2.*dy2.^2)./radius(dx2, dy2).^6;
end

function u3y = u3yF(x, y, qQuant, qPosition)
    dx1 = x-qPosition(1, 1); dx2 = x-qPosition(2, 1);
    dy1 = y-qPosition(1, 2); dy2 = y-qPosition(2, 2);
	u3y = ...
        qQuant(1).*(2*dy1.^3-6*dy1.*dx1.^2)./radius(dx1, dy1).^6 + ...
        qQuant(2).*(2*dy2.^3-6*dy2.*dx2.^2)./radius(dx2, dy2).^6;
end

function u4x = u4xF(x, y, qQuant, qPosition)
    dx1 = x-qPosition(1, 1); dx2 = x-qPosition(2, 1);
    dy1 = y-qPosition(1, 2); dy2 = y-qPosition(2, 2);
	u4x = ...
        6*qQuant(1).*(-dx1.^4+6*dx1.^2.*dy1.^2-dy1.^4)./radius(dx1, dy1).^8 + ...
        6*qQuant(2).*(-dx2.^4+6*dx2.^2.*dy2.^2-dy2.^4)./radius(dx2, dy2).^8;
end

function u4y = u4yF(x, y, qQuant, qPosition)
    dx1 = x-qPosition(1, 1); dx2 = x-qPosition(2, 1);
    dy1 = y-qPosition(1, 2); dy2 = y-qPosition(2, 2);
	u4y = ...
        6*qQuant(1).*(-dx1.^4+6*dx1.^2.*dy1.^2-dy1.^4)./radius(dx1, dy1).^8 + ...
        6*qQuant(2).*(-dx2.^4+6*dx2.^2.*dy2.^2-dy2.^4)./radius(dx2, dy2).^8;
end

function [uSolMat, uTruth] = solver(qQuant, qPosition, alpha, N, initial)
    uFunc = @(x, y) uF(x, y, qQuant, qPosition);
    fFunc = @(x, y) 0*x;
    d1Func = @(x) uFunc(x, 0*x);
    d2Func = @(y) uFunc(0*y+1, y);
    d3Func = @(x) uFunc(x, 0*x+1);
    uxFunc = @(x, y) ...
        qQuant(1).*(x-qPosition(1, 1))./radius(x-qPosition(1, 1), y-qPosition(1, 2)).^2 + ...
        qQuant(2).*(x-qPosition(2, 1))./radius(x-qPosition(2, 1), y-qPosition(2, 2)).^2;
    gFunc = @(y) alpha * uFunc(0*y, y) + uxFunc(0*y, y);

    h = 1/N;

    xVec = (0:N-1)/N;
    yVec = (1:N-1)/N;
    [interX, interY] = meshgrid(xVec, yVec);
    fMat = fFunc(interX, interY);
    d1Vec = d1Func(xVec)';
    d2Vec = d2Func(yVec)';
    d3Vec = d3Func(xVec)';
    gVec = gFunc(yVec)';

    rhs = reshape(-fMat, (N-1)*N, 1) * h^2;
    rhs(1:N-1:end) = rhs(1:N-1:end) + d1Vec;
    rhs(1:N-1) = rhs(1:N-1) - (2*h)*gVec;
    rhs(N-1:N-1:end) = rhs(N-1:N-1:end) + d3Vec;
    rhs((N-1)^2+1:end) = rhs((N-1)^2+1:end) + d2Vec;

    lhsE = generatehSparse(N, alpha*h);
    rhsM = rhs; rhsM(1:N-1) = rhsM(1:N-1) / 2;
    if size(initial, 1) == 0
        initial = zeros(size(rhs, 1), 1);
    else
        initial = reshape(initial(2:end-1, 1:end-1), size(rhs, 1), 1);
    end
    uSolVec = CG(initial, lhsE, rhsM, struct('torlerance', 1e-8));
    % uSolVec = lhsE \ rhsM;
    
    uSolMat = zeros(N+1);
    uSolMat(2:N, 1:N) = reshape(uSolVec, N-1, N);
    uSolMat(1, 1:N) = d1Vec';
    uSolMat(2:N, N+1) = d2Vec;
    uSolMat(N+1, 1:N) = d3Vec';
    uSolMat(1, N+1) = d2Func(0);
    uSolMat(N+1, N+1) = d3Func(1);

    xDisVec = (0:N)/N;
    yDisVec = (0:N)/N;
    [disX, disY] = meshgrid(xDisVec, yDisVec);
    uTruth = uFunc(disX, disY);
    figure(1); clf
    mesh(disX, disY, uSolMat-uTruth); title('Error to the ground truth')
    figure(2); clf
    subplot(1, 2, 1)
    mesh(disX, disY, uSolMat); title('Numerical solution')
    subplot(1, 2, 2)
    mesh(disX, disY, uFunc(disX, disY)); title('Ground truth')
end

function A = generateSparse(N)
    M = (N-1)^2;
    O = 5*N^2 - 14*N + 9;
    ix = zeros(1, O);
    iy = zeros(1, O);
    v = zeros(1, O);
    counter = 0;
    
    for i = 0:N-2
        ix(counter+(1:N-2)) = i*(N-1)+1+(1:N-2);
        iy(counter+(1:N-2)) = i*(N-1)+(1:N-2);
        v(counter+(1:N-2)) = -1;
        counter = counter + N-2;
        
        ix(counter+(1:N-2)) = i*(N-1)+(1:N-2);
        iy(counter+(1:N-2)) = i*(N-1)+1+(1:N-2);
        v(counter+(1:N-2)) = -1;
        counter = counter + N-2;
    end
    
    ix(counter+(1:M)) = 1:M;
    iy(counter+(1:M)) = 1:M;
    v(counter+(1:M)) = 4;
    counter = counter + M;
    
    ix(counter+(1:M-N+1)) = N:M;
    iy(counter+(1:M-N+1)) = 1:(N-1)*(N-2);
    v(counter+(1:M-N+1)) = -1;
    counter = counter + M-N+1;
    
    ix(counter+(1:M-N+1)) = 1:(N-1)*(N-2);
    iy(counter+(1:M-N+1)) = N:M;
    v(counter+(1:M-N+1)) = -1;
    
    A = sparse(ix, iy, v);
end