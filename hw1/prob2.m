function [s, u, u3x, u3y, u4x, u4y] = prob2()
    s = @solver;
    u = @uF;
    u3x = @u3xF;
    u3y = @u3yF;
    u4x = @u4xF;
    u4y = @u4yF;
end

function u = uF(x, y)
	u = exp(x .* y);
end

function u3x = u3xF(x, y)
	u3x = y.^3 .* exp(x .* y);
end

function u3y = u3yF(x, y)
    u3y = x.^3 .* exp(x .* y);
end

function u4x = u4xF(x, y)
	u4x = y.^4 .* exp(x .* y);
end

function u4y = u4yF(x, y)
	u4y = x.^4 .* exp(x .* y);
end

function [uSolMat, uTruth] = solver(N, initial)
    uFunc = @uF;
    fFunc = @(x, y) (x.^2+y.^2) .* exp(x .* y);
    gFunc = @(x, y) y .* exp(x .* y);
    d1Func = @(x, y) x*0 + 1;
    d2Func = @(x, y) exp(y);
    d3Func = @(x, y) exp(x);

    h = 1/N;

    xVec = (0:N-1)/N;
    yVec = (1:N-1)/N;
    [interX, interY] = meshgrid(xVec, yVec);
    fMat = fFunc(interX, interY);
    gVec = gFunc(0, yVec)';
    d1Vec = d1Func(xVec, 0)';
    d2Vec = d2Func(1, yVec)';
    d3Vec = d3Func(xVec, 0)';

    rhs = reshape(-fMat, (N-1)*N, 1) * h^2;
    rhs(1:N-1) = rhs(1:N-1) - gVec * (2*h);
    rhs(1:N-1:end) = rhs(1:N-1:end) + d1Vec;
    rhs(N-1:N-1:end) = rhs(N-1:N-1:end) + d3Vec;
    rhs((N-1)^2+1:end) = rhs((N-1)^2+1:end) + d2Vec;

    lhsE = generatehSparse(N, 0);
    rhsM = rhs; rhsM(1:N-1) = rhsM(1:N-1) / 2;
    if size(initial, 1) == 0
        initial = zeros(size(rhs, 1), 1);
    else
        initial = reshape(initial(2:end-1, 1:end-1), size(rhs, 1), 1);
    end
    uSolVec = CG(zeros(size(rhsM, 1), 1), lhsE, rhsM, struct('torlerance', 1e-8));
    uSolMat = zeros(N+1);
    uSolMat(2:N, 1:N) = reshape(uSolVec, N-1, N);
    uSolMat(1, 1:N) = d1Vec';
    uSolMat(N+1, 1:N) = d3Vec';
    uSolMat(:, N+1) = [d2Func(1, 0); d2Vec; d2Func(1, 1)];

    xDisVec = (0:N)/N;
    yDisVec = (0:N)/N;
    [disX, disY] = meshgrid(xDisVec, yDisVec);
    uTruth = uFunc(disX, disY);
    mesh(disX, disY, uSolMat-uTruth); title('Error to the ground truth')
end