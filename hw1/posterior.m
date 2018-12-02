clear; %clc
addpath util

% ** Change the line below to choose the problem **
prob = @prob2;
[solver, ~, ~, ~, ~, ~] = prob();

N0 = 4; N = N0;
errorHistory = [];
uCoarseSol = [];
while true
    [uSol, ~] = solver(N, linearRefine(uCoarseSol));
    
    if ~isempty(uCoarseSol)
        % Examine if the posterior error has been small enough
        errorMat = abs(uSol(1:2:end, 1:2:end) - uCoarseSol) * 4/3;
        inftyErr = max(max(errorMat));
        fprintf('Solution at size %d has estimated error of %.4e.\n', N, inftyErr);

        errorHistory = [errorHistory inftyErr];
        if inftyErr < 1e-5
            break
        end
    end
    
    uCoarseSol = uSol;
    N = N*2;
end

figure; grid on
loglog(2.^(0:numel(errorHistory)-1)*N0, errorHistory)
xlabel('N'); ylabel('\epsilon'); title('Posterior estimated error')
fitCoeff = polyfit(log(2.^(0:numel(errorHistory)-1)*N0), log(errorHistory), 1);
fprintf('The convergence order is %.3f.\n', -fitCoeff(1))