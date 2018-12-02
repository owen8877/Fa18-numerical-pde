clear; %clc
addpath util

% ** Change the line below to choose the problem **
prob = @prob2;
[~, uF, u3xF, u3yF, u4xF, u4yF] = prob();

totalTor = 1e-5;
linAlTor = 1e-6;
C = 17/16;

alpha = 1;
estN = 20;
xEstVec = (0:estN)/estN;
yEstVec = (0:estN)/estN;
[estX, estY] = meshgrid(xEstVec, yEstVec);

% First, interior estimation
u4xi = u4xF(estX(2:end-1, 2:end-1), estY(2:end-1, 2:end-1));
u4yi = u4yF(estX(2:end-1, 2:end-1), estY(2:end-1, 2:end-1));
fprintf('u4xi: %.4e, u4yi: %.4e\n', maxAbs(u4xi), maxAbs(u4yi))

% Then, for Omega4 estimation
u3xb = u3xF(estX(:, 1), estY(:, 1));
u4yb = u4yF(estX(:, 1), estY(:, 1));
fprintf('u3xb: %.4e, u4yb: %.4e\n', maxAbs(u3xb), maxAbs(u4yb))

intCoeff = (maxAbs(u4xi) + maxAbs(u4yi))/12;
borCoeff = (maxAbs(u3xb) + maxAbs(u4yb)/8)/6;

h = sqrt((totalTor/C-linAlTor)/max([intCoeff, borCoeff]));
fprintf('Max grid length: %.4e, min grid number: %d\n', h, ceil(1/h))