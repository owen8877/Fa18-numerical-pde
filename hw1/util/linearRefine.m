function matRef = linearRefine(mat)
    N = size(mat, 1);
    matRef = zeros(2*N-1);
    matRef(1:2:end, 1:2:end) = mat;
    matRef(2:2:end, 1:2:end) = (matRef(1:2:end-2, 1:2:end)+matRef(3:2:end, 1:2:end)) / 2;
    matRef(1:2:end, 2:2:end) = (matRef(1:2:end, 1:2:end-2)+matRef(1:2:end, 3:2:end)) / 2;
    matRef(2:2:end, 2:2:end) = (matRef(1:2:end-2, 1:2:end-2)+matRef(3:2:end, 3:2:end)) / 2;
end

