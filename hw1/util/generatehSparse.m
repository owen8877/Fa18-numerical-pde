function A = generatehSparse(N, ah)
    M = N*(N-1);
    O = 5*N^2 - 9*N + 2;
    ix = zeros(1, O);
    iy = zeros(1, O);
    v = zeros(1, O);
    
    ix(1:N-2) = 2:N-1;
    iy(1:N-2) = 1:N-2;
    v(1:N-2) = -0.5;
    
    ix((1:N-2)+N-2) = 1:N-2;
    iy((1:N-2)+N-2) = 2:N-1;
    v((1:N-2)+N-2) = -0.5;
    
    counter = 2*(N-2);
    for i = 1:N-1
        ix(counter+(1:N-2)) = i*(N-1)+1+(1:N-2);
        iy(counter+(1:N-2)) = i*(N-1)+(1:N-2);
        v(counter+(1:N-2)) = -1;
        counter = counter + N-2;
        
        ix(counter+(1:N-2)) = i*(N-1)+(1:N-2);
        iy(counter+(1:N-2)) = i*(N-1)+1+(1:N-2);
        v(counter+(1:N-2)) = -1;
        counter = counter + N-2;
    end
    
    ix(counter+(1:N-1)) = 1:N-1;
    iy(counter+(1:N-1)) = 1:N-1;
    v(counter+(1:N-1)) = 2-ah;
    counter = counter + N-1;
    
    ix(counter+(1:(N-1)^2)) = N:M;
    iy(counter+(1:(N-1)^2)) = N:M;
    v(counter+(1:(N-1)^2)) = 4;
    counter = counter + (N-1)^2;
    
    ix(counter+(1:(N-1)^2)) = N:M;
    iy(counter+(1:(N-1)^2)) = 1:(N-1)^2;
    v(counter+(1:(N-1)^2)) = -1;
    counter = counter + (N-1)^2;
    
    ix(counter+(1:(N-1)^2)) = 1:(N-1)^2;
    iy(counter+(1:(N-1)^2)) = N:M;
    v(counter+(1:(N-1)^2)) = -1;
    counter = counter + (N-1)^2;
    
    A = sparse(ix, iy, v);
end