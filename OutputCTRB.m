function M = OutputCTRB(A, B, C)

n = size(A, 1);
% [n1 , n2] = size(A);
% [m1 , m2] = size(B);
% [p1 , p2] = size(C);

M = C * B;
for k = 1 : n - 1
    M = cat(2, M, C*A^k*B);
end