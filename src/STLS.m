function x = STLS(b, A, l0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sequentially thresholded least-squares     %
% regression to approximate the solution to: %
% argmin|Ax-b|^2_2 + |l0 x|_0                %
%    x                                       %
%                                            %
% l0 : sparsity promoting parameter          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_norms = max(abs(A));
A = A./A_norms;
n = size(b,2);

x = pinv(A)*b;
for j=1:10
    for i=1:n
        small_inds = abs(x(:,i))<l0*max(abs(x(:,i)));
        x(small_inds,i) = 0;
        big_inds = ~small_inds;
        x(big_inds,i) = pinv(A(:,big_inds))*b(:,i);
    end
end

x = x./(A_norms');

end