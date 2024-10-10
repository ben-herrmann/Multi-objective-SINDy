function k = my_cond(A)

[~,s,~] = svd(A,'econ','vector');
k = max(s)/min(s);

end