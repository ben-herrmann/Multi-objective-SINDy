function [error,success, lambda] = performance(Xdot,Theta,Xi_true)

lambdas = logspace(-4,0,100)';
% lambdas = [0.8];
errors = zeros(length(lambdas),1);

for i=1:length(lambdas)
    Xi = STLS(Xdot,Theta,lambdas(i));
    errors(i) = norm(Xi-Xi_true,'fro')/norm(Xi_true,'fro');
end

[error,i] = min(errors);
lambda = lambdas(i);
Xi = STLS(Xdot,Theta,lambda);
success = isequal(~Xi, ~Xi_true);

end