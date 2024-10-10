%% Saddle-node hysteresys loop

clear variables; close all; clc
addpath('../src'); % Add the source files to the path
load('../data/data_SNHys.mat')

Xi_true = zeros(length(theta),1);
Xi_true(string(theta)=='u') = 1;
Xi_true(string(theta)=='x') = 1;
Xi_true(string(theta)=='x^3') = -1;

X_dot_tr = X_dot_tr(:,1);
X_dot_att = X_dot_att(:,1);
Theta_1tr = Theta_tr(X_tr(:,2)==mu_tr,:);
X_dot_1tr = X_dot_tr(X_tr(:,2)==mu_tr,:);

[error_1tr,success_1tr,l_1tr] = performance(X_dot_1tr,Theta_1tr,Xi_true);
[error_att,success_att,l_att] = performance(X_dot_att,Theta_att,Xi_true);
[error_both,success_both,l_both] = performance([X_dot_1tr; X_dot_att],[Theta_1tr; Theta_att],Xi_true);
[error_full,success_full,l_full] = performance([X_dot_tr; X_dot_att],[Theta_tr; Theta_att],Xi_true);

alphas = logspace(-20,20,11)';
ls_mo = zeros(length(alphas),1);
errors_mo = zeros(length(alphas),1);
for i=1:length(alphas)
    X_dot = [X_dot_1tr; sqrt(alphas(i))*X_dot_att];
    Theta = [Theta_1tr; sqrt(alphas(i))*Theta_att];
    [errors_mo(i),~,ls_mo(i)] = performance(X_dot,Theta,Xi_true);
end
[error_mo,i] = min(errors_mo);
alpha = alphas(i);
l_mo = ls_mo(i);
X_dot = [X_dot_1tr; sqrt(alpha)*X_dot_att];
Theta = [Theta_1tr; sqrt(alpha)*Theta_att];
Xi = STLS(X_dot,Theta,l_mo);
success_mo = isequal(~Xi, ~Xi_true);

error =[error_1tr, error_att, error_both, error_full, error_mo];
success = [success_1tr, success_att, success_both, success_full, success_mo];
lambda = [l_1tr, l_att, l_both, l_full, l_mo];
disp(error);
disp(success);
disp(lambda);
disp(alpha)

save('../data/performance_SNHys.mat',"Xi_true","error","success","lambda","alpha",'-v7.3');

%% Hopf normal form
clear variables; close all; clc
addpath('../src'); % Add the source files to the path
load('../data/data_Hopf.mat')

Xi_true = zeros(length(theta),2);
Xi_true(string(theta)=='x u',1) = 1;
Xi_true(string(theta)=='y',1) = -1;
Xi_true(string(theta)=='x^3',1) = -1;
Xi_true(string(theta)=='x y^2',1) = -1;
Xi_true(string(theta)=='x',2) = 1;
Xi_true(string(theta)=='y u',2) = 1;
Xi_true(string(theta)=='y^3',2) = -1;
Xi_true(string(theta)=='x^2 y',2) = -1;

X_dot_tr = X_dot_tr(:,1:2);
X_dot_att = X_dot_att(:,1:2);
Theta_1tr = Theta_tr(X_tr(:,3)==mu_tr,:);
X_dot_1tr = X_dot_tr(X_tr(:,3)==mu_tr,:);

[error_1tr,success_1tr,l_1tr] = performance(X_dot_1tr,Theta_1tr,Xi_true);
[error_att,success_att,l_att] = performance(X_dot_att,Theta_att,Xi_true);
[error_both,success_both,l_both] = performance([X_dot_1tr; X_dot_att],[Theta_1tr; Theta_att],Xi_true);
% [error_full,success_full,l_full] = performance([X_dot_tr; X_dot_att],[Theta_tr; Theta_att],Xi_true);
error_full = 0; success_full = 1; l_full = 0;

alphas = logspace(-20,20,11)';
ls_mo = zeros(length(alphas),1);
errors_mo = zeros(length(alphas),1);
for i=1:length(alphas)
    X_dot = [X_dot_1tr; sqrt(alphas(i))*X_dot_att];
    Theta = [Theta_1tr; sqrt(alphas(i))*Theta_att];
    [errors_mo(i),~,ls_mo(i)] = performance(X_dot,Theta,Xi_true);
end
[error_mo,i] = min(errors_mo);
alpha = alphas(i);
l_mo = ls_mo(i);
X_dot = [X_dot_1tr; sqrt(alpha)*X_dot_att];
Theta = [Theta_1tr; sqrt(alpha)*Theta_att];
Xi = STLS(X_dot,Theta,l_mo);
success_mo = isequal(~Xi, ~Xi_true);

error =[error_1tr, error_att, error_both, error_full, error_mo];
success = [success_1tr, success_att, success_both, success_full, success_mo];
lambda = [l_1tr, l_att, l_both, l_full, l_mo];
disp(error);
disp(success);
disp(lambda);
disp(alpha);

% save('../data/performance_Hopf.mat',"Xi_true","error","success","lambda","alpha",'-v7.3');

%% Coupled Stuart–Landau oscillators
clear variables; close all; clc
addpath('../src'); % Add the source files to the path
load('../Final_Results/data_SL.mat')

R1 = 1; R2 = 0.2; w1 = 1/pi; w2 =1;
Xi_true = zeros(length(theta),4);
Xi_true(string(theta)=='y1',1) = -w1;
Xi_true(string(theta)=='x1',1) = R1^2;
Xi_true(string(theta)=='x1^3',1) = -1;
Xi_true(string(theta)=='x1 y1^2',1) = -1;
Xi_true(string(theta)=='x1 x2 u',1) = 1;
Xi_true(string(theta)=='y1 y2 u',1) = 1;
Xi_true(string(theta)=='x1',2) = w1;
Xi_true(string(theta)=='y1',2) = R1^2;
Xi_true(string(theta)=='y1^3',2) = -1;
Xi_true(string(theta)=='x1^2 y1',2) = -1;
Xi_true(string(theta)=='x1 y2 u',2) = 1;
Xi_true(string(theta)=='y1 x2 u',2) = -1;
Xi_true(string(theta)=='y2',3) = -w2;
Xi_true(string(theta)=='x2',3) = R2^2;
Xi_true(string(theta)=='x2^3',3) = -1;
Xi_true(string(theta)=='x2 y2^2',3) = -1;
Xi_true(string(theta)=='x2',4) = w2;
Xi_true(string(theta)=='y2',4) = R2^2;
Xi_true(string(theta)=='y2^3',4) = -1;
Xi_true(string(theta)=='x2^2 y2',4) = -1;

X_dot_tr = X_dot_tr(:,1:4);
X_dot_att = X_dot_att(:,1:4);
Theta_1tr = Theta_tr(X_tr(:,5)==mu_tr,:);
X_dot_1tr = X_dot_tr(X_tr(:,5)==mu_tr,:);

[error_1tr,success_1tr,l_1tr] = performance(X_dot_1tr,Theta_1tr,Xi_true);
[error_att,success_att,l_att] = performance(X_dot_att,Theta_att,Xi_true);
[error_both,success_both,l_both] = performance([X_dot_1tr; X_dot_att],[Theta_1tr; Theta_att],Xi_true);
[error_full,success_full,l_full] = performance([X_dot_tr; X_dot_att],[Theta_tr; Theta_att],Xi_true);

alphas = logspace(-20,20,11)';
ls_mo = zeros(length(alphas),1);
errors_mo = zeros(length(alphas),1);
for i=1:length(alphas)
    X_dot = [X_dot_1tr; sqrt(alphas(i))*X_dot_att];
    Theta = [Theta_1tr; sqrt(alphas(i))*Theta_att];
    [errors_mo(i),~,ls_mo(i)] = performance(X_dot,Theta,Xi_true);
end
[error_mo,i] = min(errors_mo);
alpha = alphas(i);
l_mo = ls_mo(i);
X_dot = [X_dot_1tr; sqrt(alpha)*X_dot_att];
Theta = [Theta_1tr; sqrt(alpha)*Theta_att];
Xi = STLS(X_dot,Theta,l_mo);
success_mo = isequal(~Xi, ~Xi_true);

error =[error_1tr, error_att, error_both, error_full, error_mo];
success = [success_1tr, success_att, success_both, success_full, success_mo];
lambda = [l_1tr, l_att, l_both, l_full, l_mo];
disp(error);
disp(success);
disp(lambda);
disp(alpha);

save('../data/performance_SL.mat',"Xi_true","error","success","lambda","alpha",'-v7.3');

%% Lorenz system
clear variables; close all; clc
addpath('../src'); % Add the source files to the path
load('../Final_Results/data_Lorenz.mat')

sigma = 10; beta = 8/3;
Xi_true = zeros(length(theta),3);
Xi_true(string(theta)=='y',1) = sigma;
Xi_true(string(theta)=='x',1) = -sigma;
Xi_true(string(theta)=='x u',2) = 1;
Xi_true(string(theta)=='x z',2) = -1;
Xi_true(string(theta)=='y',2) = -1;
Xi_true(string(theta)=='x y',3) = 1;
Xi_true(string(theta)=='z',3) = -beta;

% f = @(x) [sigma*(x(2,:)-x(1,:));
%           x(1,:).*(x(4,:)-x(3,:))-x(2,:);
%           x(1,:).*x(2,:)-beta*x(3,:)];
% X_dot_tr = f(X_tr')';
% X_dot_att = f(X_att')';

% Theta = @(x,y,z,u) [1+0*x,x,y,z,u,x.^2,x.*y,x.*z,x.*u,...
%                    y.^2,y.*z,y.*u,z.^2,z.*u,u.^2,x.^3,...
%                    x.^2 .*y,x.^2 .*z,x.^2 .*u,x.*y.^2,x.*y.*z,...
%                    x.*y.*u,x.*z.^2,x.*z.*u,x.*u.^2,y.^3,y.^2 .*z,...
%                    y.^2 .*u,y.*z.^2,y.*z.*u,y.*u.^2,z.^3,z.^2 .*u,...
%                    z.*u.^2,u.^3,x.^4,x.^3 .*y,x.^3 .*z,x.^3 .*u,...
%                    x.^2 .*y.^2,x.^2 .*y.*z,x.^2 .*y.*u,x.^2 .*z.^2,...
%                    x.^2 .*z.*u,x.^2 .*u.^2,x.*y.^3,x.*y.^2 .*z,...
%                    x.*y.^2 .*u,x.*y.*z.^2,x.*y.*z.*u,x.*y.*u.^2,...
%                    x.*z.^3,x.*z.^2 .*u,x.*z.*u.^2,x.*u.^3,y.^4,...
%                    y.^3 .*z,y.^3 .*u,y.^2 .*z.^2,y.^2 .*z.*u,y.^2 .*u.^2,...
%                    y.*z.^3,y.*z.^2 .*u,y.*z.*u.^2,y.*u.^3,z.^4,...
%                    z.^3 .*u,z.^2 .*u.^2,z.*u.^3,u.^4];
% Theta_tr = Theta(X_tr(:,1),X_tr(:,2),X_tr(:,3),X_tr(:,4));
% Theta_att = Theta(X_att(:,1),X_att(:,2),X_att(:,3),X_att(:,4));

X_dot_tr = X_dot_tr(:,1:3);
X_dot_att = X_dot_att(:,1:3);
Theta_1tr = Theta_tr(X_tr(:,4)==mu_tr,:);
X_dot_1tr = X_dot_tr(X_tr(:,4)==mu_tr,:);

[error_1tr,success_1tr,l_1tr] = performance(X_dot_1tr,Theta_1tr,Xi_true);
[error_att,success_att,l_att] = performance(X_dot_att,Theta_att,Xi_true);
[error_both,success_both,l_both] = performance([X_dot_1tr; X_dot_att],[Theta_1tr; Theta_att],Xi_true);
[error_full,success_full,l_full] = performance([X_dot_tr; X_dot_att],[Theta_tr; Theta_att],Xi_true);

alphas = logspace(-20,20,11)';
ls_mo = zeros(length(alphas),1);
errors_mo = zeros(length(alphas),1);
for i=1:length(alphas)
    X_dot = [X_dot_1tr; sqrt(alphas(i))*X_dot_att];
    Theta = [Theta_1tr; sqrt(alphas(i))*Theta_att];
    [errors_mo(i),~,ls_mo(i)] = performance(X_dot,Theta,Xi_true);
end
[error_mo,i] = min(errors_mo);
alpha = alphas(i);
l_mo = ls_mo(i);
X_dot = [X_dot_1tr; sqrt(alpha)*X_dot_att];
Theta = [Theta_1tr; sqrt(alpha)*Theta_att];
Xi = STLS(X_dot,Theta,l_mo);
success_mo = isequal(~Xi, ~Xi_true);

error =[error_1tr, error_att, error_both, error_full, error_mo];
success = [success_1tr, success_att, success_both, success_full, success_mo];
lambda = [l_1tr, l_att, l_both, l_full, l_mo];
disp(error);
disp(success);
disp(lambda);
disp(alpha);

save('../data/performance_Lorenz.mat',"Xi_true","error","success","lambda","alpha",'-v7.3');
%%

clear variables; close all; clc
addpath('../src'); % Add the source files to the path

c_1tr = [1 0 1]; %'m';
c = [0.4940 0.1840 0.5560];
c_att = [0 0 1]; %'b';
c_full = [0 0 0]; %'k';
colors = [c_1tr; c_att; c; c_full; c; c];
lw = 1.5;

aspect =1.8;
len = 300;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 1.05*len len/aspect])

loaddir = {'../data/performance_SNhys'
           '../data/performance_Hopf'
           '../data/performance_SL'
           '../data/performance_Lorenz'};

savedir = {'../plots/performance_SNhys'
           '../plots/performance_Hopf'
           '../plots/performance_SL'
           '../plots/performance_Lorenz'};

j = 4;
load(loaddir{j});
disp(error);
disp(success);
disp(lambda);
disp(alpha);
%%
error = [error(1:4) 0 error(5)];
b = bar(error,0.65,'FaceColor','flat','EdgeColor','flat');
b.CData = colors;
hold on
plot([5 5], [1e-15 1e10], 'k--')
ylabel('Model coef. error')
set(gca,'Fontsize',16)
set(gca,'YScale','log')
axis([0.3 6.7 1e-5 1e5])
set(gcf,'Color','w')
% grid on
xticks([]);
yticks(10.^(-10:2:10))
% yticks([0,1e-5,1e-3,1e-1, 10, 1e3])
% yticks([0,1e-9,1e-6,1e-3, 1, 1e3])
print(gcf,savedir{j},'-dpng','-r600')