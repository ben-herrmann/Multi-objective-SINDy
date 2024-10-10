%% Saddle-node hysteresys loop

clear variables; close all; clc
addpath('../src'); % Add the source files to the path
load('../Final_Results/data_SNHys.mat')

c_1tr = 'm';
c = [0.4940 0.1840 0.5560];
c_att = 'b';
c_full = 'k';
lw = 1.5;

Theta_1tr = Theta_tr(X_tr(:,2)==mu_tr,:);
K_1tr = my_cond(Theta_1tr);
K_att = my_cond(Theta_att);
K_full = my_cond([Theta_tr; Theta_att]);

q = 50;
alphas = logspace(-20,20,q)';
K = zeros(q);
for i=1:q
    alpha = alphas(i);
    K(i) = my_cond([Theta_1tr; sqrt(alpha)*Theta_att]);
end

aspect =1.5;
len = 300;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 1.05*len len/aspect])

loglog(alphas,K,'color',c,'LineWidth',lw)
hold on
loglog(alphas([1,end]),[K_1tr,K_1tr],'-','color',c_1tr,'LineWidth',lw)
loglog(alphas([1,end]),[K_att,K_att],'-','color',c_att,'LineWidth',lw)
loglog(alphas([1,end]),[K_full,K_full],'-','color',c_full,'LineWidth',lw)
hold off
xlabel('$\alpha$');
ylabel('$\kappa(\Theta)$')
set(gca,'Fontsize',16)
axis([1e-20 1e20 1 1e39])
set(gcf,'Color','w')
% grid on
xticks([1e-20 1 1e20]);
yticks([1e0 1e15 1e30]);
% print(gcf,'../plots/cond_SNHys','-dpng','-r600')

%% Hopf normal form

clear variables; close all; clc
addpath('../src'); % Add the source files to the path
load('../Final_Results/data_hopf.mat')

c_1tr = 'm';
c = [0.4940 0.1840 0.5560];
c_att = 'b';
c_full = 'k';
lw = 1.5;

Theta_1tr = Theta_tr(X_tr(:,3)==mu_tr,:);
K_1tr = my_cond(Theta_1tr);
K_att = my_cond(Theta_att);
K_full = my_cond([Theta_tr; Theta_att]);

q = 50;
alphas = logspace(-20,20,q)';
K = zeros(q,1);
for i=1:q
    alpha = alphas(i);
    K(i) = my_cond([Theta_1tr; sqrt(alpha)*Theta_att]);
end

aspect =1.5;
len = 300;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 1.05*len len/aspect])
loglog(alphas([1,end]),[K_full,K_full],'-','color',c_full,'LineWidth',lw)
hold on
loglog(alphas([1,end]),[K_1tr,K_1tr],'-','color',c_1tr,'LineWidth',lw)
loglog(alphas([1,end]),[K_att,K_att],'-','color',c_att,'LineWidth',lw)
loglog(alphas,K,'color',c,'LineWidth',lw)
hold off
xlabel('$\alpha$');
ylabel('$\kappa(\Theta)$')
set(gca,'Fontsize',16)
axis([1e-20 1e20 1e2 1e18])
set(gcf,'Color','w')
% grid on
% xticks([]);
yticks([1e5 1e10 1e15]);
% print(gcf,'../plots/cond_Hopf','-dpng','-r600')

%% Stuart–Landau oscillators

clear variables; close all; clc
addpath('../src'); % Add the source files to the path
load('../Final_Results/data_SL.mat')

c_1tr = 'm';
c = [0.4940 0.1840 0.5560];
c_att = 'b';
c_full = 'k';
lw = 1.5;

Theta_1tr = Theta_tr(X_tr(:,5)==mu_tr,:);
K_1tr = my_cond(Theta_1tr);
K_att = my_cond(Theta_att);
K_full = my_cond([Theta_tr; Theta_att]);

q = 50;
alphas = logspace(-20,20,q)';
K = zeros(q);
for i=1:q
    alpha = alphas(i);
    K(i) = my_cond([Theta_1tr; sqrt(alpha)*Theta_att]);
end

aspect =1.5;
len = 300;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 1.05*len len/aspect])

loglog(alphas([1,end]),[K_full,K_full],'-','color',c_full,'LineWidth',lw)
hold on
loglog(alphas([1,end]),[K_1tr,K_1tr],'-','color',c_1tr,'LineWidth',lw)
loglog(alphas([1,end]),[K_att,K_att],'-','color',c_att,'LineWidth',lw)
loglog(alphas,K,'color',c,'LineWidth',lw)
hold off
xlabel('$\alpha$');
ylabel('$\kappa(\Theta)$')
set(gca,'Fontsize',16)
axis([1e-20 1e20 1e4 1e17])
set(gcf,'Color','w')
% grid on
xticks([1e-20 1 1e20]);
yticks([1e5 1e10 1e15]);
% print(gcf,'../plots/cond_SL','-dpng','-r600')

%% Lorenz system

clear variables; close all; clc
addpath('../src'); % Add the source files to the path
load('../Final_Results/data_lorenz.mat')

c_1tr = 'm';
c = [0.4940 0.1840 0.5560];
c_att = 'b';
c_full = 'k';
lw = 1.5;

% Theta_att = [Theta_att(1:2:16,:);Theta_att(17:end,:)];
Theta_1tr = Theta_tr(X_tr(:,4)==mu_tr,:);
% Theta_1tr = Theta_1tr(1:300,:);
K_1tr = my_cond(Theta_1tr);
K_att = my_cond(Theta_att);
K_full = my_cond([Theta_tr; Theta_att]);

q = 50;
alphas = logspace(-20,20,q)';
K = zeros(q);
for i=1:q
    alpha = alphas(i);
    K(i) = my_cond([Theta_1tr; sqrt(alpha)*Theta_att]);
end

aspect =1.5;
len = 300;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 1.05*len len/aspect])

loglog(alphas([1,end]),[K_full,K_full],'-','color',c_full,'LineWidth',lw)
hold on
loglog(alphas([1,end]),[K_1tr,K_1tr],'-','color',c_1tr,'LineWidth',lw)
loglog(alphas([1,end]),[K_att,K_att],'-','color',c_att,'LineWidth',lw)
loglog(alphas,K,'color',c,'LineWidth',lw)
hold off
xlabel('$\alpha$');
ylabel('$\kappa(\Theta)$')
set(gca,'Fontsize',16)
axis([1e-20 1e20 1e8 1e23])
set(gcf,'Color','w')
% grid on
xticks([1e-20 1 1e20]);
yticks([1e10 1e15 1e20]);
print(gcf,'../plots/cond_Lorenz','-dpng','-r600')
