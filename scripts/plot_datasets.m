%% Saddle-node hysteresis loop
clear variables; close all; clc
load('../data/data_SNHys.mat')
disp(str2double(alpha_exp));
disp(lambda);
%%
m = 1000*ones(1,length(Mu));
r = ones(1,length(Mu));


aspect =1.6;
lw = 1.6;
gray = 0.6*[1,1,1];
alpha = 0.3;
c1 = [0.4940 0.1840 0.5560];
c2 = 'b';
len = 400;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 1.05*len len/aspect])
xx = linspace(-2,2,100);
mmu = -xx+xx.^3;
muc = 1/sqrt(3);
xc = fsolve(@(x)muc+x-x^3,-1);
plot(mmu(xx<xc),xx(xx<xc),mmu(xx>-xc),xx(xx>-xc),mmu(abs(xx)<-xc),xx(abs(xx)<-xc),'--','color',[gray 1.5*alpha],'LineWidth',1.5*lw)
xlabel('$\mu$')
ylabel('$x$')
set(gca,'Fontsize',16)%,'XAxisLocation','origin')
hold on

k = length(Mu);
for i=1:k

    i0_tr = 1+sum(m(1:i-1));
    iend_tr = sum(m(1:i));
    x = X_tr(i0_tr:iend_tr,1);
    mu = X_tr(i0_tr:iend_tr,2);
    if mu == mu_tr
        plot(mu,x,'LineWidth',lw,'Color','m');
    % else
    %     plot(mu,x,'LineWidth',lw,'Color',c1);
    end

    i0_att = 1+sum(r(1:i-1));
    iend_att = sum(r(1:i));
    x = X_att(i0_att:iend_att,1);
    mu = X_att(i0_att:iend_att,2);
    plot(mu,x,'o','LineWidth',1,'Color',c2,'MarkerFaceColor',c2)
end
hold off
axis([-6.3,6.3,-2.1,2.1])
% axis off
% view(-80,35)
% camlight('headlight')
% lighting gouraud
set(gcf,'Color','w')
% grid on
xticks([]); yticks([]);
% print(gcf,'../plots/data_SNHys_1tr','-dpng','-r600')

%% Hopf normal form
clear variables; close all; clc
load('../data/data_hopf.mat')
disp(str2double(alpha_exp));
disp(lambda);
%%
m = double(m)*ones(1,length(Mu));
aspect =1.4;
lw = 1.6;
gray = 0.6*[1,1,1];
alpha = 0.3;
c1 = [0.4940 0.1840 0.5560];
c2 = 'b';
len = 400;
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 1.05*len len/aspect])

plot3([0,0],[-0.5,0],[0,0],'color',[gray 2*alpha],'LineWidth',2*lw)
xlabel('$x$')
ylabel('$\mu$')
zlabel('$y$')
set(gca,'Fontsize',16,'Ydir', 'reverse')%,'XAxisLocation','origin')
hold on

xx = linspace(-1,1,100);
yy = linspace(-1,1,100);
[XX,YY] = meshgrid(xx,yy);
MMu = XX.^2 + YY.^2;
surf(XX,MMu,YY,'faceColor',gray,'faceAlpha',alpha,'linestyle','none')

k = length(Mu);
for i=1:k

    i0_tr = 1+sum(m(1:i-1));
    iend_tr = sum(m(1:i));
    x = X_tr(i0_tr:iend_tr,1);
    y = X_tr(i0_tr:iend_tr,2);
    mu = X_tr(i0_tr:iend_tr,3);
    disp([x(1) y(1) mu(1)]);
    if mu == mu_tr
        plot3(x,mu,y,'LineWidth',lw,'Color','m');
    else 
        plot3(x,mu,y,'LineWidth',lw,'Color',[c1 0.6]);
    end

    if i<k
    i0_att = 1+sum(r(1:i-1));
    iend_att = sum(r(1:i));
    x = X_att(i0_att:iend_att,1);
    y = X_att(i0_att:iend_att,2);
    mu = X_att(i0_att:iend_att,3);
    plot3(x,mu,y,'o','LineWidth',1,'Color',c2,'MarkerFaceColor',c2)
    end
end

hold off
axis([-1,1,-0.4,0.65,-1,1])
% axis off
view(-80,35)
camlight('headlight')
lighting gouraud
set(gcf,'Color','w')
% grid on
% xticks([]); yticks([]); zticks([]);
% print(gcf,'../plots/data_Hopf_1tr','-dpng','-r600')

%% Coupled Stuart–Landau oscillators
clear variables; close all; clc;
load('../data/data_SL.mat')
disp(str2double(alpha_exp));
disp(lambda);
%%
k = length(Mu);
lw = 1.8;
gray = 0.6*[1,1,1];
alpha = 0.08;
c1 = [0.4940 0.1840 0.5560];
c2 = 'b';

is = (1:k);

aspect =1*length(is);
len = 200*length(is);
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 1.05*len len/aspect])

for i=is

R1 = 1; R2 = 0.2; w1 = 1/pi; w2 =1;
muc = abs(w2-2*w1)/(2*R2);
dt = 0.4;
t = 0:dt:30000-dt;
x0 = [R1; 0; R2; 0];
mu = Mu(i);
f = @(t,x) [-w1*x(2) + (R1^2-(x(1)^2+x(2)^2))*x(1) + mu*(x(1)*x(3)+x(2)*x(4));
    w1*x(1) + (R1^2-(x(1)^2+x(2)^2))*x(2) + mu*(x(1)*x(4)-x(2)*x(3));
    -w2*x(4) + (R2^2-(x(3)^2+x(4)^2))*x(3);
    w2*x(3) + (R2^2-(x(3)^2+x(4)^2))*x(4)];
[~,X] = ode78(f,t,x0);
X = X(end-1000:end,:);


subplot(1,length(is),i)
plot3(X(:,1),X(:,2),X(:,3),'color',[gray alpha],'LineWidth',2*lw)
xlabel('$x_1$')
ylabel('$y_1$')
zlabel('$x_2$')
set(gca,'Fontsize',16,'Ydir', 'reverse')%,'XAxisLocation','origin')
hold on

i0_tr = 1+sum(m(1:i-1));
iend_tr = sum(m(1:i));
x = X_tr(i0_tr:iend_tr,1);
y = X_tr(i0_tr:iend_tr,2);
z = X_tr(i0_tr:iend_tr,3);
mu = X_tr(i0_tr:iend_tr,5);
disp([x(1) y(1) z(1) X_tr(i0_tr,4)]);
if mu == mu_tr
    plot3(x,y,z,'LineWidth',lw,'Color','m');
else
    plot3(x,y,z,'LineWidth',lw,'Color',c1);
end

i0_att = 1+sum(r(1:i-1));
iend_att = sum(r(1:i));
x = X_att(i0_att:iend_att,1);
y = X_att(i0_att:iend_att,2);
z = X_att(i0_att:iend_att,3);
plot3(x,y,z,'o','LineWidth',1,'Color',c2,'MarkerFaceColor',c2)

hold off
axis(1.2*[-R1 R1 -R1 R1 -R2 R2])
% axis off
view(-60,55)
camlight('headlight')
lighting gouraud
set(gcf,'Color','w')
grid on
% xticks([]); yticks([]); zticks([]);
title(['$\mu=$' num2str(round(Mu(i),1))])

if i~=is(1)
axis off
end

end

% print(gcf,'../plots/data_SL','-dpng','-r1800')

%% Lorenz system
clear variables; close all; clc;
load('../data/data_lorenz.mat')
disp(str2double(alpha_exp));
disp(lambda);
%%
Mu = double(Mu);
k = length(Mu);
r = [2 2 2 2 2 2 2 2 25 25 25 25];
m = 500*ones(1,k);

lw = 1.8;
gray = 0.6*[1,1,1];
alpha = 0.1;
c1 = [0.4940 0.1840 0.5560];
c2 = 'b';

is = (1:k);

aspect =1*length(is);
len = 200*length(is);
f1 = figure('DefaultTextInterpreter','Latex','DefaultAxesTickLabelInterpreter','Latex');
set(f1,'Position',[-1800 1000 1.05*len len/aspect])

sigma = 10; beta = 8/3;
muc = 24.74;
dt = 0.02;
t = 0:dt:20000-dt;
x0 = [10; 10; 10];

for i=is

mu = Mu(i);

subplot(1,length(is),i)

if mu>muc
    f = @(t,x) [sigma*(x(2)-x(1)); x(1)*(mu-x(3))-x(2); x(1)*x(2)-beta*x(3)];
    [~,X] = ode78(f,t,x0);
    X = X(end-1000:end,:);
    plot3(X(:,1),X(:,2),X(:,3),'color',[gray alpha],'LineWidth',2*lw)
else
    plot3([-sqrt(beta*(mu-1)),sqrt(beta*(mu-1))],[-sqrt(beta*(mu-1)),sqrt(beta*(mu-1))],[mu-1,mu-1],'o','color',[gray 5*alpha],'MarkerFaceColor',gray)
end
    xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
set(gca,'Fontsize',16,'Ydir', 'reverse')%,'XAxisLocation','origin')
hold on

i0_tr = 1+sum(m(1:i-1));
iend_tr = sum(m(1:i));
x = X_tr(i0_tr:iend_tr,1);
y = X_tr(i0_tr:iend_tr,2);
z = X_tr(i0_tr:iend_tr,3);
mu = X_tr(i0_tr:iend_tr,4);
disp([x(1) y(1) z(1)]);
if mu == mu_tr
    plot3(x,y,z,'LineWidth',lw,'Color','m');
% else
%     plot3(x,y,z,'LineWidth',lw,'Color',c1);
end

i0_att = 1+sum(r(1:i-1));
iend_att = sum(r(1:i));
x = X_att(i0_att:iend_att,1);
y = X_att(i0_att:iend_att,2);
z = X_att(i0_att:iend_att,3);
plot3(x,y,z,'o','LineWidth',0.1*lw,'Color',c2,'MarkerFaceColor',c2,'MarkerSize',5)


hold off
axis([-30 30 -50 20 0 60])
view(-30,10)
camlight('headlight')
lighting gouraud
set(gcf,'Color','w')
title(['$\mu=$' num2str(round(Mu(i),1))])

if i~=is(1)
axis off
end

% grid on
xticks([]); yticks([]); zticks([]);

end

print(gcf,'../plots/data_Lorenz_1tr','-dpng','-r1800')

