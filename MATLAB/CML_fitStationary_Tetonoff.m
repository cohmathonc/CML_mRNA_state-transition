clear


%%%%%%%%%%%%%%%%%%%%%%May 5th 2023 values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [-208.51554, -43.88014, 111.85604, 138.55368]'; %rotated space only nonzero genes  
%Tet off/on
%x0 = [-135.983923, 3.012764, 112.561001]';

%remap [-129 19] interval to [0 1] with linear transformation
%[-129, 0], [19, 1]
A = 1;
B = 0; 

m = (B-A) / (x0(end)-x0(1));

x = m.*(x0-x0(1)) + A;

%%%%%%%%%%%%%%%%%%%%%%May 5th 2023 values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tetfactor = 2.5;
f = [0.002204351, 0.002422282, 0.001310297, 0.001322798]'; %rotated space only nonzero genes
%Tet off/on
%f = [0.003155247, 0.001223097, 0.003697754]';


y = -log(f);

c3 = x0(3);
%Up 6th degree polynomial
X0 = [x.^6, x.^5, x.^4, x.^3, x.^2, x, ones(size(x))]; %F
F = @(x,P)P(1)*x.^6 + P(2)*x.^5 + P(3)*x.^4 + P(4)*x.^3 + P(5)*x.^2 + P(6)*x + P(7);

%dUp 5th degree polynomial gradient spatial derivative of Up
X1 = [6*x.^5, 5*x.^4, 4*x.^3, 3*x.^2, 2*x, ones(size(y)), zeros(size(x))]; %G
G = @(x,P)6*P(1)*x.^5 + 5*P(2).*x.^4 + 4*P(3).*x.^3 + 3*P(4).*x.^2 + 2*P(5).*x + P(6);

rhs = [zeros(size(x)); y]; %G(X) = 0, F(X) = -ln(X)
a = [X1;X0]\rhs; %solve



X = linspace(-.15,1.15);

figure(1);
hold on
subplot(1,3,1)
plot(x,zeros(size(x)),'o') %plot critical points as zeros of G(X)
hold on
plot(X,G(X,a),'b') %plot identified G(x)
plot(X,zeros(size(X)),'r')
title('G(X)=0')
set(gca,'xtick',sort(x),'xticklabel',{'healthy c_1','c_3','c_4','CML c_5'})
subplot(1,3,2)
plot(x,y,'o') %plot critical points evaluated at potential
hold on
plot(X,F(X,a),'b')
title('F(X)=-ln(x)=U_p(x)')
set(gca,'xtick',sort(x),'xticklabel',{'healthy c_1','c_3','c_4','CML c_5'})

subplot(1,3,3)
plot(x,f,'o')
hold on
plot(X,exp(-F(X,a)),'b')
title('stationary distribution')
set(gca,'xtick',sort(x),'xticklabel',{'healthy c_1','c_3','c_4','CML c_5'})

%%
%% Simulation of Langevin equation
alpha = 1e-2; 
t0 = 0;
tf = 6;


t0tki = 6; 
tftki = 9;
tkidt = tftki - t0tki;


t0post = 9;
tfpost = 16;
postdt = tfpost - t0post;


Ntrials = 20;
nPeriods = 100;

D = 0.008;

Pot = @(t,X) -alpha*G(X,a);
Diff = @(t,X) sqrt(2*D)*ones(size(X));

TKIPot = @(t,X) 4*alpha*G(X,a);
TKIDiff = @(t,X) .8*sqrt(2*D)*ones(size(X));

PostPot = @(t,X) 2*alpha*G(X,a);

langevin = sde(Pot, Diff,'StartState',-.1);
[S,T] = simulate(langevin, round(nPeriods), 'DeltaTime', tf/nPeriods, 'nTrials', Ntrials);
S = squeeze(S);
T = linspace(t0,tf,nPeriods+1);

langevin = sde(TKIPot, TKIDiff,'StartState',S(end,:));
[S2,T2] = simulate(langevin, round(nPeriods), 'DeltaTime', tkidt/nPeriods, 'nTrials', Ntrials);
S2 = squeeze(S2);

S = [S; S2];
T = [T linspace(t0tki,tftki,nPeriods+1)];

langevin = sde(PostPot, Diff,'StartState',S(end,:));
[S3,T3] = simulate(langevin, round(nPeriods), 'DeltaTime', postdt/nPeriods, 'nTrials', Ntrials);
S3 = squeeze(S3);

S = [S; S3];

T = [T linspace(t0post,tfpost,nPeriods+1)];

figure
%subplot(1,2,1)
hold on; box on;
plot(T,S); 
plot([t0 tfpost],[x(1) x(1)],'k')
plot([t0 tfpost],[x(2) x(2)],'k')
plot([t0 tfpost],[x(3) x(3)],'k')
plot([t0 tfpost],[x(4) x(4)],'k')
plot([t0tki t0tki],[0 1],':k')
plot([tftki tftki],[0 1],':k')
xlabel('time (weeks)')
set(gca,'ylim',[-.25 1.25],'xlim',[t0 tfpost],'ytick',sort(x),'yticklabel',{'healthy c_1','c_3','c_4','CML c_5'})
ax = gca; 
%ax.ColorOrder = mycolors;
% 
% subplot(1,2,2)
% hold on; box on;
% plot(S(end,:),alpha*F(S(end,:),a),'o')
% plot(x,alpha*y,'o') %plot critical points evaluated at potential
% plot(X,alpha*F(X,a),'k')
% title('U_p(x)')
% set(gca,'xtick',sort(x),'xticklabel',{'healthy c_1','c_3','c_4','CML c_5'})
% colororder(gca,ax.ColorOrder)

%%
% p = plot(nan,nan);
% p.XData = t_step;
% for n = 1:NPeriods
%       p.YData = x.^n;
%       exportgraphics(gcf,'testAnimated.gif','Append',true);
% end