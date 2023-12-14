clear

 %x0 = [-129, -108, -66, -33, 19]';
 %x0 = [-129, -108,  -66, 19]'; %remove unstable point c2
 %x0 = [-129, -66, -33, 19]'; %remove unstable point c4

%x0 = [-96.16223, -66.43368, -17.11132, 22.75195, 75.45256]'; %rotated space 
%x0 = [-96.16223, -66.43368, -17.11132,  75.45256]'; %rotated space remove c2

%x0 = [-99.31529, -69.45024, -19.42626, 38.06397, 70.16891]'; %rotated space adjusted bandwidth
%x0 = [-99.31529, -69.45024, -19.42626, 70.16891]'; %rotated space adjusted bandwidth 

%%%%%%%%%%%%%%%%%%%%%%May 5th 2023 values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x0 = [-208.51554, -125.45624, -43.88014, 111.85604, 138.55368]'; %rotated space only nonzero genes 
x0 = [-208.51554, -43.88014, 111.85604, 138.55368]'; %rotated space only nonzero genes  

%TKI
%x0 = [-266.27847, -249.44004, -160.75767, 32.32295, 137.84375]';
%x0 = [-266.27847, -249.44004, -160.75767, 137.84375]';

%Tet off/on
%x0 = [-135.983923, 3.012764, 112.561001]';

%remap [-129 19] interval to [0 1] with linear transformation
%[-129, 0], [19, 1]
A = 1;
B = 0; 

m = (B-A) / (x0(end)-x0(1));

x = m.*(x0-x0(1)) + A;


%f = [0.00315, 0.00225, 0.00815, 0.00325, 0.00485]';
% outlierfactor = 1; 
%f = [0.00315, 0.00225, 0.00815,  outlierfactor*0.00485]'; %remove unstable point c2
%f = [0.00315, 0.00815, 0.00325, 0.00485]'; %remove unstable point c4

%f = [0.004354416, 0.003036103, 0.007902762, 0.002272277, 0.003237811]'; %rotated space 
%f = [0.004354416, 0.003036103, 0.007902762, 0.003237811]'; %rotated space 

%f = [0.004047431, 0.75*0.003541261, 0.006633756, 0.75*0.002427884, 0.002823594]'; %rotated space adjusted bandwidth
%f = [0.004047431, 0.003541261, 0.006633756, 0.002823594]'; %rotated space adjusted bandwidth


%%%%%%%%%%%%%%%%%%%%%%May 5th 2023 values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%f = [0.002204351, 0.002072662, 0.002422282, 0.001310297, 0.001322798]'; %rotated space only nonzero genes
f = [0.002204351, 0.002422282, 0.001310297, 0.001322798]'; %rotated space only nonzero genes

%TKI
%f = [0.0019949871, 0.0019259394, 0.0054708354, 0.0002445441, 0.0005553807]';
%f = [0.0019949871, 0.0019259394, 0.0054708354, 0.0005553807]';

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


ang = (180/(2*pi()))*acos((X0*X1')./(norm(X0)*norm(X1))); %angle between Up and dUp
 c2 = fzero(@(x) G(x,a),0.3) 
 OGc2 = ((c2 - A)./m) + x0(1)

 c4 = fzero(@(x) G(x,a),0.7) 
 OGc4 = ((c4 - A)./m) + x0(1)

X = linspace(-.15,1.15);

figure(1);
hold on
subplot(1,3,1)
plot(x,zeros(size(x)),'o') %plot critical points as zeros of G(X)
hold on
plot(X,G(X,a),'b') %plot identified G(x)
plot(X,zeros(size(X)),'r')
title('G(X)=0')
set(gca,'xtick',sort(x),'xticklabel',{'healthy c_1','c_2','c_3','CML c_5'})
subplot(1,3,2)
plot(x,y,'o') %plot critical points evaluated at potential
hold on
plot(X,F(X,a),'b')
title('F(X)=-ln(x)')
set(gca,'xtick',sort(x),'xticklabel',{'healthy c_1','c_2','c_3','CML c_5'})

subplot(1,3,3)
plot(x,f,'o')
hold on
plot(X,exp(-F(X,a)),'b')
title('stationary distribution')
set(gca,'xtick',sort(x),'xticklabel',{'healthy c_1','c_2','c_4','CML c_5'})

%%
%% Simulation of Langevin equation
alpha = 1e-2; 
t0 = 0;
tf = 20;
Ntrials = 10;
nPeriods = 100;
t_step = linspace(t0,tf,nPeriods+1);
D = 0.012;

Pot = @(t,X) -alpha*G(X,a);
Diff = @(t,X) sqrt(2*D);

langevin = sde(Pot, Diff,'StartState',0);
[S,T] = simulate(langevin, round(nPeriods), 'DeltaTime', tf/nPeriods, 'nTrials', Ntrials);

figure
subplot(1,2,1)
hold on; box on;
plot(T,squeeze(S)); 
plot([t0 tf],[x(1) x(1)],'k')
plot([t0 tf],[x(2) x(2)],'k')
plot([t0 tf],[x(3) x(3)],'k')
plot([t0 tf],[x(4) x(4)],'k')

xlabel('t')
ylabel('x')
set(gca,'ylim',[-.5 1.25],'xlim',[t0 tf],'ytick',sort(x),'yticklabel',{'healthy c_1','c_2','c_3','CML c_5'})
ax = gca; 
%ax.ColorOrder = mycolors;

subplot(1,2,2)
hold on; box on;
plot(squeeze(S(end,1,:)),alpha*F(squeeze(S(end,1,:)),a),'o')
plot(x,alpha*y,'o') %plot critical points evaluated at potential
plot(X,alpha*F(X,a),'k')
title('U_p(x)')
set(gca,'xtick',sort(x),'xticklabel',{'healthy c_1','c_2','c_3','CML c_5'})
colororder(gca,ax.ColorOrder)

%%
% p = plot(nan,nan);
% p.XData = t_step;
% for n = 1:NPeriods
%       p.YData = x.^n;
%       exportgraphics(gcf,'testAnimated.gif','Append',true);
% end