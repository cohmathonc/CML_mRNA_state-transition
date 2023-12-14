%% parameter setting
clear
%gpuDevice
global gamma c1 c2 c3 c4 c5 alpha D Dp lam1 lam2 ti tj Dc Dd m ini_max ini_min dt2 ti2 tj2 ti
sympref('HeavisideAtOrigin',1);

m = 0;

t0 = 0;
tf = 18;
c1 = 138.55368;
c2 = 111.85604;
c3 = -43.88014;
c4 = -125.45624;
c5 = -208.51554; 
ti = [6 9];
% tki half-life 17 hr=0.1 week
lam1 = log(2)/0.1; 
t = linspace(t0,tf,tf);
alpha = 1e2;
D = 0.0375;

A=1;
B=0;
mn=(B-A)/(c5-c1);
Y = @(x)(mn.*(x-c1))+A;
c1 = Y(c1); c2 = Y(c2); c3 = Y(c3); c4=Y(c4); c5=Y(c5);
x = linspace(-0.2,1.2);
dx = x(2)-x(1);

% for SDE trajectories
Ntrials = 300;
dt = 1/(24*7);
nPeriods = round(18/(dt));

cmldata = readtable("CML.mRNA_state-space_sampleInfo_table.tsv", "FileType","text",'Delimiter','\t');
cmldata_cml_space = Y(cmldata.CML_space);

%'TET_OFF_B', 'TET_OFF_NIL_ON_D', 'TET_OFF_ON_A', 'TET_ON_C'
tet_off_on_a = find(contains(cmldata.treatment,'TET_OFF_ON_A'));
tet_off_b = find(contains(cmldata.treatment,'TET_OFF_B'));
tet_on_c = find(contains(cmldata.treatment,'TET_ON_C'));
tet_off_nil_on_d = find(contains(cmldata.treatment,'TET_OFF_NIL_ON_D'));

%mid
tet_off_b_mid = unique(cmldata.mouse_id(tet_off_b));
tet_on_c_mid = unique(cmldata.mouse_id(tet_on_c));
tet_off_nil_on_d_mid = unique(cmldata.mouse_id(tet_off_nil_on_d));
tet_off_on_a_mid = unique(cmldata.mouse_id(tet_off_on_a));

%% potential function
fCML.p1 = -alpha*c1*c2*c3*c4*c5;
fCML.p2 = alpha*(c1*c3*c4*c5 + c2*c3*c4*c5 + c1*c2*(c4*c5 + c3*(c4 + c5)));
fCML.p3 = -alpha*(c1*c4*c5 + c2*c4*c5 + c3*c4*c5 + c1*c3*(c4 + c5) +  c2*c3*(c4 + c5) + c1*c2*(c3 + c4 + c5));
fCML.p4 = alpha*(c3*c4 + c3*c5 + c4*c5 + c2*(c3 + c4 + c5) +  c1*(c2 + c3 + c4 + c5));
fCML.p5 = -alpha*(c1 + c2 + c3 + c4 + c5); 
fCML.p6 = alpha;

LUp = @(x,myf)1*(myf.p1 + myf.p2*x + myf.p3*x.^2 + myf.p4*x.^3 + myf.p5*x.^4 + myf.p6*x.^5);
LU = @(x,myf)1*(myf.p1*x + (myf.p2/2)*x.^2 + (myf.p3/3)*x.^3 + (myf.p4/4)*x.^4 + (myf.p5/5)*x.^5 + (myf.p6/6)*x.^6);

gamma_inc = 1:1:15;
D_inc = 0.001:0.001:0.01;

%gamma_inc = 9;
gamma_inc = 0;
D_inc = 0.01;

Si_dist_g_d_m = zeros(size(D_inc,2),size(gamma_inc,2));
Si_dist_g_d_p = zeros(size(D_inc,2),size(gamma_inc,2));

Si_diff = [];
Si_diff_m = [];
Si_diff_e = [];

g_inc = 1;
for gamma=gamma_inc
d_inc = 1;
for Dp=D_inc

%% FP solution

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
% Extract the first solution component as u.
u = sol(:,:,1);
% scale to be PDF
u = u.*(1./sum(u,2));
%% plot all time solutions
figure(100); 
subplot(1,3,1)
plot(x,LU(x,fCML)-gamma*LU(x,fCML));
ylabel('$U_p-\gamma*U_p$','Interpreter','latex')
set(gca,'xdir','reverse','fontsize',14,'xtick',[c5 c4 c3 c2 c1 ],'xticklabel',{'c_5','c_4','c_3','c_2','c_1'})
subplot(1,3,2)
plot(x,u(1,:),'b')
hold on;
plot(x,u(end,:),'r')
ylabel('$P(x,t)$','Interpreter','latex')
set(gca,'xdir','reverse','fontsize',14,'xtick',[c5 c4 c3 c2 c1 ],'xticklabel',{'c_5','c_4','c_3','c_2','c_1'})
subplot(1,3,3)
plot(x,exp(-(LU(x,fCML)-gamma*LU(x,fCML))))
ylabel('$\exp(-(U_p- \gamma*U_p) )$','Interpreter','latex')
set(gca,'xdir','reverse','fontsize',14,'xtick',[c5 c4 c3 c2 c1 ],'xticklabel',{'c_5','c_4','c_3','c_2','c_1'})
%saveas(gcf,strcat('FP_Up_iesolutions_eUp',num2str(gamma),'.png'))
%% all time solutions
figure(110)
plot(x,u')
xlabel('x'); ylabel('P')
set(gca,'xdir','reverse','fontsize',14,'xtick',[c5 c4 c3 c2 c1 ],'xticklabel',{'c_5','c_4','c_3','c_2','c_1'})
%saveas(gcf,strcat('FP_all_time_solutions',num2str(gamma),'.png'))
%% surface plot
figure(90)
h = surfc(x,t(2:end),u(2:end,:));
c = h(2);
c.ContourZLevel = -.03;
view(-150,21)
%view(-150,41)
xlabel('CML state')
zlabel('Probability Density')
ylabel('Time (Weeks)')
pbaspect([1 1 1])
set(gca,'fontsize',18,'xtick',[c5 c4 c3 c2 c1 ],'xticklabel',{'c_5','c_4','c_3','c_2','c_1'},'ytick',0:3:tf)
%saveas(gcf,strcat('FPSurface',num2str(gamma),'.png'))
%}
%% check trajectories
% adjust beta value

% TET OFF-ON
G = @(t,X) sqrt( 2*( D +(Dp-D)*heaviside(t-ti(1)) ) );
F = @(t,X) -LUp(X,fCML)+gamma*(heaviside(t-ti(1)))*LUp(X,fCML);
% TKI ON-OFF
%F = @(t,X) -LUp(X,fCML) + gamma*exp(-lam1*(t-ti(2))*heaviside(t-ti(2)))*(heaviside(t-ti(1)))*LUp(X,fCML); 
langevin = sde(F, G,'StartState',c1);
[Si,T] = simulate(langevin, nPeriods, 'DeltaTime', dt, 'nTrials', Ntrials);
Si = Si(:,all(~isnan(Si)));
Si = Si(:,sum(Si <= 1.4 & Si > -0.4) == size(Si,1));

fprintf('P(Si):%f, D:%f, gamma:%f\n',size(Si,2)/Ntrials,Dp,gamma)
Si_dist_g_d_p(d_inc,g_inc) = size(Si,2)/Ntrials;

figure(101);
clf;
[d,ix] = min(abs(T-5));
%plot(T,Si,'Color',[.8 .8 .8]); hold on;
plot(T,Si(:,Si(ix,:) < c4),'Color',[.8 .8 .8]); hold on;

tet_off_on_a_mid = unique(cmldata.mouse_id(tet_off_on_a));
tet_off_b_mid = unique(cmldata.mouse_id(tet_off_b));
tet_on_c_mid = unique(cmldata.mouse_id(tet_on_c));
tet_off_nil_on_d_mid = unique(cmldata.mouse_id(tet_off_nil_on_d));

for i =1:size(tet_off_on_a_mid,1)
    ind = find(cmldata.mouse_id==tet_off_on_a_mid(i));
    plot(cmldata.timepoint(ind),cmldata_cml_space(ind),'b'); hold on;
end
scatter(cmldata.timepoint(tet_off_on_a),cmldata_cml_space(tet_off_on_a),'b'); hold on;
%}
%{
for i =1:size(tet_off_nil_on_d_mid,1)
    ind = find(cmldata.mouse_id==tet_off_nil_on_d_mid(i));
    plot(cmldata.timepoint(ind),cmldata_cml_space(ind),'g'); hold on;
end
scatter(cmldata.timepoint(tet_off_nil_on_d),cmldata_cml_space(tet_off_nil_on_d),'g'); hold on;
%}
for i =1:size(tet_off_b_mid,1)    
    ind = find(cmldata.mouse_id==tet_off_b_mid(i));
    plot(cmldata.timepoint(ind),cmldata_cml_space(ind),'r'); hold on;
end
scatter(cmldata.timepoint(tet_off_b),cmldata_cml_space(tet_off_b),'r'); hold on;

for i =1:size(tet_on_c_mid,1)
    ind = find(cmldata.mouse_id==tet_on_c_mid(i));
    plot(cmldata.timepoint(ind),cmldata_cml_space(ind),'k'); hold on;
end
scatter(cmldata.timepoint(tet_on_c),cmldata_cml_space(tet_on_c),'k'); hold on;

yline(c1,':k','fontweight','bold','fontsize',16,'LineWidth',3); hold on;
yline(c2,':g','fontweight','bold','fontsize',16,'LineWidth',3); hold on;
yline(c3,':c','fontweight','bold','fontsize',16,'LineWidth',3); hold on;
yline(c4,':y','fontweight','bold','fontsize',16,'LineWidth',3); hold on;
yline(c5,':m','fontweight','bold','fontsize',16,'LineWidth',3); hold on;
xline(ti(1),'k','fontweight','bold','fontsize',16,'LineWidth',1.5); hold on;
%xline(ti(2),'k','fontweight','bold','fontsize',16,'LineWidth',1.5); hold on;
xlim([0 18])
ylabel('CML state')
xlabel('Time (Weeks)')
%saveas(gcf,strcat('FP_trajectories_',num2str(gamma),'data.png'))

% TET OFF ON
cml_treatment_t = cmldata.timepoint(tet_off_on_a);
cml_treatment_s = cmldata_cml_space(tet_off_on_a);

% TKI ON OFF
%cml_treatment_t = cmldata.timepoint(tet_off_nil_on_d);
%cml_treatment_s = cmldata_cml_space(tet_off_nil_on_d);

% Remove timepoints before treatment
cml_treatment_t_timepoint = unique(cml_treatment_t);
cml_treatment_t_timepoint = cml_treatment_t_timepoint(cml_treatment_t_timepoint > 5);

% Differences between simulated and ovservved data
Si_diff_i = [];
for si_i = 1:size(cml_treatment_t_timepoint,1)
    [d,ix] = min(abs(T-cml_treatment_t_timepoint(si_i)));
    si_data = cml_treatment_s(cml_treatment_t == cml_treatment_t_timepoint(si_i));

    for sd = 1:size(si_data,1)
        dist_si_data = abs(Si(ix,:) - si_data(sd));
        Si_diff_i = [Si_diff_i; dist_si_data];
    end

end

Si_diff_m = [Si_diff_m mean2(Si_diff_i)];
Si_diff_e = [Si_diff_e std2(Si_diff_i)];

Si_dist_g_d_m(d_inc,g_inc) = mean2(Si_diff_i);
d_inc = d_inc + 1;

end
g_inc = g_inc +  1;

end

figure(113)
imagesc(gamma_inc,D_inc,Si_dist_g_d_m)
title('Average Distance |CML-X|')
xlabel('$\gamma$','Interpreter','latex')
ylabel('$D_{post}$','Interpreter','latex')
colormap(jet)
colorbar
[Dm_ind,gammam_ind] = find(Si_dist_g_d_m==min(min(Si_dist_g_d_m(:,2:end))));
fprintf('Min. Dist. D:%f, gamma:%f\n',D_inc(Dm_ind),gamma_inc(gammam_ind))

figure(114)
imagesc(gamma_inc,D_inc,Si_dist_g_d_p)
title('Probability of trajectories in range.')
xlabel('$\gamma$','Interpreter','latex')
ylabel('$D_{post}$','Interpreter','latex')
colormap(jet)
colorbar

[Dp_ind,gammap_ind] = find(Si_dist_g_d_p==max(max(Si_dist_g_d_p(:,2:end))));
fprintf('Max. Prob. D:%f, gamma:%f\n',D_inc(Dp_ind),gamma_inc(gammap_ind))

% Threshold for Prob. 70%
thres = 0.7;
[Dpt_ind,gammapt_ind] = find(Si_dist_g_d_p(:,2:end) > thres);
[Dmt_ind,gammamt_ind] = find(Si_dist_g_d_m==min(min(Si_dist_g_d_m(Dpt_ind,gammapt_ind))));
fprintf('Min. Dist. Thres. D:%f, gamma:%f thres:%f\n',D_inc(Dmt_ind),gamma_inc(gammamt_ind),thres)


%{
figure(113);
clf;
plot(gamma_inc,Si_diff_m); hold on;
errorbar(gamma_inc,Si_diff_m,Si_diff_e,'k'); hold on;
ylabel('Average Distance |CML-X|')
xlabel('$\gamma$','Interpreter','latex')
xlim([min(gamma_inc) max(gamma_inc)])
%saveas(gcf,strcat('Ave_distance_cml_x_gamma',num2str(gamma),'_',num2str(tri),'.png'))
fprintf('min gamma: %f',gamma_inc(find(Si_diff_m == min(Si_diff_m))))
%}
%%
% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)
global gamma c1 c2 c3 c4 c5 alpha D Dp lam1 lam2 ti tj Dc Dd m ini_max ini_min dt2 ti2 tj2 ti

fCML.p1 = -alpha*c1*c2*c3*c4*c5;
fCML.p2 = alpha*(c1*c3*c4*c5 + c2*c3*c4*c5 + c1*c2*(c4*c5 + c3*(c4 + c5)));
fCML.p3 = -alpha*(c1*c4*c5 + c2*c4*c5 + c3*c4*c5 + c1*c3*(c4 + c5) +  c2*c3*(c4 + c5) + c1*c2*(c3 + c4 + c5));
fCML.p4 = alpha*(c3*c4 + c3*c5 + c4*c5 + c2*(c3 + c4 + c5) +  c1*(c2 + c3 + c4 + c5));
fCML.p5 = -alpha*(c1 + c2 + c3 + c4 + c5); 
fCML.p6 = alpha;

LUp = @(x,myf)1*(myf.p1 + myf.p2*x + myf.p3*x.^2 + myf.p4*x.^3 + myf.p5*x.^4 + myf.p6*x.^5);
% TKI OFF-ON
f = (LUp(x,fCML)-gamma*(heaviside(t-ti(1)))*LUp(x,fCML))*u + (D +(Dp-D)*heaviside(t-ti(1)))*DuDx;
% TET ON-OFF
%f = (LUp(x,fCML) - gamma*exp(-lam1*(t-ti(2))*heaviside(t-ti(2)))*(heaviside(t-ti(1)))*LUp(x,fCML) )*u +D*DuDx;
%f = (LUp(x,fCML) - F(t,x))*u +D*DuDx;
c = 1;
s = 0;
end
% --------------------------------------------------------------
function u0 = pdex1ic(x)
global gamma c1 c2 c3 c4 c5 alpha D lam1 lam2 ti tj Dc Dd m ini_max ini_min dt2 ti2 tj2 ti

%dc=0.1;
dc=abs(c1-c2)/2;
if (x < c1+dc && x > c1-dc)
     u0 = 0.01;   
else
     u0= 0;
end

end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = 0;
qr = ur;
end