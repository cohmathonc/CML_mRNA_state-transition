
%% parameter setting
clear
%gpuDevice
global gamma c1 c2 c3 c4 c5 alpha D Dp lam1 lam2 ti tj Dc Dd m ini_max ini_min dt2 ti2 tj2 ti
sympref('HeavisideAtOrigin',1);

m = 0;

t0 = 0;
%tf = 18;
tf = 13;
c1 = 138.55368;
c2 = 111.85604;
c3 = -43.88014;
c4 = -125.45624;
c5 = -208.51554; 
ti = [6 9];
% tki half-life 17 hr=0.1 week
lam1 = log(2)/0.1; 
%t = linspace(t0,tf,tf);
t = linspace(t0,tf,10);

alpha_inc = 0.5:0.1:2;
alpha_inc = 1.5;

alpha_inc = 100.0;

%D = 0.086243543;
%CML = 0.086243543
%Ctrl = 0.002509334
%TOTO (pre-RX) = 0.087572573
%TKI (all time points) = 0.036646070
%TKI (pre-RX) = 0.072635

A=1;
B=0;
mn=(B-A)/(c5-c1);
Y = @(x)(mn.*(x-c1))+A;
c1 = Y(c1); c2 = Y(c2); c3 = Y(c3); c4=Y(c4); c5=Y(c5);
x = linspace(-0.2,1.2);
dx = x(2)-x(1);

% for SDE trajectories
Ntrials = 10;
dt = 1/(24*7);
nPeriods = round(18/(dt));

% candidate
%gamma_inc = 2;
%D_inc = 0.001;
% min.dist
%gamma_inc = 0;
D_inc = 0.005;


cmldata = readtable("CML.mRNA_state-space_sampleInfo_table.tsv", "FileType","text",'Delimiter','\t');
cmldata_cml_space = Y(cmldata.CML_space);

%'TET_OFF_B', 'TET_OFF_NIL_ON_D', 'TET_OFF_ON_A', 'TET_ON_C'
tet_off_on_a = find(contains(cmldata.treatment,'TET_OFF_ON_A'));
tet_off_b = find(contains(cmldata.treatment,'TET_OFF_B'));
tet_on_c = find(contains(cmldata.treatment,'TET_ON_C'));
tet_off_nil_on_d = find(contains(cmldata.treatment,'TET_OFF_NIL_ON_D'));

% remove outlier, mid = 487
%tet_off_b = tet_off_b(cmldata.mouse_id(tet_off_b) ~= 487);

%mid
tet_off_b_mid = unique(cmldata.mouse_id(tet_off_b));
tet_on_c_mid = unique(cmldata.mouse_id(tet_on_c));
tet_off_nil_on_d_mid = unique(cmldata.mouse_id(tet_off_nil_on_d));
tet_off_on_a_mid = unique(cmldata.mouse_id(tet_off_on_a));

% TKI ON OFF
cml_treatment_t = cmldata.timepoint(tet_off_b);
cml_treatment_s = cmldata_cml_space(tet_off_b);

cml_treatment_t_timepoint = unique(cml_treatment_t);
cml_treatment_t_timepoint = cml_treatment_t_timepoint(cml_treatment_t_timepoint > 5);

cml_treatment_ini = cml_treatment_s(find(cml_treatment_t == 0));
cml_treatment_ini(4) = [];

for alpha=alpha_inc 
    
   
%for alpha=alpha_inc
%    for D=D_inc
    
%% potential function
fCML.p1 = -alpha*c1*c2*c3*c4*c5;
fCML.p2 = alpha*(c1*c3*c4*c5 + c2*c3*c4*c5 + c1*c2*(c4*c5 + c3*(c4 + c5)));
fCML.p3 = -alpha*(c1*c4*c5 + c2*c4*c5 + c3*c4*c5 + c1*c3*(c4 + c5) +  c2*c3*(c4 + c5) + c1*c2*(c3 + c4 + c5));
fCML.p4 = alpha*(c3*c4 + c3*c5 + c4*c5 + c2*(c3 + c4 + c5) +  c1*(c2 + c3 + c4 + c5));
fCML.p5 = -alpha*(c1 + c2 + c3 + c4 + c5); 
fCML.p6 = alpha;

fCML.c1 = c1;
fCML.c2 = c2;
fCML.c3 = c3;
fCML.c4 = c4;
fCML.c5 = c5; 
fCML.a = alpha;

% original
%LUp = @(x,myf)1*(myf.p1 + myf.p2*x + myf.p3*x.^2 + myf.p4*x.^3 + myf.p5*x.^4 + myf.p6*x.^5);
%LU = @(x,myf)1*(myf.p1*x + (myf.p2/2)*x.^2 + (myf.p3/3)*x.^3 + (myf.p4/4)*x.^4 + (myf.p5/5)*x.^5 + (myf.p6/6)*x.^6);

% circuit
LUp = @(x,myf) 1/2*myf.a*myf.c1*myf.c2*(-myf.c3 - myf.c4)*myf.c4*myf.c5 + ...
 myf.a*(-(1/2)*myf.c1*myf.c2*(-myf.c3 - myf.c4)*myf.c4 - 1/2*myf.c1*myf.c2*(-myf.c3 - myf.c4)*myf.c5 + ...
    myf.c1*myf.c2*myf.c4*myf.c5 - 1/2*myf.c1*(-myf.c3 - myf.c4)*myf.c4*myf.c5 - ...
    1/2*myf.c2*(-myf.c3 - myf.c4)*myf.c4*myf.c5)*x + ...
 myf.a*(1/2*myf.c1*myf.c2*(-myf.c3 - myf.c4) - myf.c1*myf.c2*myf.c4 + 1/2*myf.c1*(-myf.c3 - myf.c4)*myf.c4 + ...
    1/2*myf.c2*(-myf.c3 - myf.c4)*myf.c4 - myf.c1*myf.c2*myf.c5 + 1/2*myf.c1*(-myf.c3 - myf.c4)*myf.c5 + ...
    1/2*myf.c2*(-myf.c3 - myf.c4)*myf.c5 - myf.c1*myf.c4*myf.c5 - myf.c2*myf.c4*myf.c5 + ...
    1/2*(-myf.c3 - myf.c4)*myf.c4*myf.c5)*x.^2 + ...
 myf.a*(myf.c1*myf.c2 - 1/2*myf.c1*(-myf.c3 - myf.c4) - 1/2*myf.c2*(-myf.c3 - myf.c4) + myf.c1*myf.c4 + myf.c2*myf.c4 + ...
    1/2*myf.c4*(myf.c3 + myf.c4) + myf.c1*myf.c5 + myf.c2*myf.c5 + myf.c4*myf.c5 + 1/2*(myf.c3 + myf.c4)*myf.c5)*x.^3 +...
  myf.a*(-myf.c1 - myf.c2 + 1/2*(-myf.c3 - myf.c4) - myf.c4 - myf.c5)*x.^4 + myf.a*x.^5;
  
  
LU = @(x,myf) -(1/2)*myf.a*myf.c1*myf.c2*myf.c4*(myf.c3 + myf.c4)*myf.c5*x + ...
 1/4*myf.a*(myf.c2*myf.c4*(myf.c3 + myf.c4)*myf.c5 + ...
    myf.c1*(myf.c4*(myf.c3 + myf.c4)*myf.c5 + myf.c2*myf.c3*(myf.c4 + myf.c5) + myf.c2*myf.c4*(myf.c4 + 3*myf.c5)))*x.^2 - ...
 1/6*myf.a*(myf.c4*(myf.c3 + myf.c4)*myf.c5 + myf.c2*myf.c3*(myf.c4 + myf.c5) + myf.c2*myf.c4*(myf.c4 + 3*myf.c5) + ...
    myf.c1*(myf.c3*(myf.c4 + myf.c5) + myf.c2*(myf.c3 + 3*myf.c4 + 2*myf.c5) + myf.c4*(myf.c4 + 3*myf.c5)))*x.^3 + ...
 1/8*myf.a*(myf.c3*myf.c4 + myf.c4^2 + myf.c3*myf.c5 + 3*myf.c4*myf.c5 + myf.c2*(myf.c3 + 3*myf.c4 + 2*myf.c5) + ...
    myf.c1*(2*myf.c2 + myf.c3 + 3*myf.c4 + 2*myf.c5))*x.^4 - ...
 1/10*myf.a*(2*myf.c1 + 2*myf.c2 + myf.c3 + 3*myf.c4 + 2*myf.c5)*x.^5 + (myf.a*x.^6)/6;
 


Si_diff = [];
Si_diff_m = [];
Si_diff_e = [];

g_inc = 1;
d_inc = 1;
for D=D_inc

%% FP solution

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
% Extract the first solution component as u.
u = sol(:,:,1);
% scale to be PDF
u = u.*(1./sum(u,2));
%% plot all time solutions
%{
figure(100); 
subplot(1,3,1)
%plot(x,LU(x,fCML)-gamma*LU(x,fCML));
plot(x,LU(x,fCML));
ylabel('$U_p$','Interpreter','latex')
set(gca,'xdir','reverse','fontsize',14,'xtick',[c5 c4 c3 c2 c1 ],'xticklabel',{'c_5','c_4','c_3','c_2','c_1'})
subplot(1,3,2)
plot(x,u(1,:),'b')
hold on;
plot(x,u(end,:),'r')
ylabel('$P(x,t)$','Interpreter','latex')
set(gca,'xdir','reverse','fontsize',14,'xtick',[c5 c4 c3 c2 c1 ],'xticklabel',{'c_5','c_4','c_3','c_2','c_1'})
subplot(1,3,3)
plot(x,exp(-(LU(x,fCML))))
ylabel('$\exp(-(U_p) )$','Interpreter','latex')
set(gca,'xdir','reverse','fontsize',14,'xtick',[c5 c4 c3 c2 c1 ],'xticklabel',{'c_5','c_4','c_3','c_2','c_1'})
%saveas(gcf,strcat('CML_CML_FP_Up_iesolutions_eUp',num2str(gamma),'.png'))
%}
%% all time solutions
figure(110)
plot(x,u')
xlabel('x'); ylabel('P')
set(gca,'xdir','reverse','fontsize',14,'xtick',[c5 c4 c3 c2 c1 ],'xticklabel',{'c_5','c_4','c_3','c_2','c_1'})
%saveas(gcf,strcat('CML_CML_FP_all_time_solutions',num2str(gamma),'.png'))
%% surface plot
figure(90)
h = surfc(x,t(2:end),u(2:end,:));
%h = surfc(x,t(2:end),u(2:end,:),'EdgeAlpha',0,'EdgeColor','interp');
%shading interp
c = h(2);
c.ContourZLevel = -.05;
c.LevelList = [0.0100 0.012500 0.01500 0.01750 0.0200 0.0400 0.0600 0.0800 0.1000 0.1200 0.1400];
view(-150,21)
%view(-150,41)
title('CML')
xlabel('CML state')
zlabel('Probability Density')
ylabel('Time (Weeks)')
pbaspect([1 1 1])
zlim([-0.05 0.16])
set(gca,'fontsize',18,'xtick',[c5 c4 c3 c2 c1 ],'xticklabel',{'c_5','c_4','c_3','c_2','c_1'},'ytick',0:3:tf,'ztick',-0.05:0.05:0.15)
saveas(gcf,strcat('CML_FPSurface.png'))
saveas(gcf,strcat('CML_FPSurface.fig'))
%}
%% check trajectories
% adjust beta value

G = @(t,X) sqrt( 2*( D ) );
F = @(t,X) -LUp(X,fCML) ;
Si = [];
Si_s = [];
s=3;
for ini=1:size(cml_treatment_ini,1)
langevin = sde(F, G,'StartState',cml_treatment_ini(ini));
[Sii,T] = simulate(langevin, nPeriods, 'DeltaTime', dt, 'nTrials', Ntrials);
Sii = Sii(:,all(~isnan(Sii)));
Sii = Sii(:,sum(Sii <= 1.4 & Sii > -0.4) == size(Sii,1));
Si = [Si Sii];
Si_s = [Si_s Sii(:,1:s)];
end

%save('SolveFP_CML_survival_Si.mat','Si')
%save('SolveFP_CML_survival_Si_s.mat','Si_s')


%[d,ix] = min(abs(T-6));
%Si = Si(:,Si(ix,:) < (c4 + c5)/2);

Si_dist_g_d_p(d_inc,g_inc) = size(Si,2)/Ntrials;

figure(101);
clf;
[d,ix] = min(abs(T-6));
[de,ixe] = min(abs(T-9));
%plot(T,Si(:,Si(ix,:) < (c4 + c5)/2),'Color',[.6 .6 .6]); hold on;
%plot(T,Si(:,(Si(ix,:) < (c4 + c5)/2) & (Si(ixe,:) > (c4 + c5)/2)),'Color',[.6 .6 .6]); hold on;
%plot(T,Si,'Color',[.8 .8 .8]); hold on;
%plot(T,Si,'Color',[.6 .6 .6]); hold on;
plot(T,Si_s,'Color',[.6 .6 .6]); hold on;
%plot(T,Si_s(:,Si_s(ixe,:) < c4),'Color',[.6 .6 .6]); hold on;

for i =1:size(tet_off_b_mid,1)    
    ind = find(cmldata.mouse_id==tet_off_b_mid(i));
    plot(cmldata.timepoint(ind),cmldata_cml_space(ind),'r'); hold on;
end
%scatter(cmldata.timepoint(tet_off_b),cmldata_cml_space(tet_off_b),'r'); hold on;

for i =1:size(tet_on_c_mid,1)
    ind = find(cmldata.mouse_id==tet_on_c_mid(i));
    plot(cmldata.timepoint(ind),cmldata_cml_space(ind),'k'); hold on;
end
%scatter(cmldata.timepoint(tet_on_c),cmldata_cml_space(tet_on_c),'k'); hold on;

yline(c1,':k','fontweight','bold','fontsize',16,'LineWidth',2); hold on;
yline(c2,':k','fontweight','bold','fontsize',16,'LineWidth',2); hold on;
yline(c3,':k','fontweight','bold','fontsize',16,'LineWidth',2); hold on;
yline(c4,':k','fontweight','bold','fontsize',16,'LineWidth',2); hold on;
yline(c5,':k','fontweight','bold','fontsize',16,'LineWidth',2); hold on;
%xline(ti(1),'k','fontweight','bold','fontsize',16,'LineWidth',1.5); hold on;
%xline(ti(2),'k','fontweight','bold','fontsize',16,'LineWidth',1.5); hold on;
xlim([0 18])
title('CML')
%ylabel('CML state')
yl = ylabel('CML state');
yl.Position(2) = yl.Position(2);
xlabel('Time (Weeks)')
text(-1,c1,'$c_1$','Interpreter', 'latex','fontsize',14);
text(-1,c2,'$c_2$','Interpreter', 'latex','fontsize',14);
text(-1,c3,'$c_3$','Interpreter', 'latex','fontsize',14);
text(-1,c4,'$c_4$','Interpreter', 'latex','fontsize',14);
text(-1,c5,'$c_5$','Interpreter', 'latex','fontsize',14);
set(gca,'YTick',[]);
ylabel('CML state')
xlabel('Time (Weeks)')
fname = strcat('CML_CML_FP_trajectories_a',num2str(alpha_inc),'_d',num2str(D_inc),'.png');
%saveas(gcf,fname)

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
%fprintf('P(Si):%f, D:%f, dist:%f\n',size(Si,2)/Ntrials,D,mean2(Si_diff_i))


%% calculate survival
load("SolveFP_CML_survival_Si.mat")

sv = [];
timepointsv = 0:3:tf;
for i=timepointsv
    [d,ix] = min(abs(T-i));
    if i==0
        svi = 1;
    else
        %svi = 1 - sum(sum(u(1:i,x<=c5)))/sum(sum(u(1:i,:)));
        svi = 1- size(find(sum(Si(1:ix,:) <= c5) ~=0 ),2)/size(Si,2);
    end
    sv = [sv svi];
end


% dob vs c5
status = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';
time =  [84.04167 84.04167 84.04167 84.04167 84.04167 84.04167 84.04167 84.04167 84.04167 84.04167 84.04167 84.04167 84.04167 84.04167 84.04167 84.04167 84.04167 84.04167 84.04167...
58.04167 58.04167 58.04167 58.04167 58.04167 58.04167 58.04167 58.04167 58.04167 58.04167 58.04167 58.04167 32.00000 32.00000 32.00000 32.00000 32.00000 32.00000 32.00000...
32.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000...
87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000 87.00000]';

m_suv = [time status];

figure(118)
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 3], 'PaperUnits', 'Inches', 'PaperSize', [5, 3])
%ecdf(m_suv(:,1),'Censoring',~m_suv(:,2),'function','survivor','-k');
[f1,x1]=ecdf(m_suv(:,1),'Censoring',~m_suv(:,2),'function','survivor'); 
x1 = x1/7;
x1 = [(1:3)'; x1 ];
%f1e = f1(end-1);
%f1(end) = f1e;
f1 = [ones(3,1); f1];
stairs(x1,f1,'-r','LineWidth',3); hold on;
stairs(timepointsv,sv,'-','Color',[.6 .6 .6],'LineWidth',3)
ylim([0 1])
xlim([0 tf])
title('CML survival')
xlabel('Weeks')
ylabel('Proportion')
%xlabel('Time (weeks)')
%ylabel('$S(x)$','Interpreter','latex')

%[h,p] = ttest2(f1,sv);
[h,p] = kstest2(f1,sv);
fprintf('p values:%f, alpha:%f D:%f dist:%f\n',p,alpha,D,mean2(Si_diff_i))   
%saveas(gcf,strcat('CML_survival_a',num2str(alpha_inc),'_d',num2str(D_inc),'_p',num2str(p),'.png'));
%saveas(gcf,strcat('CML_survival_a',num2str(alpha_inc),'_d',num2str(D_inc),'_p',num2str(p),'.fig'));

end

end
%% movie for trajectories + potential
%{
dir_path = '/home/luechi/projects/CML/Up_t_fig/';
[d,ix] = min(abs(T-15));
%FU = @(t,x) (LU(x,fCML) - gamma*(heaviside(t-ti(1)))*LU(x,fCML));
FU = @(t,x) (LU(x,fCML) - gamma*exp(-lam1*(t-ti(2))*heaviside(t-ti(2)))*(heaviside(t-ti(1)))*LF(x,fCML));

for i = 1:ix
       i
    clf;
    f = figure('visible','off');
    t = tiledlayout(1,2);

    nexttile

    xline(c1,':k','fontweight','bold','fontsize',16,'LineWidth',2); hold on;
    xline(c2,':k','fontweight','bold','fontsize',16,'LineWidth',2); hold on;
    xline(c3,':k','fontweight','bold','fontsize',16,'LineWidth',2); hold on;
    xline(c4,':k','fontweight','bold','fontsize',16,'LineWidth',2); hold on;
    xline(c5,':k','fontweight','bold','fontsize',16,'LineWidth',2); hold on;

    p_ref = plot(x,LU(x,fCML),'r','LineWidth',10); hold on;
    p_ref.Color(4) = 0.15;
    
    plot(x,FU(T(i),x),'r','LineWidth',4); hold on;
    
    for j=1:size(Si,2)
        if (Si(i,j) <= 1.4 & Si(i,j) > -0.4) & (~isnan(Si(i,j))) == 1
            scatter(Si(i,j),FU(T(i),Si(i,j)),50,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','k','MarkerFaceAlpha',1); hold on;
        end
    end
    
    
    xlim([min(x) max(x)]); hold on;
  
    xlabel('CML state','fontweight','bold','fontsize',16); hold on;
    ylabel('$\mathbf{U_p} - \bf{\vec{F}}$','Interpreter','latex','fontweight','bold','fontsize',16); hold on;
    set(gca,'xdir','reverse','fontsize',14,'xtick',[c5 c4 c3 c2 c1 ],'xticklabel',{'c_5','c_4','c_3','c_2','c_1'})
    pbaspect([1 1 1])
   
    nexttile

    %plot(T(1:i),Si(1:i,:),'LineWidth',1.0,'Color',[.75 .75 .75]); hold on

    for j=1:size(Si,2)
        if (Si(i,j) <= 1.4 & Si(i,j) > -0.4) & (~isnan(Si(i,j))) == 1
            plot(T(1:i),Si(1:i,j),'LineWidth',1.0,'Color',[.75 .75 .75]); hold on;
            scatter(T(i),Si(i,j),10,'MarkerFaceColor',[.75 .75 .75],'MarkerEdgeColor','k','MarkerFaceAlpha',1); hold on;
        end
    end

    yline(c1,':k','fontweight','bold','fontsize',16,'LineWidth',2); hold on;
    yline(c2,':k','fontweight','bold','fontsize',16,'LineWidth',2); hold on;
    yline(c3,':k','fontweight','bold','fontsize',16,'LineWidth',2); hold on;
    yline(c4,':k','fontweight','bold','fontsize',16,'LineWidth',2); hold on;
    yline(c5,':k','fontweight','bold','fontsize',16,'LineWidth',2); hold on;
    xline(ti(1),'k','fontweight','bold','fontsize',16,'LineWidth',1.5); hold on;
    xline(ti(2),'k','fontweight','bold','fontsize',16,'LineWidth',1.5); hold on;
    xlim([0 18])
    ylabel('CML state')
    xlabel('Time (Weeks)')

    xlim([0 T(ix)]); hold on;
    ylim([min(x) max(x)]); hold on;
    xlabel('Time'); hold on;
    ylabel('CML state'); hold on;

    set(gca,'fontweight','bold','fontsize',16,'ytick',[c5 c4 c3 c2 c1 ],'yticklabel',{'c_5','c_4','c_3','c_2','c_1'})
    pbaspect([1 1 1])
    set(f, 'PaperUnits', 'inches');
    x_width=10; y_width=5;
    set(f, 'PaperPosition', [0 0 x_width y_width]); %

    saveas(f,strcat(dir_path,'Up_t_',num2str(i),'.png'))
      
end
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

fCML.c1 = c1;
fCML.c2 = c2;
fCML.c3 = c3;
fCML.c4 = c4;
fCML.c5 = c5; 
fCML.a = alpha;

% Original
%LUp = @(x,myf)1*(myf.p1 + myf.p2*x + myf.p3*x.^2 + myf.p4*x.^3 + myf.p5*x.^4 + myf.p6*x.^5);

% 3 critical points
%LUp = @(x,myf) -myf.a*myf.c1*myf.c2*myf.c5 + myf.a*myf.c1*myf.c2*x + myf.a*myf.c1*myf.c5*x + myf.a*myf.c2*myf.c5*x - myf.a*myf.c1*x.^2 - ...
% myf.a*myf.c2*x.^2 - myf.a*myf.c5*x.^2 + myf.a*x.^3;

LUp = @(x,myf) 1/2*myf.a*myf.c1*myf.c2*(-myf.c3 - myf.c4)*myf.c4*myf.c5 + ...
 myf.a*(-(1/2)*myf.c1*myf.c2*(-myf.c3 - myf.c4)*myf.c4 - 1/2*myf.c1*myf.c2*(-myf.c3 - myf.c4)*myf.c5 + ...
    myf.c1*myf.c2*myf.c4*myf.c5 - 1/2*myf.c1*(-myf.c3 - myf.c4)*myf.c4*myf.c5 - ...
    1/2*myf.c2*(-myf.c3 - myf.c4)*myf.c4*myf.c5)*x + ...
 myf.a*(1/2*myf.c1*myf.c2*(-myf.c3 - myf.c4) - myf.c1*myf.c2*myf.c4 + 1/2*myf.c1*(-myf.c3 - myf.c4)*myf.c4 + ...
    1/2*myf.c2*(-myf.c3 - myf.c4)*myf.c4 - myf.c1*myf.c2*myf.c5 + 1/2*myf.c1*(-myf.c3 - myf.c4)*myf.c5 + ...
    1/2*myf.c2*(-myf.c3 - myf.c4)*myf.c5 - myf.c1*myf.c4*myf.c5 - myf.c2*myf.c4*myf.c5 + ...
    1/2*(-myf.c3 - myf.c4)*myf.c4*myf.c5)*x.^2 + ...
 myf.a*(myf.c1*myf.c2 - 1/2*myf.c1*(-myf.c3 - myf.c4) - 1/2*myf.c2*(-myf.c3 - myf.c4) + myf.c1*myf.c4 + myf.c2*myf.c4 + ...
    1/2*myf.c4*(myf.c3 + myf.c4) + myf.c1*myf.c5 + myf.c2*myf.c5 + myf.c4*myf.c5 + 1/2*(myf.c3 + myf.c4)*myf.c5)*x.^3 +...
  myf.a*(-myf.c1 - myf.c2 + 1/2*(-myf.c3 - myf.c4) - myf.c4 - myf.c5)*x.^4 + myf.a*x.^5; 

f = ( LUp(x,fCML) )*u + D*DuDx;

c = 1;
s = 0;
end
% --------------------------------------------------------------
function u0 = pdex1ic(x)
global gamma c1 c2 c3 c4 c5 alpha D lam1 lam2 ti tj Dc Dd m ini_max ini_min dt2 ti2 tj2 ti
%{
c1o = 138.55368;
c5o = -208.51554; 

A=1;
B=0;
mn=(B-A)/(c5o-c1o);
Y = @(x)(mn.*(x-c1o))+A;

cmldata = readtable("CML.mRNA_state-space_sampleInfo_table.tsv", "FileType","text",'Delimiter','\t');
cmldata_cml_space = Y(cmldata.CML_space);

%'TET_OFF_B', 'TET_OFF_NIL_ON_D', 'TET_OFF_ON_A', 'TET_ON_C'
tet_off_b = find(contains(cmldata.treatment,'TET_OFF_B'));
% TKI ON OFF
cml_treatment_t = cmldata.timepoint(tet_off_b);
cml_treatment_s = cmldata_cml_space(tet_off_b);
cml_treatment_ini = cml_treatment_s(find(cml_treatment_t == 0));
pd = fitdist(cml_treatment_ini,'Kernel');
u0 = pdf(pd,x);
%}
%plot(cml_treatment_ini,cml_treatment_ini,'r'); hold on;
%plot(x,pdf(pd,x),'b'); hold on;
%hist(cml_treatment_ini)

dc=abs(c1-c2)/2;
if (x < c1+dc && x > c1-dc)
     %u0 = 0.01;
     u0 = 1/(2*dc);
else
     u0= 0;
end
%}
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = 0;
qr = ur;
end