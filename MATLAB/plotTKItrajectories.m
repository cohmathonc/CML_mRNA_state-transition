LoadCMLdata
t0 = 0;

t0tki = 6; 
tftki = 9;
tkidt = tftki - t0tki;


t0post = 9;
tfpost = 16;
postdt = tfpost - t0post;


x0 = [-129 -108  -66 19]'; %remove unstable point c2
%remap [-129 19] interval to [0 1] with linear transformation
%[-129, 0], [19, 1]
A = 1;
B = 0; 

m = (B-A) / (x0(end)-x0(1));
cps = m.*(x0-x0(1)) + A;

%group = 'TET_OFF_NIL_ON_D';
%group = 'TET_OFF_B';
group = 'TET_ON_C';
%group = 'TET_OFF_ON_A';

rows = find(cmldata.treatment==group); %samples in TKI group

ids = unique(cmldata(rows,"mouse_id")); %mice ids in TKI group

N = height(ids); %number of unique mice

figure
hold on; box on;

for i = 1:N

    t = cmldata.timepoint(cmldata.mouse_id==ids.mouse_id(i));
    pc1 = cmldata.CML_space(cmldata.mouse_id==ids.mouse_id(i));
    spc1 = m.*(pc1-x0(1)) + A;
mean(spc1)
    plot(t,spc1)

end

plot([0 tfpost],[cps(1) cps(1)],'k');
plot([0 tfpost],[cps(2) cps(2)],'k');
plot([0 tfpost],[cps(3) cps(3)],'k');
plot([0 tfpost],[cps(4) cps(4)],'k');

if contains(group,'D')
    plot([t0tki t0tki],[0 1],':k')
    plot([tftki tftki],[0 1],':k')
end

if contains(group,'A')
        plot([t0tki t0tki],[0 1],':k')
end

title(['Data: ' strrep(group,'_',' ') ])
set(gca,'ylim',[-.25 1.25],'xlim',[t0 tfpost],'ytick',sort(cps),'yticklabel',{'healthy c_1','c_3','c_4','CML c_5'})