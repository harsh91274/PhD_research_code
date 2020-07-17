%%
%master1
% 1 = Incident ID
% 2 = radius type 1
% 3 = radius type 2
% 4 = radius 1 (m)
% 5 = radius 2 (m)
% 6 = x1
% 7 = y1
% 8 = z1
% 9 = x2
% 10 = y2
% 11 = z2
% 12 = Normal Stress
% 13 = Shear Stress
% 14 = Shear Max
% 15 = cycle
% 16 = peak force tensile or shear
% 17 = area of failure contact
% 18 = energy release due to failure
% 19 = magnitude

% positions(:,1)=master1(:,15);
% positions(:,3)=master1(:,6);
% positions(:,4)=master1(:,7);
% positions(:,5)=master1(:,9);
% positions(:,6)=master1(:,10);

%master2
% 1 = cycle
% 2 = no of events
% 3 = cumulative events
% 4 = fractional events
% 5 = cycle events energy
% 6 = cumulative energy
% 7 = fractional energy
% 8 = S3
% 9 = S1
% 10 = Differential/2
% 11 = Mean stress
% 12 = Movement of platens (m)
% 13 = strain
% 14 = fractional strain
% 15 = AE count/Max AE count
% 16 = cycle energy/ max energy
% 17 = eq magnitude

%bond_master2
% 1 = cycle
% 2 = strain
% 3 = no of total events
% 4 = no shear
% 5 = shear energy
% 6 = no tensile-shear
% 7 = tensile-shear energy
% 8 = no tensile
%9 = tensile energy

%positions2 = positions 1 corrected to initial positions

%newoutput
%1 = cycle
%2 = S1
%5 = S3
%8 = POROSITY

%%
%BOND output energies are in Pa
%Stress outputs from stresspor are in MPa

clc
clear all
close all

path(path,genpath('C:\export_fig\update'));
path(path,genpath('C:\Program Files\Glyph & Cog\XpdfReader-win64'));
path(path,genpath('C:\Program Files\gs\gs9.22\bin'));
path(path,genpath('C:\xpdf\xpdf-tools-win-4.00'));
% path(path,genpath('C:\Users\Harsh\cmu'));
path(path,genpath('C:\xpdf\xpdf-tools-win-4.00\bin32'));
warning('off','all')
% warning('off','all') %check for eps files

% n1=input('Enter First File Number: ');
% n2=input('Enter Last File Number: ');

n1=1;
n2=100;

mast2=[];
master=[];

for i=n1:n2
    if i<10
        file_name=strcat('BOND00',num2str(i),'.OUT');
    elseif i>=10 && i<100
        file_name=strcat('BOND0',num2str(i),'.OUT');
    else
        file_name=strcat('BOND',num2str(i),'.OUT');
    end
    sdata=importdata(file_name);
    
    if isa(sdata, 'struct')==1
        rel_data1=sdata.data;
        rel_data2=sdata.textdata(2:end,:);
        rel_data3=str2double(rel_data2);
        
        locaa=find(isnan(rel_data3(:,1))==0);
        for j=1:length(locaa)
            rel_data1(locaa(j)-1,15)=rel_data3(locaa(j),1);
            rel_data1(locaa(j)-1,16)=rel_data3(locaa(j),2);
        end
        rel_data2(locaa,:)=[];
        locaa2=find(isnan(rel_data1(:,1))==1);
        rel_data1(locaa2,:)=[];
        
        mast=rel_data2(:,1);
        mast2=[mast2 ; mast];
        
        mastaa=rel_data1;
        rel_data1(:,17)=i;
        master=[master; rel_data1];
        
    end
end
mastc=zeros(length(mast2),1);
I1=strmatch('Shear',mast2(:,1));
mastc(I1)=1;
I2=strmatch('Tens-shear',mast2(:,1));
% mastc(I2)=2;
mastc(I2)=1;        %overwriting tensile-shear to shear
I3=strmatch('Tensile',mast2(:,1));
mastc(I3)=3;

% master(:,14)=abs(master(:,14));
master1=master(:,[1,2,3,6,7,8,9,10,11,12,13,14,15,16,17]);
% master1(:,14)=abs(master1(:,14));


[sublocs,~] = find(isnan(master1(:,13))==1);
master1(sublocs,13)=0;
master1(sublocs,14)=0;

[p1,q1]=size(master1);

positions=zeros(p1,8);
positions(:,1)=master1(:,15);
positions(:,3)=master1(:,6);
positions(:,4)=master1(:,7);
positions(:,5)=master1(:,9);
positions(:,6)=master1(:,10);

max_Y1=max(positions(:,4));
max_Y2=max(positions(:,6));
pos_Ymat=[max_Y1 max_Y2];
max_Y=max(pos_Ymat);
%Ymin
min_Y1=min(positions(:,4));
min_Y2=min(positions(:,6));
pos_Ymat=[min_Y1 min_Y2];
min_Y=min(pos_Ymat);
%Rmax
r_max1=max(master1(:,4));
r_max2=max(master1(:,5));
max_r=max(r_max1,r_max2);
%Rmin
r_min1=min(master1(:,4));
r_min2=min(master1(:,5));
min_r=min(r_min1,r_min2);

%%
master1(:,q1+4)=0;
cf = input('Enter Youngs Modulus (Pa): ');
cf2=input('Enter Shear Modulus (Pa): ');
tensile_strength=input('Enter Tensile Strength of bonds (Pa): ');
shear_strength=input('Enter Shear Strength of bonds (Pa): ');
e_option=input('Calculate energy using TS and SS? (1)=Yes: ');
vs=input('Enter Shear Wave Velocity of Material (m/s): ');
% pvel=input('Enter Wall Velocity (m/step): ');
pvel=0.8;
% ts=input('Enter Time Step: ');
ts=2*10^-8;
% cycles=input('Enter number of cycles: ');
cycles=2500;

c_limit=input('Enter cycle event clumping limit (2): ');
r_limit=input('Enter distance event clumping limit (0.0012): ');

for i=1:p1
    zave=mean([master1(i,8) master1(i,11)]);    %this is z (m)
    master1(i,17)=0.5*pi*master1(i,4)^2+0.5*pi*master1(i,5)^2;          %this is area (m2)
    area=zave.*(mean([2.*master1(i, 4) 2.*master1(i, 5)]));     %this is volume (m3)
    
    if e_option==1
        if master1(i,13)<=0 && master1(i,14)<=0
            master1(i,16)=abs(master1(i,12));
            master1(i,18)=0.5.*(1/cf).*((tensile_strength).^2).*master1(i,17).*zave; %energy is in joules
        else
            master1(i,16)=abs(master1(i,14));
            master1(i,18)=0.5.*(1/cf2).*((shear_strength).^2).*master1(i,17).*zave; %energy is in joules
        end
             
    else
    
        if master1(i,13)<=0 && master1(i,14)<=0
            master1(i,16)=abs(master1(i,14));
            master1(i,18)=0.5.*(1/cf).*((master1(i,16)./master1(i,17)).^2).*master1(i,17).*zave; %energy is in joules
        else
            master1(i,16)=master1(i,14);
            master1(i,18)=0.5.*(1/cf2).*((master1(i,16)./master1(i,17)).^2).*master1(i,17).*zave; %energy is in joules
        end
    end
    
    %if ym is n/m then use formula below
    %     master1(i,18)=0.5.*(1/cf).*((master1(i,16).^2)./(master1(i,17).^1.5).*master1(i,17); energy is in joules.
    %      master1(i,19)=(log10(master1(i,18))-4.4)/1.5; %compute magnitude using richter scale
    master1(i,19)=(0.67.*log10(master1(i,18)))-2.9;
end
zave1=mean(master1(:,8));
zave2=mean(master1(:,11));
zave=mean([zave1 zave2]);

master1_shear1=master1(I1,:);
master1_ts=master1(I2,:);
master1_shear=[master1_shear1; master1_ts];
master1_tensile=master1(I3,:);
%%
master2=zeros(n2-n1,16);

for i=n1:n2
    master2(i-n1+1,1)= i;
end

cycle_index=master1(:,15);

sum_AE=0;
sum_energy=0;

for i=n1:n2
    [yloc]=find(i==cycle_index);
    e_count=sum(master1(yloc,18));
    
    sum_AE=sum_AE+length(yloc);
    sum_energy=sum_energy+e_count;
    
    master2(i-n1+1,2)=length(yloc);
    master2(i-n1+1,3)=sum_AE;
    
    master2(i-n1+1,5)=e_count;
    master2(i-n1+1,6)=sum_energy;
    
end

for i=n1:n2
    master2(:,4)=master2(:,3)./sum_AE;
    master2(:,7)=master2(:,6)./sum_energy;
end

newoutput=importdata('stresspor2D.out');
dlmwrite('s1.dat',newoutput(:,2),'delimiter','\n');

outfile=zeros(n2-n1+1,4);
%outfile
%1 = strian
%2 = mean stress
%3 = differential stress
%4 = porosity

master2(:,8)=newoutput(2:n2+1,5)./10^6;
master2(:,9)=newoutput(2:n2+1,2)./10^6;

step_disp=ts*pvel*cycles*2;
max_X=newoutput(1,4);
min_X=newoutput(1,3);
s_length=max_X-min_X;

for i=n1:n2
    master2(i,10)=(master2(i,8)+ master2(i,9))/2;
    master2(i,11)=(master2(i,9)- master2(i,8))/2;
    master2(i,12)=master2(i,1).*step_disp;
    master2(i,13)=master2(i,12)./s_length;
end

master2(:,14)=master2(:,12)./master2(end,12);

max_AE=max(master2(:,2));
max_energy=max(master2(:,5));

for i=n1:n2
    master2(i,15)=master2(i,2)./max_AE;
    master2(i,16)=master2(i,5)./max_energy;
    %     master2(i,17)=(log10(master2(i,5))-4.4)/1.5;
    master2(i,17)=(0.67.*log10(master2(i,5)))-2.9;
end

max_S1=max(master2(:,9));
axial_strain=master2(end,13);
smax=max(master2(:,11));
max_eq=max(master1(:,19));
[loc1,~]=find(master2(:,11)==smax);

mean_at_smax=master2(loc1,10);
cycles_at_smax=master2(loc1,1);
strain_at_smax=master2(loc1,13);
s1_at_smax=master2(loc1,9);
s3_at_smax=master2(loc1,8);

[in_loc1,~]=find(master2(:,2)~=0,1,'first');
s1_at_in=master2(in_loc1,9);
strain_at_in=master2(in_loc1,13);
cycles_at_in=master2(in_loc1,1);
mean_at_in=master2(in_loc1,10);
shear_at_in=master2(in_loc1,11);
mean_s3=mean(master2(:,8));

%%
%bond_master2
% 1 = cycle
% 2 = strain
% 3 = no of total events
% 4 = no shear
% 5 = shear energy
% 6 = no tensile-shear
% 7 = tensile-shear energy
% 8 = no tensile
%9 = tensile energy

bond_master2=zeros(length(master2),9);
bond_master2(:,1)=master2(:,1);
bond_master2(:,2)=master2(:,13);

for i=1:length(bond_master2)
    no_bonds=0;
    no_ten=0;
    no_shear=0;
    no_tensh=0;
    e_ten=0;
    e_shear=0;
    e_tensh=0;
    
    [bond_locs,~]=find(bond_master2(i,1)==master1(:,15));
    if isempty(bond_locs)~=1
        no_bonds=length(bond_locs);
        for j=1:length(bond_locs)
            if mastc(bond_locs(j))==1
                no_shear=no_shear+1;
                e_shear=e_shear+master1(bond_locs(j),18);
            elseif mastc(bond_locs(j))==2
                %tensile_shear
                %no_tensh=no_tensh+1;
                %e_tensh=e_tensh+master1(bond_locs(j),18);
                no_shear=no_shear+1;
                e_shear=e_shear+master1(bond_locs(j),18);
            elseif mastc(bond_locs(j))==3
                no_ten=no_ten+1;
                e_ten=e_ten+master1(bond_locs(j),18);
            end
        end
    end
    bond_master2(i,3)=no_bonds;
    bond_master2(i,4)=no_shear;
    bond_master2(i,5)=e_shear;
    bond_master2(i,6)=no_tensh;
    bond_master2(i,7)=e_tensh;
    bond_master2(i,8)=no_ten;
    bond_master2(i,9)=e_ten;
    
end
%master2
% 1 = cycle
% 2 = no of events
% 3 = cumulative events
% 4 = fractional events
% 5 = cycle events energy
% 6 = cumulative energy
% 7 = fractional energy
% 8 = S3
% 9 = S1
% 10 = Differential/2
% 11 = Mean stress
% 12 = Movement of platens (m)
% 13 = strain
% 14 = fractional strain
% 15 = AE count/Max AE count
% 16 = cycle energy/ max energy
% 17 = eq magnitude


n_shear=sum(bond_master2(:,4));
e_shear=sum(bond_master2(:,5));
n_ts=sum(bond_master2(:,6));
e_ts=sum(bond_master2(:,7));
n_tensile=sum(bond_master2(:,8));
e_tensile=sum(bond_master2(:,9));


disp('BOND DATA');
disp(strcat('Number of Shear MF: ',num2str(n_shear)));
disp(strcat('Number of Tensile-Shear MF: ',num2str(n_ts)));
disp(strcat('Number of Tensile MF: ',num2str(n_tensile)));

disp(strcat('Energy from Shear MF: ',num2str(e_shear),' Joules'));
disp(strcat('Energy from Tensile-Shear MF: ',num2str(e_ts),' Joules'));
disp(strcat('Energy from Tensile MF: ',num2str(e_tensile),' Joules'));

%%
%Figure 1
f1_color=linspecer(2);

a_option=input('Autoscale? (1)=yes, (2)=no :');
if a_option==1
    fig1=figure(1);
    f1s1=subplot(2,1,1);
    % %
    [af1,h1,h2]=plotyy(bond_master2(:,2),master2(:,2),bond_master2(:,2),master2(:,9));  %af1(1)=no of event, af1(2)=S1
    % xlimsa=get(af1(1),'xlim');
    xlimsa1=[0 max(master2(:,13))];
    set(af1(1),'xlim',xlimsa1,'ycolor','k');
    set(af1(2),'xlim',xlimsa1,'ycolor','k','ytick',[]);
    hold on
    delete(h1);
    set(h2,'LineWidth',2,'color','b');
    bb1=bar(bond_master2(:,2),[bond_master2(:,4) bond_master2(:,8)],'stack');colormap(f1_color);
    hold on
    set(bb1,'EdgeColor','k','linewidth',0.25);
    hold on
    xlim([0 max(bond_master2(:,2))]);
    set(af1(2),'ytick',[],'ycolor','b');
    % ylim(af1(1),[0 fig1_bondmax]);
    % set(af1(1),'ytick',linspace(0, fig1_bondmax, 6),'ycolor','k');
%     set(af1(1),'ycolor','b');
    % ylim(af1(2),[0 fig1_smax]);
%     set(af1(2),'ytick',[]);
    % set(af1(2),'ytick',linspace(0,fig1_smax,6));
    set(af1(2),'ytickmode','auto');
    hold on
    title(strcat('Total Microfractures: ',num2str(p1)),'fontsize',12);
    xlabel('Strain','fontsize',12);
    ylabel('Number of Microfractures','fontsize',12,'color','k');
    ylabel(af1(2),'Axial Stress(MPa)','fontsize',12,'color','b');
    hold on
    legend('Shear','Tensile','Location','NorthEast');
    hold on
    set(f1s1,'xcolor','k','ycolor','k');
%     set(af1(2),'ycolor','b');
    hold on
    box on
    hold on
    
    %------------------------------------------------------------------------%
    f1s2=subplot(2,1,2);
    
    [af1,h1,h2]=plotyy(bond_master2(:,2),master2(:,5),bond_master2(:,2),master2(:,7));
    % xlimsa=get(af1(1),'xlim');
    xlimsa1=[0 max(master2(:,13))];
    set(af1(1),'xlim',xlimsa1,'ycolor','k');
    set(af1(2),'xlim',xlimsa1,'ycolor','k','ytick',[]);
    hold on
    
    delete(h1);
    set(h2,'LineWidth',2,'color','b');
    
    bb1=bar(bond_master2(:,2),[bond_master2(:,5) bond_master2(:,9)],'stack'); colormap(f1_color);
    hold on
    set(bb1,'EdgeColor','k','linewidth',0.25);
    hold on
    xlim([0 max(bond_master2(:,2))]);
    % ylim(af1(1),[0 fig1_emax]);
    % set(af1(1),'ytick',linspace(0, fig1_emax, 6),'ycolor','k');
    %     set(af1(1),'ycolor','k');
%     set(af1(1),'ytickmode','auto');
 
    set(af1(2),'ytick',[],'ycolor','b');
%     ylim(af1(2),[0 1]);
    set(af1(2),'ytickmode','auto');
    hold on
    title(strcat('Total Energy Released: ',num2str(master2(end,6)),' Joules'),'fontsize',12);
    set(gca,'yticklabel',num2str(get(gca,'ytick')','%.3f'));
    set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'));
    xlabel('Strain','fontsize',12);
    ylabel('Energy(Joules)','fontsize',12,'color','k');
    ylabel(af1(2),'Fractional Energy','fontsize',12,'color','b');
    hold on
    hold on
    legend('Shear','Tensile','Location','East');
    hold on
    set(f1s2,'xcolor','k','ycolor','k');
    hold on
%     set(af1(2),'ycolor','b');
    hold on
    box on
    hold on
    set(gcf,'Color','w');
    if e_option==1
        saveas(gcf,'bond_type_distribution_autoscale_assigned','tiffn');
        saveas(gcf,'bond_type_distribution_autoscale_assigned','epsc');
        saveas(gcf,'bond_type_distribution_autoscale_assigned','fig');
        export_fig fig1 f1_bond_dist_autoscale_assigned -q101 -painters -nocrop -pdf -png -tiff -eps 
    else
        saveas(gcf,'bond_type_distribution_autoscale','tiffn');
        saveas(gcf,'bond_type_distribution_autoscale','epsc');
        saveas(gcf,'bond_type_distribution_autoscale','fig');
        export_fig fig1 f1_bond_dist_autoscale -q101 -painters -nocrop -pdf -png -tiff -eps
    end
    
else
    fig1_bondmax= input('Enter Max Bonds for plot (125): ');
    fig1_emax=input('Enter Max Energy for plot (0.065): ');
    fig1_smax=input('Enter Max S1 for plot (200): ');
    
    fig1=figure(1);
    f1s1=subplot(2,1,1);
    % %
    [af1,h1,h2]=plotyy(bond_master2(:,2),master2(:,2),bond_master2(:,2),master2(:,9));
    % xlimsa=get(af1(1),'xlim');
    xlimsa1=[0 max(master2(:,13))];
    set(af1(1),'xlim',xlimsa1,'ycolor','k');
    set(af1(2),'xlim',xlimsa1,'ycolor','k');
    hold on
    delete(h1);
    set(h2,'LineWidth',2,'color','b');
    bb1=bar(bond_master2(:,2),[bond_master2(:,4) bond_master2(:,8)],'stack');colormap(f1_color);
    hold on
    set(bb1,'EdgeColor','k','linewidth',0.25);
    hold on
    xlim([0 max(bond_master2(:,2))]);
    ylim(af1(1),[0 fig1_bondmax]);
    set(af1(1),'ytick',linspace(0, fig1_bondmax, 5),'ycolor','k');
    set(af1(1),'ycolor','b');
    set(af1(2),'ytick',[]);
    ylim(af1(2),[0 fig1_smax]);
    set(af1(2),'ytick',linspace(0,fig1_smax,5));
    %    set(af1(2),'ytickmode','auto');
    hold on
    title(strcat('Total Microfractures: ',num2str(p1)),'fontsize',12);
    xlabel('Strain','fontsize',12);
    ylabel('Number of Microfractures','fontsize',12,'color','k');
    ylabel(af1(2),'Axial Stress (MPa)','fontsize',12,'color','b');
    hold on
    legend('Shear','Tensile','Location','NorthEast');
    hold on
    set(f1s1,'xcolor','k','ycolor','k');
    set(af1(2),'ycolor','b');
    hold on
    box on
    hold on
    
    %------------------------------------------------------------------------%
    f1s2=subplot(2,1,2);
    
    [af1,h1,h2]=plotyy(bond_master2(:,2),master2(:,5),bond_master2(:,2),master2(:,7));
    %     xlimsa=get(af1(1),'xlim');
    xlimsa1=[0 max(master2(:,13))];
    set(af1(1),'xlim',xlimsa1,'ycolor','k');
    set(af1(2),'xlim',xlimsa1,'ycolor','k');
    hold on
    
    delete(h1);
    set(h2,'LineWidth',2,'color','b');
    
    bb1=bar(bond_master2(:,2),[bond_master2(:,5) bond_master2(:,9)],'stack'); colormap(f1_color);
    hold on
    set(bb1,'EdgeColor','k','linewidth',0.25);
    hold on
    xlim([0 max(bond_master2(:,2))]);
    ylim(af1(1),[0 fig1_emax]);
    set(af1(1),'ytick',linspace(0, fig1_emax, 5),'ycolor','k');
    set(af1(2),'ytick',[]);
    ylim(af1(2),[0 1]);
    set(af1(2),'ytick',linspace(0,1,5));
    hold on
    title(strcat('Total Energy Released: ',num2str(master2(end,6)),' Joules'),'fontsize',12);
    xlabel('Strain','fontsize',12);
    ylabel('Energy(Joules)','fontsize',12,'color','k');
    set(gca,'yticklabel',num2str(get(gca,'ytick')','%.3f'));
    set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'));
    ylabel(af1(2),'Fractional Energy','fontsize',12,'color','b');
    hold on
    %     set(af1(1),'ytickmode','auto');
    hold on
    legend('Shear','Tensile','Location','East');
    hold on
    set(f1s2,'xcolor','k','ycolor','k');
    set(af1(2),'ycolor','b');
    hold on
    box on
    hold on
    
    set(gcf,'Color','w');
    if e_option==1
        saveas(gcf,'bond_type_distribution_assigned','tiffn');
        saveas(gcf,'bond_type_distribution_assigned','epsc');
        saveas(gcf,'bond_type_distribution_assigned','fig');
        export_fig fig1 f1_bond_dist_assigned -q101 -painters -nocrop -pdf -png -tiff -eps 
    else
        saveas(gcf,'bond_type_distribution','tiffn');
        saveas(gcf,'bond_type_distribution','epsc');
        saveas(gcf,'bond_type_distribution','fig');
        export_fig fig1 f1_bond_dist -q101 -painters -nocrop -pdf -png -tiff -eps 
    end
    
end

%%
%correlate cycles to strain for positional correlation

for i=1:p1
    [pos_row,~]=find(positions(i,1)==master2(:,1),1,'first');
    positions(i,2)=master2(pos_row,13);     %record strain
end
positions(:,7)=master1(:,18);

disp('STRESS AND STRAIN RESULTS');
disp(strcat('Number of bonds broken: ', num2str(p1)));
disp(strcat('Max S1: ',num2str(max_S1)));
disp(strcat('Max Axial Strain: ',num2str(axial_strain)));
disp(strcat('Max Shear: ',num2str(smax)));
disp(strcat('Mean Confining Stress (S3): ',num2str(mean_s3)));
disp(strcat('Mean Stress at Max Shear: ',num2str(mean_at_smax)));
disp(strcat('Cycles at Max Shear: ',num2str(cycles_at_smax)));
disp(strcat('Strain at Max Shear: ',num2str(strain_at_smax)));
disp(strcat('S3 at Max Shear: ',num2str(s3_at_smax)));
disp(strcat('Max Earthquake Magnitude (Events not clumped): ',num2str(max_eq)));
disp(strcat('S1 at Inelastic deformation onset: ',num2str(s1_at_in)));
disp(strcat('Cycles at Inelastic deformation onset: ',num2str(cycles_at_in)));
disp(strcat('Strain at Inelastic deformation onset: ',num2str(strain_at_in)));
disp(strcat('Mean Stress at Inelastic deformation onset: ',num2str(mean_at_in)));

%%
%Figure 2
fig2=figure(2);
% subplot(2,1,1)
% plot(master2(:,1),master2(:,9),'r','LineWidth',2)
% hold on
% plot(master2(:,1),master2(:,8),'LineWidth',2)
% hold on
% plot(cycles_at_smax,s1_at_smax,'ro');
% hold on
% % plot(cycles_at_in,s1_at_in,'rx');
% hold on
% axis tight
% xlabel('Cycles','fontweight','bold');ylabel('S1 and S3 (MPa)','fontweight','bold');
% legend('Applied Stress','Confining Stress','Max Stress');
%
% subplot(2,1,2)
plot(master2(:,13),master2(:,9),'r','LineWidth',2);
hold on
plot(master2(:,13),master2(:,8),'LineWidth',2)
hold on
plot(strain_at_smax,s1_at_smax,'ro');
hold on
% plot(strain_at_in,s1_at_in,'rx');
hold on
xlabel('Axial Strain','fontweight','bold');ylabel('S1 and S3','fontweight','bold');
axis tight

satisfied=0;
while satisfied==0
    disp('Select x-range for Bulk Youngs Modulus line, lower limit first');
    [x_range_ym,~]=ginput(2);
    
    [locs1_ym]=find(master2(:,13)>=x_range_ym(1));
    [locs2_ym]=find(master2(:,13)<=x_range_ym(2));
    [locs_ym]=intersect(locs1_ym, locs2_ym);
    
    xfit_ym=master2(locs_ym,13);
    yfit_ym=master2(locs_ym,9);
    ym_fit=polyfit(xfit_ym,yfit_ym,1);
    ym_fit_plot=polyval(ym_fit,xfit_ym);
    ym_plot=plot(xfit_ym,ym_fit_plot,'k--','linewidth',3);
    hold on
    legend('Applied Stress','Confining Stress','Max Stress','best fit','location','best');
    YM_val=ym_fit(1);
    disp(strcat('Bulk Youngs Modulus: ',num2str(YM_val),' MPa'));
    
    
    %     bfit_N_c=N_vector_c(AE_locs_c);
    %     bfit_A_c=A_vector_c(AE_locs_c);
    %     P_c=polyfit(bfit_A_c,bfit_N_c,1);
    %     plotP_c=polyval(P_c,bfit_A_c);
    %     bfp_c=plot(bfit_A_c,plotP_c,'r--','LineWidth',3);
    %     hold on
    %     legend('Magnitude-Frequency relationship','Best fit line','Location','SouthWest');
    %     hold on
    %
    %     b_value_c=-P_c(1);
    %     disp(strcat('b value: ',num2str(b_value_c)));
    
    satisfied=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
    if satisfied==0
        delete(ym_plot);
    end
end

set(gcf,'Color','w');
% set(gcf,'Position',get(0,'Screensize'))
saveas(gcf,'Stress-path','tiffn');
saveas(gcf,'Stress-path','epsc');
export_fig fig2 f2_stress_path -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
fig3=figure(3);
subplot(2,1,1)
plot(master2(:,1),master2(:,10),'r');
hold on
plot(cycles_at_smax,mean_at_smax,'ro');
hold on
% plot(cycles_at_in,mean_at_in,'rx');
hold on
plot(master2(:,1),master2(:,11));
hold on
plot(cycles_at_smax,smax,'bo');
hold on
% plot(cycles_at_in,shear_at_in,'bx');
hold on
xlabel('Cycles');ylabel('Mean and Differential/2 (MPa)');
legend('Mean','Mean@Failure','Differential','Differential/2@Failure');

subplot(2,1,2)
plot(master2(:,13),master2(:,10),'r');
hold on
plot(strain_at_smax,mean_at_smax,'ro');
hold on
% plot(strain_at_in,mean_at_in,'rx');
hold on
plot(master2(:,13),master2(:,11));
hold on
plot(strain_at_smax,smax,'bo');
hold on
% plot(strain_at_in,shear_at_in,'bx');
hold on
axis tight
xlabel('Axial Strain');ylabel('Mean and Differential/2 (MPa)');
legend('Mean','Mean@Failure','Differential','Differential/2@Failure');

set(gcf,'Color','w');
saveas(gcf,'Mean_and_Diff','tiffn');
% saveas(gcf,'Mean_and_Diff','eps');
% export_fig fig3 f3_bond_dist -pdf -png -q101

%%
figure
plot(master2(:,10),master2(:,11))
xlabel('Mean Stress');ylabel('1/2*Differential Stress');
set(gcf,'Color','w');
saveas(gcf,'Mean_vs_Diff','tiffn');
%%
%m-value calculations
sigma_c=input('Enter Unconfined Compressive Strength (MPa): ');
m=((s1_at_smax-s3_at_smax)^2-sigma_c^2)/(s3_at_smax*sigma_c);
disp(strcat('m: ',num2str(m)));
%%
for i=n1:n2
    outfile(i,1)=master2(i,13);     %strain
    outfile(i,2)=master2(i,11);     %mean
    outfile(i,3)=master2(i,10).*2;  %differential
    outfile(i,4)=newoutput(i,8);    %porosity
end

fig15=figure(15);
plot(outfile(:,4),outfile(:,2),'linewidth',2);
hold on;
plot(outfile(:,4),outfile(:,3),'r','linewidth',2);
hold on;
legend('Mean vs phi','Differential vs phi');
xlabel('Porosity'); ylabel('Stress (MPa)');
% xlswrite('mdp.xls',outfile);
dlmwrite('mdp.dat',outfile,'delimiter','\t');
set(gcf,'Color','w');

% saveas(gcf,'MD_vs_phi','tiffn');
% saveas(gcf,'MD_vs_phi','epsc');
saveas(gcf,'MD_vs_phi','fig');
export_fig fig15 f12_MD_phi -q101 -painters -nocrop -pdf -png -tiff -eps 

%%
%Plotting spatial distributions of broken bonds

min_x1=min(min(master1(:,6)));          %find min and max positions
min_x2=min(min(master1(:,9)));
min_y1=min(min(master1(:,7)));
min_y2=min(min(master1(:,10)));
min_x3=min(min_x1,min_x2);
min_y3=min(min_y1,min_y2);

max_x1=max(max(master1(:,6)));          %find min and max positions
max_x2=max(max(master1(:,9)));
max_y1=min(max(master1(:,7)));
max_y2=min(max(master1(:,10)));
max_x3=min(max_x1,max_x2);
max_y3=min(max_y1,max_y2);

x_range3=max_x3-min_x3;
y_range3=max_y3-min_y3;

a_ratio3=x_range3/y_range3;

plot_buff=input('Ready for Spatial Distribution plot ?: ');
fig14=figure(14);
aviobj=avifile('Rock_Damage.avi','compression','None','fps',8);

cm=jet(n2-n1+1);

for i=1:length(master1)     %modified for 5.8a 
    if mastc(i)==1
        positions(i,8)=1;           %1 for shear, 0 for tensile
    end
end


for i=n1:n2
    [pos_x,~]=find(i==positions(:,1));
    axis([min_x3 max_x3 min_y3 max_y3]); pbaspect([a_ratio3 1 1]); box on;
    set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'));
    set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'));
    title(strcat('Strain: ',num2str(master2(i,13))));
    hold on
    
      if pos_x~=0
        pos_pm=positions(pos_x,:);
        pos_spos=find(pos_pm(:,8)==1);
        pos_tpos=find(pos_pm(:,8)==0);
        pos_pm_s=pos_pm(pos_spos,:);
        pos_pm_t=pos_pm(pos_tpos,:);
     
        if pos_tpos~=0      %tensile
            plot(((pos_pm_t(:,3)+pos_pm_t(:,5))./2),((pos_pm_t(:,4)+pos_pm_t(:,6))./2),'.','color',cm(i,:),'linewidth',1);
            hold on
        end
        if pos_spos~=0  %shear
            plot(((pos_pm_s(:,3)+pos_pm_s(:,5))./2),((pos_pm_s(:,4)+pos_pm_s(:,6))./2),'o','color',cm(i,:),'linewidth',0.8);
            hold on
        end
      end
        
%     if pos_x~=0
%         pos_pm=positions(pos_x,:);
%         plot(((pos_pm(:,3)+pos_pm(:,5))./2),((pos_pm(:,4)+pos_pm(:,6))./2),'x','color',cm(i,:),'LineWidth',2);
%         hold on
%     end
    F=getframe(fig14);
    aviobj=addframe(aviobj,F);
end
caxis([0 axial_strain]);
colorbar('SouthOutside');
% cth=get(hcb,'Title');
% ts='Strain';
% set(cth,'String',ts);
F=getframe(fig14);
aviobj=addframe(aviobj,F);
xlabel('Sample length (m)'); ylabel('Sample width (m)');
xlabel('Sample length (m)');ylabel('Sample width (m)');
set(gcf,'Color','w');
saveas(gcf,'Damage_distribution','tiffn');
aviobj=close(aviobj);

%%
%Plot corrected positions

bdata=importdata('S000.OUT');
bdata=bdata(:,1:4);

positions2=zeros(length(master1),8);
%positions 1 and 2
%1 = cycle
%2 = strain
%3 = x1 (initial for positions2)
%4 = y1 (initial for positions2)
%5 = x2 (initial for positions2)
%6 = y2 (initial for positions2)
%7 = energy
%8 = counter (1 is for shear, 0 for tensile)

positions2(:,1)=master1(:,1);
positions2(:,7)=master1(:,18);

for i=1:length(master1)
    positions2(i,1)=master1(i,15);
    positions2(i,7)=master1(i,18);
    
    particle1=master1(i,2);
    particle2=master1(i,3);
    
    [ploc1,~]=find(bdata(:,1)==particle1);
    [ploc2,~]=find(bdata(:,1)==particle2);
    
    positions2(i,3)=bdata(ploc1,2);
    positions2(i,4)=bdata(ploc1,3);
    positions2(i,5)=bdata(ploc2,2);
    positions2(i,6)=bdata(ploc2,3);
    
    if mastc(i)==1
        positions2(i,8)=1;
    end
end

for i=1:p1
    [pos_row,~]=find(positions2(i,1)==master2(:,1),1,'first');
    positions2(i,2)=master2(pos_row,13);     %record strain
end

min_x4=min(min(positions2(:,3)));          %find min and max positions
min_x5=min(min(positions2(:,5)));
min_y4=min(min(positions2(:,4)));
min_y5=min(min(positions2(:,6)));
min_x6=min(min_x4,min_x5);
min_y6=min(min_y4,min_y5);

max_x4=max(max(positions2(:,3)));          %find min and max positions
max_x5=max(max(positions2(:,5)));
max_y4=min(max(positions2(:,4)));
max_y5=min(max(positions2(:,6)));
max_x6=min(max_x4,max_x5);
max_y6=min(max_y4,max_y5);

x_range6=max_x6-min_x6;
y_range6=max_y6-min_y6;

a_ratio6=x_range6/y_range6;

%plotting
plot_buff_c=input('Ready for Corrected Spatial Distribution plot ?: ');
fig=figure(5);
aviobj_c=avifile('Rock_Damage_Corrected.avi','compression','None','fps',8);
mkdir damage_figs

cm=jet(n2-n1+1);
for i=n1:n2
    [pos_x,~]=find(i==positions2(:,1));
    axis([min_x6 max_x6 min_y6 max_y6]); pbaspect([a_ratio6 1 1]); box on;
    set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'));
    set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'));
    %     title({(strcat('Strain: ',num2str(master2(i,13)))),'fontsize',12,'fontweight','bold';'/n . = tensile | o = shear'});
    title(strcat('Strain: ',num2str(master2(i,13))),'fontsize',12,'fontweight','bold');
    hold on
    if pos_x~=0
        pos_pm=positions2(pos_x,:);
        pos_spos=find(pos_pm(:,8)==1);
        pos_tpos=find(pos_pm(:,8)==0);
        pos_pm_s=pos_pm(pos_spos,:);
        pos_pm_t=pos_pm(pos_tpos,:);
        
        %
        %         [pos_l1,~]=size(pos_pm);
        %         pos_s=0;
        %         pos_t=0;
        %         for j=1:pos_l1
        %             if pos_pm(j,8)==1
        %                 pos_s=pos_s+1;
        %                 pos_pm_s(pos_s,:)=pos_pm(:,j);
        %             else
        %                 pos_t=pos_t+1;
        %                 pos_pm_t(pos_t,:)=pos_pm(:,j);
        %             end
        %         end
        if pos_tpos~=0      %tensile
            plot(((pos_pm_t(:,3)+pos_pm_t(:,5))./2),((pos_pm_t(:,4)+pos_pm_t(:,6))./2),'.','color',cm(i,:),'linewidth',1);
            hold on
        end
        if pos_spos~=0  %shear
            plot(((pos_pm_s(:,3)+pos_pm_s(:,5))./2),((pos_pm_s(:,4)+pos_pm_s(:,6))./2),'o','color',cm(i,:),'linewidth',0.8);
            hold on
            
        end
    end
    
    F=getframe(fig);
    aviobj_c=addframe(aviobj_c,F);
    
    if rem(i,2)==0
        fname=strcat(pwd,'\damage_figs\','damage_',num2str(i));
        saveas(gcf,fname,'tiffn');
        saveas(gcf,fname,'epsc');
        saveas(gcf,fname,'fig');
        saveas(gcf, fname, 'png');
        saveas(gcf, fname, 'jpg');
    end
    
end


caxis([0 axial_strain]);
colorbar('SouthOutside');
% title(hcb,'Strain','Location','South');
F=getframe(fig);
aviobj_c=addframe(aviobj_c,F);
xlabel('Sample length (m)','fontweight','bold'); ylabel('Sample width (m)','fontweight','bold');
axis tight
box on
set(gcf,'Color','w');
saveas(gcf,'Damage_distribution_Corrected','tiffn');
saveas(gcf,'Damage_distribution_Corrected','epsc');
saveas(gcf,'Damage_distribution_Corrected','fig');
export_fig F f5_damage_dist_corr -q101 -painters -nocrop -pdf -png -tiff -eps 
aviobj_c=close(aviobj_c);


% axis tight
% % corrnew=copyobj(gcf,0);
% xlabel(gca'');ylabel(gca,'');
% colorbar off
%%
fig6=open('Damage_distribution_Corrected.fig');
% oldf=gcf;
% newf=copyobj(oldf);
title('');
xlabel('');
ylabel('');
colorbar off;
axis tight
set(gcf,'Color','w');
saveas(gcf,'Damage_distribution_Corrected_clean','tiffn');
saveas(gcf,'Damage_distribution_Corrected_clean','epsc');
export_fig fig6 f6_damage_dist_corr_clean -q101 -painters -nocrop -pdf -png -tiff -eps 

%%
disp('ACOUSTIC EMISSION AND BRITTLENESS RESULTS');
figure
subplot(2,1,1)
f1=bar(master2(:,13),master2(:,15));
set(f1,'FaceColor','b','EdgeColor','k')
axis tight
hold on
[af1,h1,h2]=plotyy(master2(:,13),master2(:,4),master2(:,13),master2(:,2));
axis tight
set(af1(2),'XTickLabel',[],'XTick',[]);
delete(h2);
set(h1,'Color','r','LineWidth',2);
set(af1,{'ycolor'},{'k';'k'});
xlabel('Axial Strain');ylabel('Events Ratio');
legend('AE','Cumulative AE','Location','NorthEast','Orientation','horizontal');
title(strcat('Total AE event count: ',num2str(p1)));
hold on

subplot(2,1,2)
f2=bar(master2(:,13),master2(:,16));
set(f2,'FaceColor','g','EdgeColor','k')
axis tight
hold on
% plot(master2(:,13),master2(:,7),'r','LineWidth',2);
[af2,h1,h2]=plotyy(master2(:,13),master2(:,7),master2(:,13),master2(:,5));
axis tight
set(af2(2),'XTickLabel',[],'XTick',[]);
delete(h2);
set(h1,'Color','r','LineWidth',2);
set(af2,{'ycolor'},{'k';'k'});
xlabel('Axial Strain');ylabel('Energy Ratio');
legend('Fractional E ratio','Cumulative Energy','Location','NorthEast','Orientation','horizontal');
title(strcat('Total Energy Released: ',num2str(master2(end,6)),' Joules'));
hold on
set(gcf,'Color','w');
saveas(gcf,'Acoustic_Emissions','tiffn');
% saveas(gcf,'Acoustic_Emissions','eps');

%%
%b-value calculations

% AE_vector=master1(:,19);
% max_AE=max(AE_vector);
% min_AE=min(AE_vector);
% bins=100;
% stepsize=(max_AE-min_AE)/bins; %100 is the number of bins
% N_vector=zeros(bins,1);
% A_vector=(linspace(min_AE,max_AE,bins))';
% threshold=min_AE;
% 
% for i=1:bins
%     [loc_AE]=find(AE_vector(:)>=threshold);
%     N_vector(i)=log10(length(loc_AE));
%     threshold=threshold+stepsize;
% end
% 
% f87=figure(87);
% plot(A_vector,N_vector,'o-')
% xlabel('AE Magnitude'); ylabel('N(>Amax)');
% title('b-value plot');
% % set(gcf,'Position',get(0,'Screensize'))
% hold on
% 
% satisfied=0;
% while satisfied==0
%     disp('Select x-range for best fit line, lower limit first');
%     [x_range,~]=ginput(2);
% 
%     [AE_locs1]=find(A_vector(:)>=x_range(1));
%     [AE_locs2]=find(A_vector(:)<=x_range(2));
%     [AE_locs]=intersect(AE_locs1,AE_locs2);
%     bfit_N=N_vector(AE_locs);
%     bfit_A=A_vector(AE_locs);
% 
%     P=polyfit(bfit_A,bfit_N,1);
%     plotP=polyval(P,bfit_A);
%     bfp=plot(bfit_A,plotP,'r--','LineWidth',3);
%     hold on
%     legend('Magnitude-Frequency relationship','Best fit line','Location','SouthWest');
%     hold on
% 
%     b_value=-P(1);
%     disp(strcat('b value: ',num2str(b_value)));
%     satisfied=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
%     if satisfied==0
%         delete(bfp);
%     end
% end
% saveas(gcf,'b_value_unclumped','tiffn');

%%
%D-value
r_vector=(linspace(0,max_X,100))';   %change number of bins for speed
cr=zeros(length(r_vector),1);
n_count=zeros(length(r_vector),1);

for i=1:length(cr)
    for j=1:p1
        for k=1:p1
            r_len1=sqrt((master1(j,6)-master1(k,6))^2+(master1(j,7)-master1(k,7))^2);
            r_len2=sqrt((master1(j,9)-master1(k,9))^2+(master1(j,10)-master1(k,10))^2);
            r_len3=sqrt((master1(j,6)-master1(k,6))^2+(master1(j,10)-master1(k,10))^2);
            r_len4=sqrt((master1(j,9)-master1(k,9))^2+(master1(j,7)-master1(k,7))^2);
            if r_len1 <= r_vector(i) && r_len2 <= r_vector(i) && r_len3 <= r_vector(i) && r_len4 <= r_vector(i) && j~=k
                n_count(i)=n_count(i)+1;
            end
        end
    end
    cr(i)=(2*n_count(i))/(p1*(p1-1));
end

log_cr=log10(cr);
log_r=log10(r_vector);
% diff_cr=diff(log_cr);
% diff_r=diff(r_vector);

% use log_cr and r_vector for plotting
%use log_cr and log_r for regression

f13=figure(13);
semilogx(r_vector,log_cr,'LineWidth',3);
hold on
axis tight

% subplot(1,2,2)
% semilogx(diff_r,diff_cr);
% hold on

r_satisfied=0;
while r_satisfied==0
    disp('Select x-range for Constant Slope, lower limit first');
    [r_x_range,~]=ginput(2);
    
    [r_locs1]=find(r_vector(:)>=r_x_range(1));
    [r_locs2]=find(r_vector(:)<=r_x_range(2));
    [r_locs]=intersect(r_locs1,r_locs2);
    r_bfit=log_r(r_locs);
    r_bfit_plot=r_vector(r_locs);
    cr_bfit=log_cr(r_locs);
    
    r_P=polyfit(r_bfit,cr_bfit,1);
    r_plotP=polyval(r_P,r_bfit);
    %     log_r_plot=log10(r_plotP);
    r_bfp=semilogx(r_bfit_plot,r_plotP,'r--','LineWidth',5);
    %     r_bfp=semilogx(r_bfit,r_plotP,'r--');
    hold on
    legend('Spatial Damage Correlation relationship','Best fit line','Location','SouthEast');
    hold on
    %     legend('Magnitude-Frequency relationship','Best fit line');
    %     hold on
    xlabel('radius'); ylabel('log(C(r))');
    title('D-value plot');
    
    d_value=r_P(1);
    disp(strcat('D-value: ',num2str(d_value)));
    r_satisfied=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
    if r_satisfied==0
        delete(r_bfp);
    end
end
saveas(gcf,'D_value','tiffn');
saveas(gcf,'D_value','epsc');
export_fig fig13 f13_D_value_nonclumped -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%D-value tensile
r_vector=(linspace(0,max_X,100))';   %change number of bins for speed
cr=zeros(length(r_vector),1);
n_count=zeros(length(r_vector),1);

for i=1:length(cr)
    for j=1:n_tensile
        for k=1:n_tensile
            r_len1=sqrt((master1_tensile(j,6)-master1_tensile(k,6))^2+(master1_tensile(j,7)-master1_tensile(k,7))^2);
            r_len2=sqrt((master1_tensile(j,9)-master1_tensile(k,9))^2+(master1_tensile(j,10)-master1_tensile(k,10))^2);
            r_len3=sqrt((master1_tensile(j,6)-master1_tensile(k,6))^2+(master1_tensile(j,10)-master1_tensile(k,10))^2);
            r_len4=sqrt((master1_tensile(j,9)-master1_tensile(k,9))^2+(master1_tensile(j,7)-master1_tensile(k,7))^2);
            if r_len1 <= r_vector(i) && r_len2 <= r_vector(i) && r_len3 <= r_vector(i) && r_len4 <= r_vector(i) && j~=k
                n_count(i)=n_count(i)+1;
            end
        end
    end
    cr(i)=(2*n_count(i))/(n_tensile*(n_tensile-1));
end

log_cr=log10(cr);
log_r=log10(r_vector);
% diff_cr=diff(log_cr);
% diff_r=diff(r_vector);

% use log_cr and r_vector for plotting
%use log_cr and log_r for regression

f19=figure(19);
semilogx(r_vector,log_cr,'LineWidth',3);
hold on
axis tight

r_satisfied=0;
while r_satisfied==0
    disp('Select x-range for Constant Slope, lower limit first');
    [r_x_range,~]=ginput(2);
    
    [r_locs1]=find(r_vector(:)>=r_x_range(1));
    [r_locs2]=find(r_vector(:)<=r_x_range(2));
    [r_locs]=intersect(r_locs1,r_locs2);
    r_bfit=log_r(r_locs);
    r_bfit_plot=r_vector(r_locs);
    cr_bfit=log_cr(r_locs);
    
    r_P=polyfit(r_bfit,cr_bfit,1);
    r_plotP=polyval(r_P,r_bfit);
    %     log_r_plot=log10(r_plotP);
    r_bfp=semilogx(r_bfit_plot,r_plotP,'r--','LineWidth',5);
    %     r_bfp=semilogx(r_bfit,r_plotP,'r--');
    hold on
    legend('Spatial Damage Correlation relationship','Best fit line','Location','SouthEast');
    hold on
    %     legend('Magnitude-Frequency relationship','Best fit line');
    %     hold on
    xlabel('radius'); ylabel('log(C(r))');
    title('D-value plot- TENSILE ONLY');
    
    d_value_t=r_P(1);
    disp(strcat('D-value: ',num2str(d_value_t)));
    r_satisfied=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
    if r_satisfied==0
        delete(r_bfp);
    end
end
saveas(gcf,'D_value_tensile','tiffn');
saveas(gcf,'D_value_tensile','epsc');
export_fig fig19 f19_D_value_nonclumped_tensile -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%D-value shear
r_vector=(linspace(0,max_X,100))';   %change number of bins for speed
cr=zeros(length(r_vector),1);
n_count=zeros(length(r_vector),1);

for i=1:length(cr)
    for j=1:n_shear
        for k=1:n_shear
            r_len1=sqrt((master1_shear(j,6)-master1_shear(k,6))^2+(master1_shear(j,7)-master1_shear(k,7))^2);
            r_len2=sqrt((master1_shear(j,9)-master1_shear(k,9))^2+(master1_shear(j,10)-master1_shear(k,10))^2);
            r_len3=sqrt((master1_shear(j,6)-master1_shear(k,6))^2+(master1_shear(j,10)-master1_shear(k,10))^2);
            r_len4=sqrt((master1_shear(j,9)-master1_shear(k,9))^2+(master1_shear(j,7)-master1_shear(k,7))^2);
            if r_len1 <= r_vector(i) && r_len2 <= r_vector(i) && r_len3 <= r_vector(i) && r_len4 <= r_vector(i) && j~=k
                n_count(i)=n_count(i)+1;
            end
        end
    end
    cr(i)=(2*n_count(i))/(n_shear*(n_shear-1));
end

log_cr=log10(cr);
log_r=log10(r_vector);
% diff_cr=diff(log_cr);
% diff_r=diff(r_vector);

% use log_cr and r_vector for plotting
%use log_cr and log_r for regression

f20=figure(20);
semilogx(r_vector,log_cr,'LineWidth',3);
hold on
axis tight

r_satisfied=0;
while r_satisfied==0
    disp('Select x-range for Constant Slope, lower limit first');
    [r_x_range,~]=ginput(2);
    
    [r_locs1]=find(r_vector(:)>=r_x_range(1));
    [r_locs2]=find(r_vector(:)<=r_x_range(2));
    [r_locs]=intersect(r_locs1,r_locs2);
    r_bfit=log_r(r_locs);
    r_bfit_plot=r_vector(r_locs);
    cr_bfit=log_cr(r_locs);
    
    r_P=polyfit(r_bfit,cr_bfit,1);
    r_plotP=polyval(r_P,r_bfit);
    %     log_r_plot=log10(r_plotP);
    r_bfp=semilogx(r_bfit_plot,r_plotP,'r--','LineWidth',5);
    %     r_bfp=semilogx(r_bfit,r_plotP,'r--');
    hold on
    legend('Spatial Damage Correlation relationship','Best fit line','Location','SouthEast');
    hold on
    %     legend('Magnitude-Frequency relationship','Best fit line');
    %     hold on
    xlabel('radius'); ylabel('log(C(r))');
    title('D-value plot- SHEAR ONLY');
    
    d_value_s=r_P(1);
    disp(strcat('D-value: ',num2str(d_value_s)));
    r_satisfied=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
    if r_satisfied==0
        delete(r_bfp);
    end
end
saveas(gcf,'D_value_shear','tiffn');
saveas(gcf,'D_value_shear','epsc');
export_fig fig19 f19_D_value_nonclumped_shear -q101 -painters -nocrop -pdf -png -tiff -eps 

keyboard
%%
%clumping of events

disp('CLUMP EVENTS');
disp(strcat('Min Cycle: ',num2str(n1)));
disp(strcat('Max Cycle: ',num2str(n2)));
% c_limit=input('Enter cycle event clumping limit: ');

X_diff=max_X-min_X;

disp(strcat('Min X: ',num2str(min_X)));
disp(strcat('Max X: ',num2str(max_X)));
disp(strcat('X range: ',num2str(X_diff)));
disp(strcat('Max particle radius: ',num2str(max_r)));
disp(strcat('Min particle radius: ',num2str(min_r)));
% r_limit=input('Enter distance event clumping limit: ');

master1_buff=master1;
master1_nonclumped=master1;

[len,~]=size(master1);

dist_count=0;
clump_events=zeros(len,15);
nonclump_index=zeros(2*len,1);
clump_count=0;

% clump_events
% 1 = x1 of first event
% 2 = y1 of fitst event
% 3 = x2 of first event
% 4 = y2 of first event
% 5 = event cycle number for event 1
% 6 = energy for event 1
% 7 = x1 of second event
% 8 = y1 of second event
% 9 = x2 of second event
% 10 = y2 of second event
% 11 = event cycle number for event 2
% 12 = energy for event 2
% 13 = difference in step
% 14 = cumulative energy
% 15 = megaevent number
% 16=incident 1 ID
% 17 = incident 2 ID

% PAIR MICROFRACTURES
c_ev=zeros(1,6);

for i=1:length(master1_buff)
    for j=1:length(master1_buff)
        r_len1_c=sqrt((master1_buff(i,6)-master1_buff(j,6))^2+(master1_buff(i,7)-master1_buff(j,7))^2);
        r_len2_c=sqrt((master1_buff(i,9)-master1_buff(j,9))^2+(master1_buff(i,10)-master1_buff(j,10))^2);
        r_len3_c=sqrt((master1_buff(i,6)-master1_buff(j,6))^2+(master1_buff(i,10)-master1_buff(j,10))^2);
        r_len4_c=sqrt((master1_buff(i,9)-master1_buff(j,9))^2+(master1_buff(i,7)-master1_buff(j,7))^2);
        cdiff=abs(master1_buff(i,15)-master1_buff(j,15));
        if i~=j && cdiff <=c_limit
            if r_len1_c<=r_limit || r_len2_c<=r_limit || r_len3_c<=r_limit || r_len4_c<=r_limit  %check if pair falls within radius and cycle limits
                
                clump_count=clump_count+1;
                nonclump_index(clump_count)=i;
                
                c_ev(clump_count,1)=master1_buff(i,6);
                c_ev(clump_count,2)=master1_buff(i,7);
                c_ev(clump_count,3)=master1_buff(i,9);
                c_ev(clump_count,4)=master1_buff(i,10);
                c_ev(clump_count,5)=master1_buff(i,15);
                c_ev(clump_count,6)=master1_buff(i,18);
                
                clump_count=clump_count+1;
                nonclump_index(clump_count)=j;
                
                c_ev(clump_count,1)=master1_buff(j,6);
                c_ev(clump_count,2)=master1_buff(j,7);
                c_ev(clump_count,3)=master1_buff(j,9);
                c_ev(clump_count,4)=master1_buff(j,10);
                c_ev(clump_count,5)=master1_buff(j,15);
                c_ev(clump_count,6)=master1_buff(j,18);
                
                dist_count=dist_count+1;
                
                clump_events(dist_count,1)=master1_buff(i,6);
                clump_events(dist_count,2)=master1_buff(i,7);
                clump_events(dist_count,3)=master1_buff(i,9);
                clump_events(dist_count,4)=master1_buff(i,10);
                clump_events(dist_count,5)=master1_buff(i,15);
                clump_events(dist_count,6)=master1_buff(i,18);
                
                clump_events(dist_count,7)=master1_buff(j,6);
                clump_events(dist_count,8)=master1_buff(j,7);
                clump_events(dist_count,9)=master1_buff(j,9);
                clump_events(dist_count,10)=master1_buff(j,10);
                clump_events(dist_count,11)=master1_buff(j,15);
                clump_events(dist_count,12)=master1_buff(j,18);
                
                clump_events(dist_count,13)=abs(clump_events(dist_count,5)-clump_events(dist_count,11));
                clump_events(dist_count,14)=clump_events(dist_count,6)+clump_events(dist_count,12);
                
                clump_events(dist_count,15)=0;
                
                clump_events(dist_count,16)=master1_buff(i,1);
                clump_events(dist_count,17)=master1_buff(j,1);
            end
        end
    end
    %master1_buff is array of non-clumped bond breakages
end

% forclump=zeros(length(master1,10));
%1=incident id
%2=x1
%3=y1
%4=y2
%5=x2
%6=energy
%7=cycle
%8=strain
%9= clumped counter (1 if clumped, 0 if not)
%10= megaevent number

nonclump_index=unique(nonclump_index);
master1_nonclumped(nonclump_index,:)=[];

m1_nc=zeros(length(master1_nonclumped),6);
m1_nc(:,1)=master1_nonclumped(:,6);
m1_nc(:,2)=master1_nonclumped(:,7);
m1_nc(:,3)=master1_nonclumped(:,9);
m1_nc(:,4)=master1_nonclumped(:,10);
m1_nc(:,5)=master1_nonclumped(:,15);
m1_nc(:,6)=master1_nonclumped(:,18);

clump_events=unique(clump_events,'rows');

c_check1=intersect(c_ev,m1_nc,'rows');

if isempty(c_check1)==1
    disp('No overlap between clumping and non clumping events!');
end

clump_events2=clump_events;
[cl2,~]=size(clump_events2);
% clump_events 2 is a matrix of all the bonds broken within the specified
% radius and within the specified cycle limit

clump_master=zeros(cl2,6);


%clump_master
%1 = megaevent
%2 = no of individual events
%3 = megaevent energy
%4 = megaevent avg event x
%5 = megaevent avg event y
%6 = megaevent avg event cycle
%7 to 19 = same as clump_event2

total_megaevents=0;
c_index=zeros(length(clump_events),1);
cce=0;
eventchecker=master1(:,1);
eventchecker(:,2)=0;

for i=1:length(clump_events)
    if c_index(i)==0
        cce=cce+2;
        newid=[clump_events(i,16), clump_events(i,17)];
        if i==1
            newid_master=newid;
        else
            newid_master=[newid_master,newid];
        end
        
        newbuff=[];
        
        total_megaevents=total_megaevents+1;
        me_ec=2;
        megaevents_energy=clump_events2(i,14);
        xsum=mean([clump_events(i,1),clump_events(i,3)])+mean([clump_events(i,7),clump_events(i,9)]);
        ysum=mean([clump_events(i,2),clump_events(i,4)])+mean([clump_events(i,8),clump_events(i,10)]);
        csum=min(clump_events(i,5),clump_events(i,11));
        
        while isempty(newid)~=1
            newid_master=[newid_master,newid];
            newbuff2=[];
            for j=1:length(newid)
                [a1,~]=find(clump_events(:,16)==newid(j));
                [a2,~]=find(clump_events(:,17)==newid(j));
                
                if length(a1)==1
                    ise1=0;
                else
                    ise1=1;
                end
                if length(a2)==1
                    ise2=0;
                else
                    ise2=1;
                end
                
                
                if ise1~=0
                    for k=1:length(a1)
                        if isempty(intersect(clump_events(a1(k),17),newid_master))==1 && a1(k)~=i && c_index(a1(k))==0 && isempty(intersect(clump_events(a1(k),17),newbuff2))==1
                            cce=cce+1;
                            me_ec=me_ec+1;
                            megaevents_energy=megaevents_energy+clump_events(a1(k),12);
                            xsum=xsum+mean([clump_events(a1(k),7),clump_events(a1(k),9)]);
                            ysum=ysum+mean([clump_events(a1(k),8),clump_events(a1(k),10)]);
                            minmat=[csum,clump_events2(a1(k),11)];
                            csum=min(minmat);
                            newbuff2=[newbuff2,clump_events(a1(k),17)];
                            newbuff=[newbuff,clump_events(a1(k),17)];
                        end
                    end
                end
                c_index(a1)=1;
                clump_events(a1,15)=total_megaevents;
                
                if ise2~=0
                    for k=1:length(a2)
                        if isempty(intersect(clump_events(a2(k),16),newid_master))==1 && a2(k)~=i && c_index(a2(k))==0 && isempty(intersect(clump_events(a2(k),16),newbuff2))==1
                            cce=cce+1;
                            me_ec=me_ec+1;
                            megaevents_energy=megaevents_energy+clump_events(a2(k),12);
                            xsum=xsum+mean([clump_events(a2(k),7),clump_events(a2(k),9)]);
                            ysum=ysum+mean([clump_events(a2(k),8),clump_events(a2(k),10)]);
                            minmat=[csum,clump_events(a2(k),11)];
                            csum=min(minmat);
                            newbuff2=[newbuff2,clump_events(a2(k),16)];
                            newbuff=[newbuff,clump_events(a2(k),16)];
                        end
                    end
                end
                c_index(a2)=1;
                clump_events(a2,15)=total_megaevents;
                
                
            end
            
            newbuff=unique(newbuff);
            newid_buff=setdiff(newbuff,newid_master);
            newid=newid_buff;
        end
        
        clump_master(i,1)=total_megaevents;
        clump_master(i,2)=me_ec;
        clump_master(i,3)=megaevents_energy;
        clump_master(i,4)=xsum./me_ec;
        clump_master(i,5)=ysum./me_ec;
        clump_master(i,6)=round(csum);
    end
end
%%
todelete_clumped=clump_master(:,1)==0;
clump_master(todelete_clumped,:)=[];
tomod=clump_events(:,15)==0;

for i=1:length(clump_master)
    %     clump_master(i,7)=(log10(clump_master(i,3))-4.4)./1.5;       %Magnitude calculation
    clump_master(i,7)=(0.67.*log10(clump_master(i,3)))-2.9;
end
no_megaevents=length(clump_master);

% %clump_master
% %1 = megaevent
% %2 = no of individual events
% %3 = megaevent energy
% %4 = megaevent avg event x
% %5 = megaevent avg event y
% %6 = megaevent avg event cycle
% %7 = magnitude

no_singleevents=length(master1_nonclumped);

checksum=no_megaevents+no_singleevents;
disp(strcat('No of total events: ',num2str(p1)));
disp(strcat('No of events after clumping: ',num2str(checksum)));
disp(strcat('Single Events: ',num2str(no_singleevents)));
disp(strcat('Clumped Events: ',num2str(no_megaevents)));

clump_master(no_megaevents+1:no_megaevents+no_singleevents,:)=0;

for i=1:no_singleevents
    clump_master(no_megaevents+i,1)=no_megaevents+i;
    master1_nonclumped(i,20)=clump_master(no_megaevents+i,1);
    clump_master(no_megaevents+i,2)=1;
    clump_master(no_megaevents+i,3)=master1_nonclumped(i,18);
    clump_master(no_megaevents+i,4)=(master1_nonclumped(i,6)+master1_nonclumped(i,9))./2;
    clump_master(no_megaevents+i,5)=(master1_nonclumped(i,7)+master1_nonclumped(i,10))./2;
    clump_master(no_megaevents+i,6)=master1_nonclumped(i,15);
    clump_master(no_megaevents+i,7)=master1_nonclumped(i,19);
    [yloc_a3]=find(eventchecker(:,1)==master1_nonclumped(i,1),1);
    eventchecker(yloc_a3,2)=1;
end


%%
max_eq_clumped=max(clump_master(:,7));
disp(strcat('Max Eq Magnitude: ',num2str(max_eq_clumped)));

[~,sorted_inds]=sort(clump_master(:,6));
clump_master1=clump_master(sorted_inds,:);

for i=length(clump_master1)
    cyc=clump_master1(i,6);
    [cycrow,~]=find(cyc==master2(:,1));
    clump_master1(:,8)=master2(cycrow,13);
end

%clump_master, clump_master1 is sorted version by cycle

%1 = megaevent
%2 = no of individual events
%3 = megaevent energy
%4 = megaevent avg event x
%5 = megaevent avg event y
%6 = megaevent avg event cycle
%7 = eq magnitude

clump_master2=zeros(n2-n1+1,11);
for i=1:length(clump_master2)
    clump_master2(i)=i;
end

clump_master2(:,8)=master2(:,13);

% sum_AE_c=0;
% sum_energy_c=0;

for i=1:length(clump_master2)
    [yloc_c]=find(i==clump_master1(:,6));
    if isempty(yloc_c)~=1
        clump_master2(i,2)=sum(clump_master1(yloc_c,2));
        clump_master2(i,5)=sum(clump_master1(yloc_c,3));
    end
end

for i=1:length(clump_master2)
    if i==1
        clump_master2(i,3)=clump_master2(i,2);
        clump_master2(i,6)=clump_master2(i,5);
    else
        clump_master2(i,3)=clump_master2(i,2)+clump_master2(i-1,3);
        clump_master2(i,6)=clump_master2(i,5)+clump_master2(i-1,6);
    end
end
tot_events=clump_master2(end,3);
tot_energy=clump_master2(end,6);

max_AE_c=max(clump_master2(:,2));
max_energy_c=max(clump_master2(:,5));

for i=1:length(clump_master2)
    clump_master2(i,4)=clump_master2(i,3)./tot_events;
    clump_master2(i,7)=clump_master2(i,6)./tot_energy;
end

for i=n1:n2
    clump_master2(i,9)=clump_master2(i,2)./max_AE_c;
    clump_master2(i,10)=clump_master2(i,5)./max_energy_c;
    %     clump_master2(i,11)=(log10(clump_master2(i,5))-4.4)/1.5;
    clump_master2(i,11)=(0.67.*log10(clump_master2(i,5)))-2.9;
    %     clump_master2(i,11)=log10(clump_master2(i,5));
end

[p1_c,~]=size(clump_master1);
%%
% clump_master2
% 1 = cycle
% 2 = no of events
% 3 = cumulative events
% 4 = fractional events
% 5 = cycle events energy
% 6 = cumulative energy
% 7 = fractional energy
% 8 = strain
% 9 = AE count/Max AE count
% 10 = cycle energy/ max energy
% 11 = eq magnitude

figure
subplot(2,1,1)
f1=bar(clump_master2(:,8),clump_master2(:,9));
set(f1,'FaceColor','b','EdgeColor','k')
axis tight
hold on
% plot(clump_master2(:,8),clump_master2(:,4),'r','LineWidth',2);
[af1,h1,h2]=plotyy(clump_master2(:,8),clump_master2(:,4),clump_master2(:,8),clump_master2(:,2));
axis tight
set(af1(2),'XTickLabel',[],'XTick',[]);
delete(h2);
set(h1,'Color','r','LineWidth',2);
set(af1,{'ycolor'},{'k';'k'});

xlabel('Axial Strain');ylabel('AE');
legend('AE','Cumulative AE','Location','NorthEast','Orientation','horizontal');
title(strcat('Total AE event count: ',num2str(p1_c)));
hold on

subplot(2,1,2)
f2=bar(clump_master2(:,8),clump_master2(:,10));
set(f2,'FaceColor','g','EdgeColor','k')
axis tight
hold on
% plot(clump_master2(:,8),clump_master2(:,7),'r','LineWidth',2);
[af2,h1,h2]=plotyy(clump_master2(:,8),clump_master2(:,7),clump_master2(:,8),clump_master2(:,5));
axis tight
set(af2(2),'XTickLabel',[],'XTick',[]);
delete(h2);
set(h1,'Color','r','LineWidth',2);
set(af2,{'ycolor'},{'k';'k'});
xlabel('Axial Strain');ylabel('Energy Ratio');
legend('Fractional E ratio','Cumulative Energy','Location','NorthEast','Orientation','horizontal');
title(strcat('Total Energy Released: ',num2str(clump_master2(end,6)),' Joules'));
hold on
set(gcf,'Color','w');
saveas(gcf,'Acoustic_Emissions-CLUMPED','tiffn');
% saveas(gcf,'Acoustic_Emissions-CLUMPED','eps');
%%
% clumped_master1_buff=zeros(p1,7);
% clumped_master1_buff
% C1 = step
% C2 = x1 avg
% C3 = y1 avg
% C4 = x2 avg
% C5 = y2 avg
% C6 = energy
% C7 = Magnitude
%%
%b-value calculations - clumped

AE_vector_c=clump_master1(:,7);
max_AE_c=max(AE_vector_c);
min_AE_c=min(AE_vector_c);
bins_c=100;
stepsize_c=(max_AE_c-min_AE_c)/bins_c; %100 is the number of bins
N_vector_c=zeros(bins_c,1);
A_vector_c=(linspace(min_AE_c,max_AE_c,bins_c))';
threshold_c=min_AE_c;

for i=1:bins_c
    [loc_AE_c]=find(AE_vector_c(:)>=threshold_c);
    N_vector_c(i)=log10(length(loc_AE_c));
    threshold_c=threshold_c+stepsize_c;
end

fig9=figure(9);
plot(A_vector_c,N_vector_c,'o-')
xlabel('AE Magnitude'); ylabel('N(>Amax)');
title('b-value plot - CLUMPED');
hold on

satisfied=0;
while satisfied==0
    disp('Select x-range for best fit line, lower limit first');
    [x_range_c,~]=ginput(2);
    
    [AE_locs1_c]=find(A_vector_c(:)>=x_range_c(1));
    [AE_locs2_c]=find(A_vector_c(:)<=x_range_c(2));
    [AE_locs_c]=intersect(AE_locs1_c,AE_locs2_c);
    bfit_N_c=N_vector_c(AE_locs_c);
    bfit_A_c=A_vector_c(AE_locs_c);
    P_c=polyfit(bfit_A_c,bfit_N_c,1);
    plotP_c=polyval(P_c,bfit_A_c);
    bfp_c=plot(bfit_A_c,plotP_c,'r--','LineWidth',3);
    hold on
    legend('Magnitude-Frequency relationship','Best fit line','Location','SouthWest');
    hold on
    
    b_value_c=-P_c(1);
    disp(strcat('b value: ',num2str(b_value_c)));
    satisfied=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
    if satisfied==0
        delete(bfp_c);
    end
end
set(gcf,'Color','w');
saveas(gcf,'b_value_clumped','tiffn');
saveas(gcf,'b_value_clumped','epsc');
export_fig fig9 f9_bvalue -pdf -q101 -painters -nocrop -pdf -png -tiff -eps 

%%

% D-value

r_vector_c=(linspace(0,max_X,100))';
cr_c=zeros(length(r_vector_c),1);
n_count_c=zeros(length(r_vector_c),1);

for i=1:length(cr_c)
    for j=1:p1_c
        for k=1:p1_c
            r_len1_c=sqrt((clump_master1(j,4)-clump_master1(k,4))^2+(clump_master1(j,5)-clump_master1(k,5))^2);
            %             r_len2_c=sqrt((clumped_master1(j,4)-clumped_master1(k,4))^2+(clumped_master1(j,5)-clumped_master1(k,5))^2);
            %             r_len3_c=sqrt((clumped_master1(j,2)-clumped_master1(k,2))^2+(clumped_master1(j,5)-clumped_master1(k,5))^2);
            %             r_len4_c=sqrt((clumped_master1(j,4)-clumped_master1(k,4))^2+(clumped_master1(j,2)-clumped_master1(k,2))^2);
            if r_len1_c <= r_vector_c(i) && j~=k
                n_count_c(i)=n_count_c(i)+1;
            end
        end
    end
    cr_c(i)=(2*n_count_c(i))/(p1_c*(p1_c-1));
end

log_cr_c=log10(cr_c);
log_r_c=log10(r_vector_c);
% diff_cr=diff(log_cr);
% diff_r=diff(r_vector);

% use log_cr and r_vector for plotting
% use log_cr and log_r for regression

log_cr_c=log10(cr_c);
log_r_c=log10(r_vector_c);
% diff_cr=diff(log_cr);
% diff_r=diff(r_vector);

% use log_cr and r_vector for plotting
%use log_cr and log_r for regression

fig10=figure(10);
semilogx(r_vector_c,log_cr_c,'LineWidth',3);
hold on
axis tight

% subplot(1,2,2)
% semilogx(diff_r,diff_cr);
% hold on

r_satisfied_c=0;
while r_satisfied_c==0
    disp('Select x-range for Constant Slope, lower limit first');
    [r_x_range_c,~]=ginput(2);
    
    [r_locs1_c]=find(r_vector_c(:)>=r_x_range_c(1));
    [r_locs2_c]=find(r_vector_c(:)<=r_x_range_c(2));
    [r_locs_c]=intersect(r_locs1_c,r_locs2_c);
    r_bfit_c=log_r_c(r_locs_c);
    r_bfit_plot_c=r_vector_c(r_locs_c);
    cr_bfit_c=log_cr_c(r_locs_c);
    
    r_P_c=polyfit(r_bfit_c,cr_bfit_c,1);
    r_plotP_c=polyval(r_P_c,r_bfit_c);
    r_bfp_c=semilogx(r_bfit_plot_c,r_plotP_c,'r--','LineWidth',5);
    hold on
    legend('Spatial Damage Correlation relationship','Best fit line','Location','SouthEast');
    hold on
    xlabel('radius (m)'); ylabel('log(C(r))');
    title('D-value plot - CLUMPED');
    
    d_value_c=r_P_c(1);
    disp(strcat('D-value Coupled: ',num2str(d_value_c)));
    r_satisfied_c=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
    if r_satisfied_c==0
        delete(r_bfp_c);
    end
end
set(gcf,'Color','w');
saveas(gcf,'D_value_clumped','tiffn');
saveas(gcf,'D_value_clumped','epsc');
export_fig fig10 f10_dvalue -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%histogram of fracture lengths

%1 = megaevent
%2 = no of individual events
%3 = megaevent energy
%4 = megaevent avg event x
%5 = megaevent avg event y
%6 = megaevent avg event cycle
%7 = eq magnitude

% [count,bins]=hist(clump_master1(:,2),20);
% figure
% ff=bar(bins,count);
% set(ff,'FaceColor','y','EdgeColor','k')
% ylabel('Number of Events');xlabel('Clumping extent (events)');
% hold on
% saveas(gcf,'events_histogram','tiffn');

%%
%clump_vid
%1 = x1
%2 = y1
%3 = x2
%4 = y2
%5 = megaevent number
%6 = megaevent cycle
%7 = megaevent strain
%8 = event counter, 1= megaevent and 0 = individual event

v_m1=clump_events(:,[1,2,3,4,15]);
v_m2=clump_events(:,[7,8,9,10,15]);

clump_vid=[v_m1;v_m2];
clump_vid=unique(clump_vid,'rows');

for i=1:length(clump_vid)
    [yloc9,~]=find(clump_vid(i,5)==clump_master(:,1),1);
    clump_vid(i,6)=clump_master(yloc9,6);
end

lint=length(clump_vid);

for i=1:no_singleevents
    clump_master(no_megaevents+i,1)=no_megaevents+i;
    clump_master(no_megaevents+i,2)=1;
    clump_master(no_megaevents+i,3)=master1_nonclumped(i,18);
    clump_master(no_megaevents+i,4)=(master1_nonclumped(i,6)+master1_nonclumped(i,9))./2;
    clump_master(no_megaevents+i,5)=(master1_nonclumped(i,7)+master1_nonclumped(i,10))./2;
    clump_master(no_megaevents+i,6)=master1_nonclumped(i,15);
    clump_master(no_megaevents+i,7)=master1_nonclumped(i,19);
    
    clump_vid(lint+i,1)=master1_nonclumped(i,6);
    clump_vid(lint+i,2)=master1_nonclumped(i,7);
    clump_vid(lint+i,3)=master1_nonclumped(i,9);
    clump_vid(lint+i,4)=master1_nonclumped(i,10);
    clump_vid(lint+i,5)=master1_nonclumped(i,20);
    clump_vid(lint+i,6)=master1_nonclumped(i,15);
end

for i=1:length(clump_vid)
    [yloc10,~]=find(clump_vid(i,6)==master2(:,1),1);
    clump_vid(i,7)=master2(yloc10,13);
end

clump_vid(1:lint,8)=1;
clump_vid(lint+1:end,8)=0;

clump_vid=sortrows(clump_vid,6);

% plot_buff23=input('Ready for Spatial Distribution plot ?: ');
% fig11=figure(11);
% % aviobj23=avifile('Clumped_Damage.avi','compression','None','fps',8);
% 
% cm2=jet(n2-n1+1);
% for i=n1:n2
%     [pos_x23,~]=find(i==clump_vid(:,6));
%     axis([min_X max_X min_Y max_Y]);
%     title('Cumulative Damage','fontweight','bold');
%     hold on
%     if pos_x23~=0
%         pos_pm23=clump_vid(pos_x23,:);
%         [lp,~]=size(pos_pm23);
%         for j=1:lp
%             if pos_pm23(j,8)==1
%                 plot(((pos_pm23(j,1)+pos_pm23(j,3))./2),((pos_pm23(j,2)+pos_pm23(j,4))./2),'o','color',cm2(i,:),'LineWidth',1);
%                 %                 plot(pos_pm23(j,3),pos_pm23(j,4),'o','color',cm2(i,:),'LineWidth',3);
%                 hold on
%             else
%                 plot(((pos_pm23(j,1)+pos_pm23(j,3))./2),((pos_pm23(j,2)+pos_pm23(j,4))./2),'x','color',cm2(i,:),'LineWidth',2);
%                 %                 plot(pos_pm23(j,3),pos_pm23(j,4),'x','color',cm2(i,:),'LineWidth',2);
%                 hold on
%             end
%         end
%     end
% end
% 
% caxis([0 axial_strain]);
% box on
% colorbar('SouthOutside')
% xlabel('Sample length (m)','fontweight','bold');ylabel('Sample width (m)','fontweight','bold');
% saveas(gcf,'Clumped_distribution','tiffn');
% saveas(gcf,'Clumped_distribution','eps');

%Re-animation with incremental plot
plot_buff24=input('Ready for Spatial Distribution plot ?: ');
fig=figure;
aviobj24=avifile('Clumped_Damage_with_incremental.avi','compression','None','fps',8);

cm4=jet(n2-n1+1);
for i=n1:n2
    [pos_x24,~]=find(i==clump_vid(:,6));
    sp1=subplot(2,1,1);
    axis([min_x3 max_x3 min_y3 max_y3]); pbaspect([a_ratio3 1 1]); box on;
    hold on
    if pos_x24~=0
        pos_pm24=clump_vid(pos_x24,:);
        [lp,~]=size(pos_pm24);
        for j=1:lp
            if pos_pm24(j,8)==1
                axis([min_x3 max_x3 min_y3 max_y3]); 
                hold on
                plot(((pos_pm24(j,1)+pos_pm24(j,3))./2),((pos_pm24(j,2)+pos_pm24(j,4))./2),'.','color',cm4(i,:));
                %                 plot(pos_pm24(j,3),pos_pm24(j,4),                ,'color',cm4(i,:),'LineWidth',3);
                hold on
            else
                axis([min_x3 max_x3 min_y3 max_y3]); 
                hold on
                plot(((pos_pm24(j,1)+pos_pm24(j,3))./2),((pos_pm24(j,2)+pos_pm24(j,4))./2),'.','color',cm4(i,:));
                %                 plot(pos_pm24(j,3),pos_pm24(j,4),'x','color',cm4(i,:),'LineWidth',2);
                hold on
            end
        end
    end
    sp2=subplot(2,1,2); pbaspect([a_ratio3 1 1]); box on;
    title(strcat('Damage || Cycle: ',num2str(i),' , Strain: ',num2str(master2(i,13))));
    axis([min_x3 max_x3 min_y3 max_y3]); pbaspect([a_ratio3 1 1]); box on;
    hold on
    if pos_x24~=0
        pos_pm24=clump_vid(pos_x24,:);
        [lp,~]=size(pos_pm24);
        for j=1:lp
            if pos_pm24(j,8)==1
                plot(((pos_pm24(j,1)+pos_pm24(j,3))./2),((pos_pm24(j,2)+pos_pm24(j,4))./2),'.','color',cm4(i,:));
                %                 plot(pos_pm24(j,3),pos_pm24(j,4),'o','color',cm4(i,:),'LineWidth',3);
                hold on
            else
                plot(((pos_pm24(j,1)+pos_pm24(j,3))./2),((pos_pm24(j,2)+pos_pm24(j,4))./2),'.','color',cm4(i,:));
                %                 plot(pos_pm24(j,3),pos_pm24(j,4),'x','color',cm4(i,:),'LineWidth',2);
                hold on
            end
        end
    end
    F24=getframe(fig);
    aviobj24=addframe(aviobj24,F24);
    if i~=n2
        clf(sp2);
    end
end
caxis([0 axial_strain]);
sp2_pos=get(sp2,'position');
hcb=colorbar('NorthOutside');
set(sp2,'position',sp2_pos);
% cth=get(hcb,'Title');
% ts='Strain';
% set(cth,'String',ts);

F24=getframe(fig);
aviobj24=addframe(aviobj24,F24);
set(gcf,'Color','w');
saveas(gcf,'Clumped_distribution_with_incremental','tiffn');
aviobj24=close(aviobj24);
