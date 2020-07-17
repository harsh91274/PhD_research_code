%%
%for p3
%add tensile-shear?
%uncomment event rate and energy statistics


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
% % 11 = z2
% 12 = Normal Force
% 13 = Shear Force
% 14 = Shear Force
% 15 = cycle
% 16 = peak force tensile or shear
% 17 = area of failure contact
% 18 = energy release due to failure
% 19 = magnitude
% 20 = energy/unit area 

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
% 13 = axial strain
% 14 = fractional strain
% 15 = AE count/Max AE count
% 16 = cycle energy/ max energy
% 17 = eq magnitude
% 18 = energy/unit area

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
%3 = x-left
%4 = x-right
%5 = S3
%6 = y-top
%7 = y_bot
%8 = POROSITY
%9 = area of sample

%user_input
%1 = first file
%2 = second file
%3 = YM of bonds
%4 = SM of bonds
%5 = TS of bonds
%6 = C of bonds
%7 = Energy calculation choice 
%8 = Shear wave velocity
%9 = cycle clumping limit
%10 = radial clumping limit
%11 = input strain at peak S1

%output
%1 = confining pressure
%2 = stage of fracturing
%3 = file starting number
%4 = file ending number
%5 = input energy from S1 (Joules)
%6 = input energy from S3 (Joules)
%7 = total input enegy (Joules)
%8 = total input energy (J/m2)
%9 = Number of shear MF
%10 = Number of tensile shear MF
%11 = Number of tensile MF
%12 = Total MF
%13 = Fraction of MF in shear
%14 = Energy in Shear (Joules)
%15 = Energy in tensile-shear (Joules)
%16 = Energy in tensile (Joules)
%18 = Total MF energy (Joules)
%19 = Shear MF energy fraction
%20 = fracture energy/Input energy fraction
%21= fracture energy (J/m2)
%22 = Event rate std deviation ---------deleted from p2
%23 = event rate variance
%24 = event rate kurtosis
%25 = event rate foreshocks
%26 = event rate aftershocks
%27 = Event energy std deviation
%28 = event energy variance
%29 = event energy kurtosis
%30 = event energy foreshocks
%31 = event energy aftershocks ------------ till here
%32 = max S1 (MPa)
%32 = Max axial strain
%33 = strain at max shear stress 
%34 = Young's Modulus (MPa)
%35 = Youngs Modulus input for v0 calcs(MPa)
%36 = poissons ratio input for v0 calcs 
%37 = strain_energy
%38 = strain energy/ input energy
%39 = d-value clumped
%40 = average seismic moment
%41 = maximum seismic moment
%42 = max-min seismic moment
%43 = moment standard deviation
%44 = moment variance
%45 = moment kurtosis
%46 = b-value
%47 = D-value clumped

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
user_input=[]; inp_index=0;

diary_n1=cd;
diary_n2=diary_n1(end-3:end);
diary diary_n2;

output0=[];
output=[];

cp=input('Enter Confining Pressure: ');
stage=input('Enter Stage: ');

n1=input('Enter First File Number: ');
inp_index=inp_index+1; user_input(inp_index)=n1;
n2=input('Enter Last File Number: ');
inp_index=inp_index+1; user_input(inp_index)=n2;

output=[cp, stage, n1, n2];
%%
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
master1(:,q1+5)=0;
cf = input('Enter Youngs Modulus (Pa): ');
inp_index=inp_index+1; user_input(inp_index)=cf;

cf2=input('Enter Shear Modulus (Pa): ');
inp_index=inp_index+1; user_input(inp_index)=cf2;

tensile_strength=input('Enter Tensile Strength of bonds (Pa): ');
inp_index=inp_index+1; user_input(inp_index)=tensile_strength;

shear_strength=input('Enter Shear Strength of bonds (Pa): ');
inp_index=inp_index+1; user_input(inp_index)=shear_strength;

e_option=input('Calculate energy using TS and SS? (1)=Yes: ');
inp_index=inp_index+1; user_input(inp_index)=e_option;

vs=input('Enter Shear Wave Velocity of Material (m/s): ');
inp_index=inp_index+1; user_input(inp_index)=vs;

% pvel=input('Enter Wall Velocity (10^-4 m/step): ');
pvel=0.8;
% ts=input('Enter Time Step: ');
ts=2*10^-8;
% cycles=input('Enter number of cycles: ');
cycles=2500;

c_limit=input('Enter cycle event clumping limit (2): ');
inp_index=inp_index+1; user_input(inp_index)=c_limit;

r_limit=input('Enter distance event clumping limit (0.0012): ');
inp_index=inp_index+1; user_input(inp_index)=r_limit;

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
            master1(i,20)=0.5.*(1/cf).*((master1(i,16)./master1(i,17)).^2).*zave; %energy is in joules/m2
        else
            master1(i,16)=master1(i,14);
            master1(i,18)=0.5.*(1/cf2).*((master1(i,16)./master1(i,17)).^2).*master1(i,17).*zave; %energy is in joules
            master1(i,20)=0.5.*(1/cf2).*((master1(i,16)./master1(i,17)).^2).*zave; %energy is in joules/m2
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
master2=zeros(n2-n1,18);

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

master2(:,8)=newoutput(n1+1:n2+1,5)./10^6;
master2(:,9)=newoutput(n1+1:n2+1,2)./10^6;

step_disp=ts*pvel*cycles*2;
max_X=newoutput(1,4);
min_X=newoutput(1,3);
s_length=max_X-min_X;

for i=1:n2-n1+1
    master2(i,10)=(master2(i,8)+ master2(i,9))/2;
    master2(i,11)=(master2(i,9)- master2(i,8))/2;
    master2(i,12)=master2(i,1).*step_disp;
    master2(i,13)=master2(i,12)./s_length;
end

master2(:,14)=master2(:,12)./master2(end,12);

max_AE=max(master2(:,2));
max_energy=max(master2(:,5));

for i=1:n2-n1+1
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
% cycles_at_s1max=master2(loc1,1);
%%
%energy density calculation
sum_energyden=0;
e_count2=0;
for i=n1:n2
    [yloc]=find(i==cycle_index);
    e_count2=sum(master1(yloc,20));
    
    master2(i-n1+1,18)=e_count2;        %J/m2
end
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

[ll,qq]=size(master2);
bond_master2=zeros(ll,9);
bond_master2(:,1)=master2(:,1);
bond_master2(:,2)=master2(:,13);

[lll,~]=size(bond_master2);

for i=1:lll
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
%                 no_tensh=no_tensh+1;
%                 e_tensh=e_tensh+master1(bond_locs(j),18);
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
%%
%input energy calculations
newoutput_i=newoutput(1,:);
newoutput=newoutput(2:end,:);

%S3 energy
sample_area=newoutput(n1:n2,9);
sample_vol=sample_area.*zave;
sample_length=newoutput(n1:n2,4)-newoutput(n1:n2,3);
area_s3=sample_length.*zave;

y_top_initial=newoutput_i(1,6);
y_bot_initial=newoutput_i(1,7);
r_initial=y_top_initial-y_bot_initial;
y_top=newoutput(n1:n2,6);
y_bot=newoutput(n1:n2,7);
r_change=y_top-y_bot;

area_s1=zave*r_initial;

r_strain=(r_change-r_initial)./r_initial;

%S3 energy
pe=master2(:,8);
% pe2=master2(:,8).*r_strain;
pe3=trapz(r_strain,pe);

%S1 energy
sed=master2(:,9); 
% sed2=master2(:,9).*master2(:,13);
se=trapz(master2(:,13),sed);
% 
% ie=sed2+pe2;
% ie2=ie'*sample_vol;

pe4=master2(:,8).*10^6.*area_s3;    %energy from S3
pe5=trapz(r_strain,pe4);


se4=master2(:,9).*10^6.*area_s1;    %energy from S1
se5=trapz(master2(:,13),se4);

ie0=se+pe3;
ie=(pe5+se5);

disp('INPUT ENERGY');
disp(strcat('Input Strain Energy (from S1): ',num2str(se),' MJ/m2'));
disp(strcat('Input Strain Energy (from S1): ',num2str(se5),' Joules'));
disp(strcat('Input Potential Energy (from S3): ',num2str(pe3),' MJ/m2'));
disp(strcat('Input Potential Energy (from S3): ',num2str(pe5),' Joules'));
disp(strcat('Total Input Energy: ',num2str(ie0),' MJ/m2'));
disp(strcat('Total Input Energy: ',num2str(ie),' Joules'));
% disp(strcat('Total Input Energy: ',num2str(ie2),' Joules'));
% output=[output, se pe3 se+pe3]; 
output=[output, se5 pe5 ie ie0];
output0=[se5 pe5 ie ie0];
%%
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
n_sum=n_shear+n_ts+n_tensile; n_sf=(n_shear+n_ts)/n_sum;

disp(strcat('Energy from Shear MF: ',num2str(e_shear),' J'));
disp(strcat('Energy from Tensile-Shear MF: ',num2str(e_ts),' J'));
disp(strcat('Energy from Tensile MF: ',num2str(e_tensile),' J'));
e_sum=e_shear+e_ts+e_tensile; e_sf=(e_shear+e_ts)/e_sum;
e_frac0=e_sum/ie;
disp(strcat('Total Microcracking Energy: ',num2str(e_sum),' J'));
disp(strcat('Microcracking Energy (J)/Input Energy (J): ',num2str(e_frac0)));

fracture_energy1=master2(:,5)./(sample_vol);
G=sum(fracture_energy1); efrac1=G/(ie0*10^6);
% G2=sum(master2(:,18)); efrac2=G2/(ie0*10^6);

% disp(strcat('Total Input Energy: ',num2str(ie0),' MJ/m2'));
% disp(strcat('Total Input Energy: ',num2str(ie),' Joules')ie0);

disp(strcat('Fracture Energy from area calculations (G1): ',num2str(G),' J/m2'));
disp(strcat('Input/Output energy fraction (G1/IE): ',num2str(efrac1)));
% disp(strcat('Fracture Energy from microcrack calculations (G2): ',num2str(G2),' J/m2'));
% disp(strcat('Input/Output energy fraction (G2/IE): ',num2str(efrac2)));

output=[output, n_shear, n_ts, n_tensile, n_sum, n_sf, e_shear, e_ts, e_tensile, e_sum, e_sf, e_frac0, G];
output0=[output0, n_shear, n_ts, n_tensile, n_sum, n_sf, e_shear, e_ts, e_tensile, e_sum, e_sf, e_frac0, G];
%%
%Figure 1
f_color=linspecer(3);
f1_color=f_color;
f1_color(2,:)=f_color(3,:);
f1_color(3,:)=f_color(2,:);

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
    bb1=bar(bond_master2(:,2),[bond_master2(:,4) bond_master2(:,6) bond_master2(:,8)],'stack');colormap(f1_color);
    hold on
    set(bb1,'EdgeColor','k','linewidth',0.25);
    hold on
    xlim([0 max(bond_master2(:,2))]);
    set(af1(2),'ytick',[],'ycolor','b');

    set(af1(2),'ytickmode','auto');
    hold on
    title(strcat('Total Microfractures: ',num2str(p1)),'fontsize',12);
    xlabel('Strain','fontsize',12);
    ylabel('Number of Microfractures','fontsize',12,'color','k');
    ylabel(af1(2),'Axial Stress(MPa)','fontsize',12,'color','b');
    hold on
    legend('Shear', 'Tensile-Shear', 'Tensile','Location','NorthEast');
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
    
    bb1=bar(bond_master2(:,2),[bond_master2(:,5) bond_master2(:,7) bond_master2(:,9)],'stack'); colormap(f1_color);
    hold on
    set(bb1,'EdgeColor','k','linewidth',0.25);
    hold on
    xlim([0 max(bond_master2(:,2))]);
 
    set(af1(2),'ytick',[],'ycolor','b');
%     ylim(af1(2),[0 1]);
    set(af1(2),'ytickmode','auto');
    hold on
    title(strcat('Total Energy Released: ',num2str(master2(end,6)),' Joules'),'fontsize',12);
    set(gca,'yticklabel',num2str(get(gca,'ytick')','%.3f'));
%     set(gca,'xtick',[],'xcolor','b');
    xlabel('Strain','fontsize',12);
    ylabel('Energy(Joules)','fontsize',12,'color','k');
    ylabel(af1(2),'Fractional Energy','fontsize',12,'color','b');
    hold on
    hold on
    legend('Shear','Tensile-Shear', 'Tensile','Location','East');
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
        export_fig fig1 f1_bond_dist_autoscale_assigned -q101 -painters -nocrop -pdf -jpeg -eps 
    else
        saveas(gcf,'bond_type_distribution_autoscale','tiffn');
        saveas(gcf,'bond_type_distribution_autoscale','epsc');
        saveas(gcf,'bond_type_distribution_autoscale','fig');
        export_fig fig1 f1_bond_dist_autoscale -q101 -painters -nocrop -pdf -jpeg -eps
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
    bb1=bar(bond_master2(:,2),[bond_master2(:,4) bond_master2(:,6) bond_master2(:,8)],'stack');colormap(f1_color);
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
    
    bb1=bar(bond_master2(:,2),[bond_master2(:,5) bond_master2(:,7) bond_master2(:,9)],'stack'); colormap(f1_color);
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
% event rate
er_std=std(master2(:,2));
disp(strcat('Event Rate Standard Deviation: ',num2str(er_std)));
er_var=var(master2(:,2));
disp(strcat('Event Rate Variance: ',num2str(er_var)));
er_kur=kurtosis(master2(:,2));
disp(strcat('Event Rate Kurtosis: ',num2str(er_kur)));

output=[output, er_std, er_var, er_kur];
    
fig23=figure(23);
plot(master2(:,13),master2(:,2),'ro','LineWidth',2);
hold on
xlabel('Axial Strain','fontweight','bold');ylabel('Number of microfractures','fontweight','bold');
axis tight

satisfied_er=0;
while satisfied_er==0
    disp('Select x-range for Event Rate line, lower limit first');
    [x_range_er,~]=ginput(2);
    
    [locs1_er]=find(master2(:,13)>=x_range_er(1));
    [locs2_er]=find(master2(:,13)<=x_range_er(2));
    [locs_er]=intersect(locs1_er, locs2_er);
    
    xfit_er=master2(locs_er,13);
    yfit_er=master2(locs_er,2);
    er_fit=polyfit(xfit_er,yfit_er,1);
    er_fit_plot=polyval(er_fit,xfit_er);
    er_plot=plot(xfit_er,er_fit_plot,'k--','linewidth',3);
    hold on
    legend('Number of Microfractures','best fit line', 'location','best');
    er_val=er_fit(1);
    disp(strcat('Event Rate: ',num2str(er_val),' events/ 1% strain'));
    
    satisfied_er=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
    if satisfied_er==0
        delete(er_plot);
    end
end
set(gcf,'Color','w');
% set(gcf,'Position',get(0,'Screensize'))
saveas(gcf,'Event_rate','tiffn');
saveas(gcf,'Event_rate','epsc');
export_fig fig23 f23_microfracture_rate -q101 -painters -nocrop -pdf -png -tiff -eps 

output=[output, er_val, 0];
%%
% event energy statistics
ee_std=std(master2(:,5));
disp(strcat('Event Eenergy Standard Deviation: ',num2str(ee_std)));
ee_var=var(master2(:,5));
disp(strcat('Event Energy Variance: ',num2str(ee_var)));
ee_kur=kurtosis(master2(:,5));
disp(strcat('Event Energy Kurtosis: ',num2str(ee_kur)));

output=[output, ee_std, ee_var, ee_kur];

fig24=figure(24);
plot(master2(:,13),master2(:,5),'ro','LineWidth',2);
hold on
xlabel('Axial Strain','fontweight','bold');ylabel('Microfracture energy','fontweight','bold');
axis tight

satisfied_ee=0;
while satisfied_ee==0
    disp('Select x-range for Event Energy Rate line, lower limit first');
    [x_range_ee,~]=ginput(2);
    
    [locs1_ee]=find(master2(:,13)>=x_range_ee(1));  %strain
    [locs2_ee]=find(master2(:,13)<=x_range_ee(2));
    [locs_ee]=intersect(locs1_ee, locs2_ee);
    
    xfit_ee=master2(locs_ee,13);
    yfit_ee=master2(locs_ee,5);
    ee_fit=polyfit(xfit_ee,yfit_ee,1);
    ee_fit_plot=polyval(ee_fit,xfit_ee);
    ee_plot=plot(xfit_ee,ee_fit_plot,'k--','linewidth',3);
    hold on
    legend('Microfracture Energy','best fit line', 'location','best');
    ee_val=ee_fit(1);
    disp(strcat('Event Energy Rate: ',num2str(ee_val),' Joules/ 1% strain'));
    
    satisfied_ee=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
    if satisfied_ee==0
        delete(ee_plot);
    end
end

set(gcf,'Color','w');
% set(gcf,'Position',get(0,'Screensize'))
saveas(gcf,'Event_energy_rate','tiffn');
saveas(gcf,'Event_energy_rate','epsc');
export_fig fig24 f24_energy_rate -q101 -painters -nocrop -pdf -png -tiff -eps 

output=[output, ee_val, 0];
%%
mfwrite=master2(:,[13,2,5]);
fname=strcat('mf_statistics_',num2str(n1),'_',num2str(n2),'.dat');
dlmwrite(fname,mfwrite,'delimiter','\t');
%%
%correlate cycles to strain for positional correlation

for i=1:p1
    [pos_row,~]=find(positions(i,1)==master2(:,1),1,'first');
    positions(i,2)=master2(pos_row,13);     %record strain
end
positions(:,7)=master1(:,18);

disp('STRESS AND STRAIN RESULTS');
disp(strcat('Number of bonds broken: ', num2str(p1)));
disp(strcat('Max S1: ',num2str(max_S1)));  output=[output, max_S1];
disp(strcat('Max Axial Strain: ',num2str(axial_strain))); output=[output, axial_strain];
disp(strcat('Max Shear: ',num2str(smax)));
disp(strcat('Mean Confining Stress (S3): ',num2str(mean_s3)));
disp(strcat('Mean Stress at Max Shear: ',num2str(mean_at_smax)));
disp(strcat('Cycles at Max Shear: ',num2str(cycles_at_smax)));
disp(strcat('Strain at Max Shear: ',num2str(strain_at_smax))); output=[output, strain_at_smax];
disp(strcat('S3 at Max Shear: ',num2str(s3_at_smax)));
disp(strcat('Max Earthquake Magnitude (Events not clumped): ',num2str(max_eq)));
disp(strcat('S1 at Inelastic deformation onset: ',num2str(s1_at_in)));
disp(strcat('Cycles at Inelastic deformation onset: ',num2str(cycles_at_in)));
disp(strcat('Strain at Inelastic deformation onset: ',num2str(strain_at_in)));
disp(strcat('Mean Stress at Inelastic deformation onset: ',num2str(mean_at_in)));

%%
fig2=figure(2);
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

output=[output, YM_val];
%%
%strain energy calculations
v=r_strain./master2(:,13);
fig25=figure(25);
plot(master2(:,13),v,'linewidth',3); hold on;
xlabel('Axial Strain');ylabel('Poissons Ratio');


v_satisfied=0;
while v_satisfied==0
    disp('Select x-range for Poissons Ratio line, lower limit first');
    [x_range_v,~]=ginput(2);
    
    [locs1_v]=find(master2(:,13)>=x_range_v(1));
    [locs2_v]=find(master2(:,13)<=x_range_v(2));
    [locs_v]=intersect(locs1_v, locs2_v);
    
    
    v_range=v(locs_v);
%     
%     xfit_v=master2(locs_v,13);
%     yfit_v=master2(locs_v,9);
%     v_fit=polyfit(xfit_v,yfit_v,1);
%     v_fit_plot=polyval(v_fit,xfit_v);
%     v_plot=plot(xfit_v,v_fit_plot,'k--','linewidth',3);
%     hold on
%     legend('Poissons Ratio','best fit','location','best');
%     v_val=v_fit(1);
%     
    v_val=mean(v_range);

    disp(strcat('Poissons Ratio: ',num2str(v_val)));
    
    v_satisfied=input('Satisfied with fit ? (Yes = 1, No = 0) : ');

end

set(gcf,'Color','w');
saveas(gcf,'poissons_ratio','tiffn');
saveas(gcf,'poissons_ratio','epsc');
export_fig fig25 f25_poissons_ratio -q101 -painters -nocrop -pdf -png -tiff -eps 

%%
%Strain energy calculation

ym_in=input('Enter YM for v0 calculation (MPa): ');
v_in=input('Enter Poissons Ratio for v0 calculation: ');

c1= (1-v_in^2)/(2*ym_in);
c2=(1+v_in)/ym_in;

v0=c1.*(master2(:,9).^2+master2(:,8).^2)+c2.*(master2(:,10).^2+v_in.*master2(:,9).*master2(:,8));  %v0 in MPa

strain_energy=(trapz(sample_area,v0))*10^6;    %Pa
strain_energy_ratio=strain_energy/ie;

output0=[output0, ym_in, v_in, strain_energy, strain_energy_ratio];
output=[output, ym_in, v_in, strain_energy, strain_energy_ratio];
dlmwrite('output0.dat',output,'delimiter','\t');
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
saveas(gcf,'Mean_and_Diff','eps');
% export_fig fig3 f3_bond_dist -pdf -png -q101

%%
fig4=figure (4);
plot(master2(:,10),master2(:,11))
xlabel('Mean Stress');ylabel('1/2*Differential Stress');
set(gcf,'Color','w');
saveas(gcf,'Mean_vs_Diff','tiffn');
%%
for i=1:1:n2-n1+1
    outfile(i,1)=master2(i,13);     %strain
    outfile(i,2)=master2(i,11);     %mean
    outfile(i,3)=master2(i,10).*2;  %differential
    outfile(i,4)=newoutput(i,8);    %porosity
end

fig5=figure(5);
plot(outfile(:,4),outfile(:,2),'linewidth',2);
hold on;
plot(outfile(:,4),outfile(:,3),'r','linewidth',2);
hold on;
legend('Mean vs phi','Differential vs phi');
xlabel('Porosity'); ylabel('Stress (MPa)');
dlmwrite('mdp.dat',outfile,'delimiter','\t');
set(gcf,'Color','w');

% saveas(gcf,'MD_vs_phi','tiffn');
saveas(gcf,'MD_vs_phi','epsc');
saveas(gcf,'MD_vs_phi','fig');
export_fig fig5 f5_MD_phi -q101 -painters -nocrop -pdf -jpeg -eps 

%%
%Plotting spatial distributions of broken bonds

mkdir damage_figs
mkdir damage_figs_uncorrected

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
fig6=figure(6);
aviobj=avifile('Rock_Damage.avi','compression','None','fps',8);

cm=jet(n2-n1+1);

[lla,~]=size(master1);

for i=1:lla     %modified for 5.8a 
    if mastc(i)==1
        positions(i,8)=1;           %1 for shear, 0 for tensile
    end
end

for i=n1:n2
    [pos_x,~]=find(i==positions(:,1));
    axis([min_x3 max_x3 min_y3 max_y3]); pbaspect([a_ratio3 1 1]); box on;
    set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'));
    set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'));
    title(strcat('Strain: ',num2str(master2(i-n1+1,13))));  %this changed
    hold on
    
    if pos_x~=0
        pos_pm=positions(pos_x,:);
        pos_spos=find(pos_pm(:,8)==1);
        pos_tpos=find(pos_pm(:,8)==0);
        pos_pm_s=pos_pm(pos_spos,:);
        pos_pm_t=pos_pm(pos_tpos,:);
        
        if pos_tpos~=0      %tensile
            plot(((pos_pm_t(:,3)+pos_pm_t(:,5))./2),((pos_pm_t(:,4)+pos_pm_t(:,6))./2),'x','color',cm(i-n1+1,:),'markersize',2.5,'linewidth',0.3); %this changed
            hold on
        end
        if pos_spos~=0  %shear
            plot(((pos_pm_s(:,3)+pos_pm_s(:,5))./2),((pos_pm_s(:,4)+pos_pm_s(:,6))./2),'o','color',cm(i-n1+1,:),'markersize',3,'linewidth',0.7);   %this chnaged
            hold on
        end
    end
    
    F=getframe(fig6);
    aviobj=addframe(aviobj,F);
    
    if rem(i,2)==0
        fname=strcat(pwd,'\damage_figs_uncorrected\','damage_',num2str(i));
        saveas(gcf,fname,'epsc');
        saveas(gcf, fname, 'jpg');
    end
    
end

caxis([0 axial_strain]);
colorbar('SouthOutside');
F=getframe(fig6);
aviobj=addframe(aviobj,F);
xlabel('Sample length (m)'); ylabel('Sample width (m)');
xlabel('Sample length (m)');ylabel('Sample width (m)');
set(gcf,'Color','w');
saveas(gcf,'Damage_distribution','tiffn');
saveas(gcf,'Damage_distribution','fig');
aviobj=close(aviobj);

%%
fig7=open('Damage_distribution.fig');
% oldf=gcf;
% newf=copyobj(oldf);
title('');
xlabel('');
ylabel('');
colorbar off;
axis tight
set(gcf,'Color','w');
saveas(gcf,'Damage_distribution_clean','epsc');
export_fig fig7 f7_damage_dist_clean -q101 -painters -nocrop -pdf -jpeg -eps 

%%
%Plot corrected positions

bdata=importdata('S000.OUT');
bdata=bdata(:,1:4);

positions2=zeros(lla,8);
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

for i=1:lla
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
fig8=figure(8);
aviobj_c=avifile('Rock_Damage_Corrected.avi','compression','None','fps',8);

cm=jet(n2-n1+1);
for i=n1:n2
    [pos_x,~]=find(i==positions2(:,1));
    axis([min_x6 max_x6 min_y6 max_y6]); pbaspect([a_ratio6 1 1]); box on;
    set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'));
    set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'));
    title(strcat('Strain: ',num2str(master2(i-n1+1,13))),'fontsize',12,'fontweight','bold');
    hold on
    if pos_x~=0
        pos_pm=positions2(pos_x,:);
        pos_spos=find(pos_pm(:,8)==1);
        pos_tpos=find(pos_pm(:,8)==0);
        pos_pm_s=pos_pm(pos_spos,:);
        pos_pm_t=pos_pm(pos_tpos,:);
        
        if pos_tpos~=0      %tensile
            plot(((pos_pm_t(:,3)+pos_pm_t(:,5))./2),((pos_pm_t(:,4)+pos_pm_t(:,6))./2),'x','color',cm(i-n1+1,:),'markersize',2.5,'linewidth',0.3);
            hold on
        end
        if pos_spos~=0  %shear
            plot(((pos_pm_s(:,3)+pos_pm_s(:,5))./2),((pos_pm_s(:,4)+pos_pm_s(:,6))./2),'o','color',cm(i-n1+1,:),'markersize',3,'linewidth',0.7);
            hold on
            
        end
    end
    
    F=getframe(fig8);
    aviobj_c=addframe(aviobj_c,F);
    
    if rem(i,2)==0
        fname=strcat(pwd,'\damage_figs\','damage_',num2str(i));
        saveas(gcf,fname,'epsc');
        saveas(gcf, fname, 'jpg');
    end  
end


caxis([0 axial_strain]);
colorbar('SouthOutside');
% title(hcb,'Strain','Location','South');
F=getframe(fig8);
aviobj_c=addframe(aviobj_c,F);
xlabel('Sample length (m)','fontweight','bold'); ylabel('Sample width (m)','fontweight','bold');
axis tight
box on
set(gcf,'Color','w');
saveas(gcf,'Damage_distribution_Corrected','epsc');
saveas(gcf,'Damage_distribution_Corrected','fig');
export_fig fig8 f8_damage_dist_corr -q101 -painters -nocrop -pdf -jpeg -eps 
aviobj_c=close(aviobj_c);
%%
fig9=open('Damage_distribution_Corrected.fig');
% oldf=gcf;
% newf=copyobj(oldf);
title('');
xlabel('');
ylabel('');
colorbar off;
axis tight
set(gcf,'Color','w');
saveas(gcf,'Damage_distribution_Corrected_clean','epsc');
export_fig fig9 f9_damage_dist_corr_clean -q101 -painters -nocrop -pdf -jpeg -eps 
%%
%p-value
% disp(strcat('Cycles at S1 max; ',num2str(cycles_at_smax)));
% disp(strcat('Axial at S1 max; ',num2str(strain_at_smax)));
% p_s1=input('Enter Strain at Peak S1: ');
% inp_index=inp_index+1; user_input(inp_index)=p_s1;
% dlmwrite('user_input.dat',user_input,'delimiter','\n');
% p_under=master2(master2(:,13)<=p_s1,[2,13]);
% p_over=master2(master2(:,13)>p_s1,[2,13]);
% p_matrix=master2(:,[4,13]);
% 
% pind_und=find(p_under(:,1)==0);
% p_under(pind_und,:)=[];
% pind_over=find(p_over(:,1)==0);
% p_over(pind_over,:)=[];
% 
% [l_under,~]=size(p_under);
% [l_over,~]=size(p_over);
% 
% fig16=figure(16);
% 
% if isempty(p_over)==1
%     disp('ALL FORESHOCKS!');
%     for i=1:l_under
%         p_under(i,3)=p_s1-p_under(i,2);
%     end
%     loglog(p_under(:,3),p_under(:,1),'bx-','linewidth',3); hold on
%     disp('Interpret FORESHOCKS');
%     p_vector=p_under(:,3);
%     py_vector=p_under(:,1);
%     
%     p_satisfied=0;
%     while p_satisfied==0
%         disp('Select x-range for Constant Slope, lower limit first');
%         [p_x_range,~]=ginput(2);
%         
%         [p_locs1]=find(p_vector(:)>=p_x_range(1));
%         [p_locs2]=find(p_vector(:)<=p_x_range(2));
%         [p_locs]=intersect(p_locs1,p_locs2);
%         
%         p_xpoints=p_vector(p_locs);
%         p_ypoints=py_vector(p_locs);
%         
%         logp_xpoints=log10(p_xpoints);
%         logp_ypoints=log10(p_ypoints);
%         
%         p_linefit=polyfit(logp_xpoints,logp_ypoints,1);
%         p_lineval=polyval(p_linefit,logp_xpoints);
%         p_yy=10.^p_lineval;
%         p_lineplot=loglog(p_xpoints,p_yy,'r','linewidth',3);
%         
%         p_value=p_linefit(1);
%         disp(strcat('p-value: ',num2str(p_value)));
%         p_satisfied=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
%         if p_satisfied==0
%             delete(p_lineplot);
%         end
%     end
%     
%     xlabel('Axial strain difference to peak rupture'); hold on
%     ylabel('Event Rate, dN/dstrain'); hold on
%     legend('Foreshocks','Best fit','location','southwest');
%     
%     output=[output, p_value, 0];
% 
% elseif isempty(p_under)==1
%     
%     disp('ALL AFTERSHOCKS!');
%     for i=1:l_over
%         p_over(i,3)=abs(p_s1-p_over(i,2));
%     end
%     loglog(p_over(:,3),p_over(:,1),'go-','linewidth',2); hold on;
%     disp('Interpret AFTERSHOCKS');
%     
%     p2_vector=p_over(:,3);
%     py2_vector=p_over(:,1);
%     
%     p2_satisfied=0;
%     while p2_satisfied==0
%         disp('Select x-range for Constant Slope, lower limit first');
%         [p2_x_range,~]=ginput(2);
%         
%         [p2_locs1]=find(p2_vector(:)>=p2_x_range(1));
%         [p2_locs2]=find(p2_vector(:)<=p2_x_range(2));
%         [p2_locs]=intersect(p2_locs1,p2_locs2);
%         
%         p2_xpoints=p2_vector(p2_locs);
%         p2_ypoints=py2_vector(p2_locs);
%         
%         logp2_xpoints=log10(p2_xpoints);
%         logp2_ypoints=log10(p2_ypoints);
%         
%         p2_linefit=polyfit(logp2_xpoints,logp2_ypoints,1);
%         p2_lineval=polyval(p2_linefit,logp2_xpoints);
%         p2_yy=10.^p2_lineval;
%         p2_lineplot=loglog(p2_xpoints,p2_yy,'k','linewidth',2);
%         
%         p2_value=p2_linefit(1);
%         disp(strcat('p-value: ',num2str(p2_value)));
%         p2_satisfied=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
%         if p2_satisfied==0
%             delete(p2_lineplot);
%         end
%     end
% 
%     xlabel('Axial strain difference to peak rupture'); hold on
%     ylabel('Event Rate, dN/dstrain'); hold on
%     legend('Aftershocks','Best fit','location','southwest');
%     output=[output, 0, p2_value];
% 
% else
%     disp('FORESHOCKS and AFTERSHOCKS!');
%     for i=1:l_under
%         p_under(i,3)=p_s1-p_under(i,2);
%     end
%     for i=1:l_over
%         p_over(i,3)=abs(p_s1-p_over(i,2));
%     end
%     
%     loglog(p_under(:,3),p_under(:,1),'bx-','linewidth',2); hold on
%     disp('Interpret FORESHOCKS');
%     p_vector=p_under(:,3);
%     py_vector=p_under(:,1);
%     
%     p_satisfied=0;
%     while p_satisfied==0
%         disp('Select x-range for Constant Slope, lower limit first');
%         [p_x_range,~]=ginput(2);
%         
%         [p_locs1]=find(p_vector(:)>=p_x_range(1));
%         [p_locs2]=find(p_vector(:)<=p_x_range(2));
%         [p_locs]=intersect(p_locs1,p_locs2);
%         
%         p_xpoints=p_vector(p_locs);
%         p_ypoints=py_vector(p_locs);
%         
%         logp_xpoints=log10(p_xpoints);
%         logp_ypoints=log10(p_ypoints);
%         
%         p_linefit=polyfit(logp_xpoints,logp_ypoints,1);
%         p_lineval=polyval(p_linefit,logp_xpoints);
%         p_yy=10.^p_lineval;
%         p_lineplot=loglog(p_xpoints,p_yy,'r','linewidth',2);
%         
%         p_value=p_linefit(1);
%         disp(strcat('p-value: ',num2str(p_value)));
%         p_satisfied=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
%         if p_satisfied==0
%             delete(p_lineplot);
%         end
%     end
%     
%     loglog(p_over(:,3),p_over(:,1),'go-','linewidth',2); hold on;
%     disp('Interpret AFTERSHOCKS');
%     
%     p2_vector=p_over(:,3);
%     py2_vector=p_over(:,1);
%     
%     p2_satisfied=0;
%     while p2_satisfied==0
%         disp('Select x-range for Constant Slope, lower limit first');
%         [p2_x_range,~]=ginput(2);
%         
%         [p2_locs1]=find(p2_vector(:)>=p2_x_range(1));
%         [p2_locs2]=find(p2_vector(:)<=p2_x_range(2));
%         [p2_locs]=intersect(p2_locs1,p2_locs2);
%         
%         p2_xpoints=p2_vector(p2_locs);
%         p2_ypoints=py2_vector(p2_locs);
%         
%         logp2_xpoints=log10(p2_xpoints);
%         logp2_ypoints=log10(p2_ypoints);
%         
%         p2_linefit=polyfit(logp2_xpoints,logp2_ypoints,1);
%         p2_lineval=polyval(p2_linefit,logp2_xpoints);
%         p2_yy=10.^p2_lineval;
%         p2_lineplot=loglog(p2_xpoints,p2_yy,'k','linewidth',3);
%         
%         p2_value=p2_linefit(1);
%         disp(strcat('p-value: ',num2str(p2_value)));
%         p2_satisfied=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
%         if p2_satisfied==0
%             delete(p2_lineplot);
%         end
%     end
%     xlabel('Axial strain difference to peak rupture'); hold on
%     ylabel('Event Rate, dN/dstrain'); hold on
%     legend('Foreshocks','Best fit - foreshocks','Aftershocks','Best fit Aftershocks', 'location','southwest');
%     output=[output, p_value, p2_value];
% end
% 
% saveas(gcf,'p_value','epsc');
% export_fig fig16 f16_p_value -q101 -painters -nocrop -pdf -jpeg -eps 
%%
%Damage Index - Renaud

damage_index=zeros(cycles_at_smax,1);
stress_ratio=zeros(cycles_at_smax,1);
initial_porosity=newoutput(1,8);

for i=n1:cycles_at_smax
    damage_index(i-n1+1)=(newoutput(i,8)-initial_porosity)/(1-initial_porosity);
    stress_ratio(i-n1+1)=master2(i-n1+1)/max_S1;
end

%write output file for Damage Index
damage_index_write=[stress_ratio,damage_index];
dlmwrite('damage_index.dat',damage_index_write,'delimiter','\t');

fig22=figure(22);
plot(stress_ratio,damage_index,'linewidth',2); 
hold on

ylabel('Damage Index (\phi - \phi_{i})/(1-\phi_{i})'); xlabel('\sigma/\sigma_{failure}');
axis tight
set(gcf,'Color','w');

% saveas(gcf,'MD_vs_phi','tiffn');
saveas(gcf,'damage_parameter','epsc');
saveas(gcf,'damage_parameter','fig');
export_fig fig22 f22_damage_parameter -q101 -painters -nocrop -pdf -jpeg -eps 
%%
%Damage Index 2- Renaud modified

damage2_index=zeros(n2-n1+1,1);
as_new=master2(:,13);
initial2_porosity=newoutput(1,8);

for i=n1:n2
    damage2_index(i-n1+1)=(newoutput(i,8)-initial_porosity)/(1-initial_porosity);
end

di2_smax=damage2_index(cycles_at_smax-n1+1);

%write output file for Damage Index
damage2_index_write=[as_new,damage2_index];
dlmwrite('damage2_index.dat',damage2_index_write,'delimiter','\t');

fig29=figure(29);
plot(as_new,damage2_index,'linewidth',2); 
hold on
plot(strain_at_smax,di2_smax,'ro','linewidth',2);
hold on

ylabel('Damage Index (\phi - \phi_{i})/(1-\phi_{i})'); xlabel('Axial Strain');
axis tight
set(gcf,'Color','w');

% saveas(gcf,'MD_vs_phi','tiffn');
saveas(gcf,'damage2_parameter','epsc');
saveas(gcf,'damage2_parameter','fig');
export_fig fig29 f29_damage_parameter -q101 -painters -nocrop -pdf -jpeg -eps 
%%
% D-value
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

%write output file for D-value
dwrite=[r_vector,log_cr];
dlmwrite('dvalue.dat',dwrite,'delimiter','\t');

fig10=figure(10);
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
saveas(gcf,'D_value','epsc');
export_fig fig10 f10_D_value_nonclumped -q101 -painters -nocrop -pdf -jpeg -eps 

output=[output, d_value];
%%
% %D-value tensile
% r_vector=(linspace(0,max_X,100))';   %change number of bins for speed
% cr=zeros(length(r_vector),1);
% n_count=zeros(length(r_vector),1);
% 
% for i=1:length(cr)
%     for j=1:n_tensile
%         for k=1:n_tensile
%             r_len1=sqrt((master1_tensile(j,6)-master1_tensile(k,6))^2+(master1_tensile(j,7)-master1_tensile(k,7))^2);
%             r_len2=sqrt((master1_tensile(j,9)-master1_tensile(k,9))^2+(master1_tensile(j,10)-master1_tensile(k,10))^2);
%             r_len3=sqrt((master1_tensile(j,6)-master1_tensile(k,6))^2+(master1_tensile(j,10)-master1_tensile(k,10))^2);
%             r_len4=sqrt((master1_tensile(j,9)-master1_tensile(k,9))^2+(master1_tensile(j,7)-master1_tensile(k,7))^2);
%             if r_len1 <= r_vector(i) && r_len2 <= r_vector(i) && r_len3 <= r_vector(i) && r_len4 <= r_vector(i) && j~=k
%                 n_count(i)=n_count(i)+1;
%             end
%         end
%     end
%     cr(i)=(2*n_count(i))/(n_tensile*(n_tensile-1));
% end
% 
% log_cr=log10(cr);
% log_r=log10(r_vector);
% 
% fig11=figure(11);
% semilogx(r_vector,log_cr,'LineWidth',3);
% hold on
% axis tight
% 
% r_satisfied=0;
% while r_satisfied==0
%     disp('Select x-range for Constant Slope, lower limit first');
%     [r_x_range,~]=ginput(2);
%     
%     [r_locs1]=find(r_vector(:)>=r_x_range(1));
%     [r_locs2]=find(r_vector(:)<=r_x_range(2));
%     [r_locs]=intersect(r_locs1,r_locs2);
%     r_bfit=log_r(r_locs);
%     r_bfit_plot=r_vector(r_locs);
%     cr_bfit=log_cr(r_locs);
%     
%     r_P=polyfit(r_bfit,cr_bfit,1);
%     r_plotP=polyval(r_P,r_bfit);
%     %     log_r_plot=log10(r_plotP);
%     r_bfp=semilogx(r_bfit_plot,r_plotP,'r--','LineWidth',5);
%     %     r_bfp=semilogx(r_bfit,r_plotP,'r--');
%     hold on
%     legend('Spatial Damage Correlation relationship','Best fit line','Location','SouthEast');
%     hold on
%     %     legend('Magnitude-Frequency relationship','Best fit line');
%     %     hold on
%     xlabel('radius'); ylabel('log(C(r))');
%     title('D-value plot- TENSILE ONLY');
%     
%     d_value_t=r_P(1);
%     disp(strcat('D-value: ',num2str(d_value_t)));
%     r_satisfied=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
%     if r_satisfied==0
%         delete(r_bfp);
%     end
% end
% 
% saveas(gcf,'D_value_tensile','epsc');
% export_fig fig11 f11_D_value_nonclumped_tensile -q101 -painters -nocrop -pdf -jpeg -eps 
% %%
% %D-value shear
% r_vector=(linspace(0,max_X,100))';   %change number of bins for speed
% cr=zeros(length(r_vector),1);
% n_count=zeros(length(r_vector),1);
% 
% for i=1:length(cr)
%     for j=1:n_shear
%         for k=1:n_shear
%             r_len1=sqrt((master1_shear(j,6)-master1_shear(k,6))^2+(master1_shear(j,7)-master1_shear(k,7))^2);
%             r_len2=sqrt((master1_shear(j,9)-master1_shear(k,9))^2+(master1_shear(j,10)-master1_shear(k,10))^2);
%             r_len3=sqrt((master1_shear(j,6)-master1_shear(k,6))^2+(master1_shear(j,10)-master1_shear(k,10))^2);
%             r_len4=sqrt((master1_shear(j,9)-master1_shear(k,9))^2+(master1_shear(j,7)-master1_shear(k,7))^2);
%             if r_len1 <= r_vector(i) && r_len2 <= r_vector(i) && r_len3 <= r_vector(i) && r_len4 <= r_vector(i) && j~=k
%                 n_count(i)=n_count(i)+1;
%             end
%         end
%     end
%     cr(i)=(2*n_count(i))/(n_shear*(n_shear-1));
% end
% 
% log_cr=log10(cr);
% log_r=log10(r_vector);
% 
% fig12=figure(12);
% semilogx(r_vector,log_cr,'LineWidth',3);
% hold on
% axis tight
% 
% r_satisfied=0;
% while r_satisfied==0
%     disp('Select x-range for Constant Slope, lower limit first');
%     [r_x_range,~]=ginput(2);
%     
%     [r_locs1]=find(r_vector(:)>=r_x_range(1));
%     [r_locs2]=find(r_vector(:)<=r_x_range(2));
%     [r_locs]=intersect(r_locs1,r_locs2);
%     r_bfit=log_r(r_locs);
%     r_bfit_plot=r_vector(r_locs);
%     cr_bfit=log_cr(r_locs);
%     
%     r_P=polyfit(r_bfit,cr_bfit,1);
%     r_plotP=polyval(r_P,r_bfit);
%     %     log_r_plot=log10(r_plotP);
%     r_bfp=semilogx(r_bfit_plot,r_plotP,'r--','LineWidth',5);
%     %     r_bfp=semilogx(r_bfit,r_plotP,'r--');
%     hold on
%     legend('Spatial Damage Correlation relationship','Best fit line','Location','SouthEast');
%     hold on
%     %     legend('Magnitude-Frequency relationship','Best fit line');
%     %     hold on
%     xlabel('radius'); ylabel('log(C(r))');
%     title('D-value plot- SHEAR ONLY');
%     
%     d_value_s=r_P(1);
%     disp(strcat('D-value: ',num2str(d_value_s)));
%     r_satisfied=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
%     if r_satisfied==0
%         delete(r_bfp);
%     end
% end
% saveas(gcf,'D_value_shear','epsc');
% export_fig fig12 f12_D_value_nonclumped_shear -q101 -painters -nocrop -pdf -jpeg -eps
% %%
%clustering
disp('CLUMP EVENTS');
disp(strcat('Min Cycle: ',num2str(n1)));
disp(strcat('Max Cycle: ',num2str(n2)));

X_diff=max_X-min_X;

disp(strcat('Min X: ',num2str(min_X)));
disp(strcat('Max X: ',num2str(max_X)));
disp(strcat('X range: ',num2str(X_diff)));
disp(strcat('Max particle radius: ',num2str(max_r)));
disp(strcat('Min particle radius: ',num2str(min_r)));

%clump_events
%1 clump number
%2 number of individual events
%3 cycle
%4 strain
%5 Average X
%6 Average Y
%7 Average Z 
%8 total energy
%9 moment

[len,~]=size(master1);

master1_buff=master1;
ccheck=zeros(len,1);
ccheck(:)=-1;
no_of_clusters=0;
%ccheck 
%0 = beginning of clump
%-1 = nonclumped
%any other number = cycle of clumping 


for i=1:lla
    if ccheck(i)==-1
        no_of_clusters=no_of_clusters+1;
        ccheck(i)=no_of_clusters;
        ind_cluster_count=1;                                                                  %reset individual event count to zero
        
        sumx=master1_buff(i,6)+master1_buff(i,9)/2;                 %reset individual counters
        sumy=master1_buff(i,7)+master1_buff(i,10)/2;
        sumz=master1_buff(i,8)+master1_buff(i,11)/2;
        sum_energy=master1_buff(i,18);
        start_cycle=master1_buff(i,15);
        start_strain=start_cycle.*step_disp;
        
        for j=i:lla %check events within cycle limit
            r_len1_c=sqrt((master1_buff(i,6)-master1_buff(j,6))^2+(master1_buff(i,7)-master1_buff(j,7))^2);
            r_len2_c=sqrt((master1_buff(i,9)-master1_buff(j,9))^2+(master1_buff(i,10)-master1_buff(j,10))^2);
            r_len3_c=sqrt((master1_buff(i,6)-master1_buff(j,6))^2+(master1_buff(i,10)-master1_buff(j,10))^2);
            r_len4_c=sqrt((master1_buff(i,9)-master1_buff(j,9))^2+(master1_buff(i,7)-master1_buff(j,7))^2);
            cdiff=abs(master1_buff(i,15)-master1_buff(j,15));
            if i~=j && cdiff <=c_limit && ccheck(j)==-1 && (r_len1_c<=r_limit || r_len2_c<=r_limit || r_len3_c<=r_limit || r_len4_c<=r_limit) %check within radius
                
                ind_cluster_count=ind_cluster_count+1;   %no of events in each cluster
                ccheck(j)=no_of_clusters;         %make sure counter indicates which clump each event belongs to
                
                sumx=sumx+((master1_buff(j,6)+master1_buff(j,9))/2);
                sumy=sumy+((master1_buff(j,7)+master1_buff(j,10))/2);
                sum_energy=sum_energy+master1_buff(j,18);
            end
        end
        clump_events(no_of_clusters,1)=no_of_clusters;
        clump_events(no_of_clusters,2)=ind_cluster_count;
        clump_events(no_of_clusters,3)= start_cycle;
        clump_events(no_of_clusters,4)= start_strain;
        
        clump_events(no_of_clusters,5)= sumx/ind_cluster_count;
        clump_events(no_of_clusters,6)= sumy/ind_cluster_count;
        clump_events(no_of_clusters,7)= sumz;
        
        clump_events(no_of_clusters,8)= sum_energy;
        clump_events(no_of_clusters,9)=(0.67.*log10(clump_events(no_of_clusters,8)))-2.9;
        
    end
end

disp(strcat('No. of individual microcracks: ',num2str(len)));
disp(strcat('No. of events accounted for during clustering: ',num2str(sum(clump_events(:,2)))));

%%
avg_moment=mean(clump_events(:,9));
max_moment=max(clump_events(:,9));
moment_differential=abs(max(clump_events(:,9)))-abs(min(clump_events(:,9)));
disp(strcat('Mean Moment: ',num2str(avg_moment)));
disp(strcat('Max Moment: ',num2str(max_moment)));
disp(strcat('Max Moment Difference: ',num2str(moment_differential)));

output=[output, avg_moment, max_moment, moment_differential];
%%

%clump_events
%1 clump number
%2 number of individual events
%3 cycle
%4 strain
%5 Average X
%6 Average Y
%7 Average Z 
%8 total energy
%9 moment


%clump_master2
%1=cycle
%2=strain
%3=no of macroruptures
%4=total energy

for i=n1:n2
    macros=find(clump_events(:,3)==i);
    macro_count=length(macros);
    clump_master2(i-n1+1,1)=i;    %cycle
    clump_master2(i-n1+1,2)=clump_master2(i-n1+1,1).*step_disp;
    clump_master2(i-n1+1,3)=macro_count;    %no of macrofractures
    clump_master2(i-n1+1,4)=sum(clump_events(macros,8));      %energy
    
%     [macros2,~]=find(clump_events(:,3)==i;   
end
[c_len,~]=size(clump_events);
%%
fig18=figure(18);
subplot(2,1,1)
f1=bar(clump_master2(:,2),clump_master2(:,3));
set(f1,'FaceColor','b','EdgeColor','k','linewidth',0.25);
hold on
xlabel('Axial Strain');ylabel('Number of Microfractures');
title(strcat('Total Clustered AE Events: ',num2str(c_len)));

subplot(2,1,2)
f2=bar(clump_master2(:,2),clump_master2(:,4));
set(f1,'FaceColor','g','EdgeColor','k','linewidth',0.25);
hold on
xlabel('Axial Strain');ylabel('Energy (Mega Joules)');
title(strcat('Total Energy Released: ',num2str(sum(clump_events(:,8))),' Joules'));

set(gcf,'Color','w');
saveas(gcf,'AE_stats_clumped','epsc');
export_fig fig18 f18_ae_stats_clumped -q101 -painters -nocrop -pdf -jpeg -eps 

%%
%clump_events
%1 clump number
%2 number of individual events
%3 cycle
%4 strain
%5 Average X
%6 Average Y
%7 Average Z 
%8 total energy
%9 moment

cm_std=std(clump_events(:,9));
disp(strcat('Moment Standard Deviation: ',num2str(cm_std)));
cm_var=var(clump_events(:,9));
disp(strcat('Moment Variance: ',num2str(cm_var)));
cm_kur=kurtosis(clump_events(:,9));
disp(strcat('Moment Kurtosis: ',num2str(cm_kur)));

output=[output, cm_std, cm_var, cm_kur];

mmwrite=clump_events(:,[1,4,9]);
fname2=strcat('moment_statistics_',num2str(n1),'_',num2str(n2),'.dat');
dlmwrite(fname2,mmwrite,'delimiter','\t');
%%
%clumped video

%clump_master2
%1=cycle
%2=strain
%3=no of macroruptures
%4=total energy

plot_buff24=input('Ready for Spatial Distribution plot ?: ');
fig19=figure(19);
aviobj24=avifile('Clumped_Damage.avi','compression','None','fps',8);
cm4=jet(n2-n1+1);

for i=n1:n2
    nm=clump_master2(i-n1+1,3);
    axis([min_X max_X min_Y max_Y]); pbaspect([a_ratio6 1 1]); box on;
    title('Cumulative Damage Clustered','fontweight','bold');
    set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'));
    set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'));
    xlabel('Sample length (m)','fontweight','bold'); ylabel('Sample width (m)','fontweight','bold'); hold on
    title(strcat('Strain: ',num2str(master2(i-n1+1,13))),'fontsize',12,'fontweight','bold');
    hold on
    if nm~=0
        ind_ind=find(i==clump_events(:,3));
        ind=sum(clump_events(ind_ind,2));   %no of individual events
        cnos=clump_events(ind_ind,1);       %record macrofracture number
        cpos=zeros(ind,2);                 %create empty matrix
        cpos_no=0;
        for j=1:length(cnos)
            prox=find(cnos(j)==ccheck(:));
            for l=1:length(prox)
                cpos_no=cpos_no+1;
                cpos(cpos_no,1)=(master1(prox(l),6)+master1(prox(l),9))/2;
                cpos(cpos_no,2)=(master1(prox(l),7)+master1(prox(l),10))/2;
            end
        end        
        for k=1:ind
            plot(cpos(k,1),cpos(k,2),'.','color',cm4(i-n1+1,:)); hold on 
        end
    end
    F=getframe(fig19);
    aviobj24=addframe(aviobj24,F);  
end

caxis([0 axial_strain]);
colorbar('SouthOutside');
F=getframe(fig19);
aviobj24=addframe(aviobj24,F);
axis tight
box on
set(gcf,'Color','w');
saveas(gcf,'Damage_distribution_Clustered','epsc');
export_fig fig19 f19_damage_dist_clustered -q101 -painters -nocrop -pdf -jpeg -eps 
aviobj24=close(aviobj24);
%%
%b-value clumped
fig13=figure(13);
hist(clump_events(:,9)); hold on
xlabel('Moment'); ylabel ('Frequency'); hold on
saveas(gcf,'moment_hist','epsc');
export_fig fig13 f13_moment_hist -q101 -painters -nocrop -pdf -jpeg -eps 
%%
AE_vector_c=clump_events(:,9);
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

%write output file for b-value
bwrite=[A_vector_c,N_vector_c];
fnameb=strcat('b_value_',num2str(n1),'_',num2str(n2),'.dat');
dlmwrite(fnameb,bwrite,'delimiter','\t');
keyboard
fig14=figure(14);
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
saveas(gcf,'b_value_clumped','epsc');
export_fig fig14 f14_bvalue -q101 -painters -nocrop -pdf -jpeg -eps 

output=[output, b_value_c];
%%
% D-value

[p1_c,~]=size(clump_events);
r_vector_c=(linspace(0,max_X,100))';
cr_c=zeros(length(r_vector_c),1);
n_count_c=zeros(length(r_vector_c),1);

for i=1:length(cr_c)
    for j=1:p1_c
        for k=1:p1_c
            r_len1_c=sqrt((clump_events(j,5)-clump_events(k,5))^2+(clump_events(j,6)-clump_events(k,6))^2);
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
%use log_cr and log_r for regression

fig15=figure(15);
semilogx(r_vector_c,log_cr_c,'LineWidth',3);
hold on
axis tight

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
export_fig fig15 f15_dvalue -q101 -painters -nocrop -pdf -jpeg -eps 

output=[output, d_value_c];
%%
dlmwrite('output.dat',output,'delimiter','\t');




            
            
        


