clc
clf
clear all
close all

path(path,genpath('C:\export_fig\update'));
path(path,genpath('C:\Program Files\Glyph & Cog\XpdfReader-win64'));
path(path,genpath('C:\Program Files\gs\gs9.22\bin'));
path(path,genpath('C:\Program Files\MATLAB\R2010a\toolbox\linspecer'));
path(path,genpath('C:\xpdf\xpdf-tools-win-4.00'));
path(path,genpath('C:\xpdf\xpdf-tools-win-4.00\bin32'));
warning('off','all')

[num_data,txt1,raw1]=xlsread('berea_ldb_data.xlsx','number_mf');
[shear_f,txt2,raw2]=xlsread('berea_ldb_data.xlsx','shear_fraction');
[energy,txt3,raw3]=xlsread('berea_ldb_data.xlsx','energy');
[moment,txt4,raw4]=xlsread('berea_ldb_data.xlsx','moment');
[bvalue,txt5,raw5]=xlsread('berea_ldb_data.xlsx','b_value');
[dvalue,txt6,raw6]=xlsread('berea_ldb_data.xlsx','d_value');
[bd_litt,txt7,raw7]=xlsread('berea_ldb_data.xlsx','bd');
bd_data=xlsread('berea_ldb_data.xlsx','data','ad11:ae30');

cp_data=xlsread('berea_ldb_data.xlsx','data','i11:i20');

berea_s=xlsread('berea_ldb_data.xlsx','data','j11:j20');
berea_t=xlsread('berea_ldb_data.xlsx','data','k11:k20');
ldb_s=xlsread('berea_ldb_data.xlsx','data','j21:j30');
ldb_t=xlsread('berea_ldb_data.xlsx','data','k21:k30');

berea_energy_s=xlsread('berea_ldb_data.xlsx','data','n11:n20');
berea_energy_t=xlsread('berea_ldb_data.xlsx','data','o11:o20');
ldb_energy_s=xlsread('berea_ldb_data.xlsx','data','n21:n30');
ldb_energy_t=xlsread('berea_ldb_data.xlsx','data','o21:o30');

b15_data=xlsread('berea_ldb_data.xlsx','bplots15');
d15_data=xlsread('berea_ldb_data.xlsx','Dplots15');


cp_proxy=[0 10 20 30 40 50 60 70 80 90];

[berea_s1,txt8,raw8]=xlsread('berea_ldb_data.xlsx','berea_splots');
[ldb_s1,txt9,raw9]=xlsread('berea_ldb_data.xlsx','ldb_splots');

[berea_di,txt10,raw10]=xlsread('berea_ldb_data.xlsx','di_berea');
[ldb_di,txt11,raw11]=xlsread('berea_ldb_data.xlsx','di_ldb');

%circle datapoints, dotted lines - berea
%square datapoints, solid lines - ldb

%%
%number of microfractures
f1=figure(1);
num_berea=plot(num_data(:,1), num_data(:,2),'k:','linewidth',2);
hold on
num_ldb=plot(num_data(:,1), num_data(:,3),'k','linewidth',2);
hold on
set(num_berea,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5]);
set(num_ldb,'Marker','s','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9]);
hold on

f1_x=xlabel('Confining Pressure (MPa)'); f1_y=ylabel ('Number of Microfractures'); ylim([2000 6500]); hold on
f1_legend=legend('Berea Sandstone','Lac Du Bonnet Granite','location','southeast');
set(f1_legend,'FontSize',9);
set(f1_legend,'FontName','Helvetica');
set([f1_x f1_y],'FontName','Helvetica');
% grid on
set(gcf,'Color','w'); 

saveas(gcf,'number_mf','mmat');
export_fig f1 f1_number -pdf -png -eps -tiff -q101 -nocrop
%%
%energy
f2=figure(2);
e_berea=plot(energy(:,1), energy(:,2),'k:','linewidth',2);
hold on
e_ldb=plot(energy(:,1), energy(:,3),'k','linewidth',2);
hold on
set(e_berea,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5]);
set(e_ldb,'Marker','s','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9]);
hold on

f2_x=xlabel('Confining Pressure (MPa)'); f2_y=ylabel ('Energy (Joules)'); hold on
% xlim([-23 -12]); ylim([0 0.9]);

f2_legend=legend('Berea Sandstone','Lac du Bonnet Granite','location','northwest');

set(f2_legend,'FontSize',9);
set(f2_legend,'FontName','Helvetica');
set([f2_x f2_y],'FontName','Helvetica');
set(gcf,'Color','w');

saveas(gcf,'energy','mmat');
export_fig f2 f2_energy -pdf -png -eps -tiff -q101 -nocrop

%%
f3=figure(3);
% xlim([-23 -12]); 
m_berea=plot(moment(:,1), moment(:,2),'k:','linewidth',2);
hold on
m_ldb=plot(moment(:,1), moment(:,3),'k','linewidth',2);
hold on
set(m_berea,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5]);
set(m_ldb,'Marker','s','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9]);
hold on

% ylim([-4.2 -2.2]);
hold on
f3_x=xlabel('Confining Pressure (MPa)'); f3_y=ylabel ('Largest Moment Magnitude'); hold on

set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'));
f3_legend=legend('Berea Sandstone', 'Lac du Bonnet Granite','location','southeast');
set(f3_legend,'FontSize',9);
set(f3_legend,'FontName','Helvetica');
set([f3_x f3_y],'FontName','Helvetica');
% grid on
set(gcf,'Color','w');

saveas(gcf,'moment','mmat');
export_fig f3 f3_moment -pdf -png -eps -tiff -q101 -nocrop

%%
f4=figure(4);
b_berea=plot(bvalue(:,1), bvalue(:,2),'k:','linewidth',2);
hold on
b_ldb=plot(bvalue(:,1), bvalue(:,3),'k','linewidth',2);
hold on
set(b_berea,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5]);
set(b_ldb,'Marker','s','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9]);
hold on

set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'));
f4_x=xlabel('Confining Pressure (MPa)'); f4_y=ylabel ('b-value'); hold on
f4_legend=legend('Berea Sandstone', 'Lac du Bonnet Granite','location','northeast');
set(f4_legend,'FontSize',9);
set(f4_legend,'FontName','Helvetica');
set([f4_x f4_y],'FontName','Helvetica');
% grid on
set(gcf,'Color','w');

saveas(gcf,'bvalue','mmat');
export_fig f4 f4_bvalue -pdf -png -eps -tiff -q101 -nocrop
%%
f5=figure(5);
d_berea=plot(dvalue(:,1), dvalue(:,2),'k:','linewidth',2);
hold on
d_ldb=plot(dvalue(:,1), dvalue(:,3),'k','linewidth',2);
hold on
set(d_berea,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5]);
set(d_ldb,'Marker','s','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9]);
hold on

set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'));
f5_x=xlabel('Confining Pressure (MPa)'); f5_y=ylabel ('D-value'); hold on
f5_legend=legend('Berea Sandstone', 'Lac du Bonnet Granite','location','southeast');
set(f5_legend,'FontSize',9);
set(f5_legend,'FontName','Helvetica');
set([f5_x f5_y],'FontName','Helvetica');
set(gcf,'Color','w');
saveas(gcf,'Dvalue','mmat');
export_fig f5 f5_dvalue -pdf -png -eps -tiff -q101 -nocrop
%%
f6=figure(6);
f_cmap2=linspecer(6);
% 
xlim([0.9 1.8]);
ylim([0.5 3.5]);
hold on
bd_berea=plot(bd_data(1:10,1), bd_data(1:10,2),'o','MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5]);
hold on
bd_ldb=plot(bd_data(11:20,1), bd_data(11:20,2),'s','MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9]);
hold on
[bd_fit,gof]=fit(bd_data(:,1),bd_data(:,2),'poly1');
bd_fitplot=plot(bd_data(:,1),bd_fit(bd_data(:,1)),'k','linewidth',2); hold on

bd_1=plot(bd_litt(1:2,1), bd_litt(1:2,2),'color',f_cmap2(1,:),'linewidth',2);
hold on
bd_2=plot(bd_litt(3:4,1), bd_litt(3:4,2),':','color',f_cmap2(2,:),'linewidth',2);
hold on
bd_3=plot(bd_litt(5:6,1), bd_litt(5:6,2),':','color',f_cmap2(3,:),'linewidth',2);
hold on
bd_4=plot(bd_litt(7:8,1), bd_litt(7:8,2),':','color',f_cmap2(4,:),'linewidth',2);
hold on
bd_5=plot(bd_litt(9:10,1), bd_litt(9:10,2),':','color',f_cmap2(5,:),'linewidth',2);
hold on
% bd_6=plot(bd_litt(11:12,1), bd_litt(11:12,2),'color',f_cmap2(6,:),'linewidth',2);
% hold on

% ylim([-4.2 -2.2]);

set(gca,'yticklabel',num2str(get(gca,'ytick')','%.2f'));
set(gca,'xticklabel',num2str(get(gca,'xtick')','%.2f'));

f6_x=xlabel('b-value'); f6_y=ylabel ('D-value'); hold on
% xlim([-23 -12]); ylim([0 0.9]);

f6_legend=legend([bd_berea bd_ldb bd_fitplot bd_1 bd_2 bd_3 bd_4 bd_5 ],{'Berea Sandstone data','Lac du Bonnet Granite data', 'Best fit line', 'Amitrano, 2003','Wyss et al., 2004', 'Hirata, 1986', 'Oncel et al., 1996', 'Aki, 1981'},'location','northeast');
set(f6_legend,'FontSize',9);
set(f6_legend,'FontName','Helvetica');
set([f6_x f6_y],'FontName','Helvetica');
% text(0.55, 0.5,'D = -0.737 b + 2.419','units','normalized','font'); 

text(0.95, 2.1,'D = -0.94 b + 2.73','fontname','helvetica','fontsize',15,'fontweight','bold');
hold on
box on
set(gcf,'Color','w');

saveas(gcf,'bd','mmat');
export_fig f6 f6_bd -pdf -png -eps -tiff -q101 -nocrop
%%
%Number bar plot - Berea
f7=figure(7);
f_cmap3=linspecer(2);
bb1=bar(cp_proxy,[berea_s berea_t],0.9,'stack');colormap(f_cmap3); 
set(bb1,'EdgeColor','k','linewidth',1);
hold on
set(gca, 'xtickLabel', cp_data); hold on
axis tight; ylim([0 6500]);
f7_x=xlabel('Confining Pressure (MPa)'); f7_y=ylabel ('Number of Microfractures'); hold on
f7_legend=legend('Shear', 'Tensile','location','northeast');
set(f7_legend,'FontSize',9);
set(f7_legend,'FontName','Helvetica');
set([f7_x f7_y],'FontName','Helvetica');
set(gcf,'Color','w');
saveas(gcf,'berea_cp_mf','mmat');
export_fig f7 f7_berea_cpmf -painters -pdf -png -eps -tiff -q101 -nocrop 
%%
%number bar plots - LDB
f8=figure(8);
bb2=bar(cp_proxy,[ldb_s ldb_t],0.9,'stack');colormap(f_cmap3);
set(bb2,'EdgeColor','k','linewidth',1);
hold on
set(gca, 'xtickLabel', cp_data); hold on
axis tight; ylim([0 4000]);
f8_x=xlabel('Confining Pressure (MPa)'); f8_y=ylabel ('Number of Microfractures'); hold on
f8_legend=legend('Shear', 'Tensile','location','northwest');
set(f8_legend,'FontSize',9);
set(f8_legend,'FontName','Helvetica');
set([f8_x f8_y],'FontName','Helvetica');
set(gcf,'Color','w');
saveas(gcf,'ldb_cp_mf','mmat');
export_fig f8 f8_ldb_cpmf -painters -pdf -png -eps -tiff -q101 -nocrop
%%
%berea-splots
f_cmap3=linspecer(10);
f9=figure(9);
berea_cp0=plot(berea_s1(1:100,1), berea_s1(1:100,2),'color',f_cmap3(1,:),'linewidth',2); hold on
berea_cp2=plot(berea_s1(1:100,1), berea_s1(1:100,3),'color',f_cmap3(2,:),'linewidth',2); hold on
berea_cp5=plot(berea_s1(1:100,1), berea_s1(1:100,4),'color',f_cmap3(3,:),'linewidth',2); hold on
berea_cp10=plot(berea_s1(1:100,1), berea_s1(1:100,5),'color',f_cmap3(4,:),'linewidth',2); hold on
berea_cp15=plot(berea_s1(1:100,1), berea_s1(1:100,6),'color',f_cmap3(5,:),'linewidth',2); hold on
berea_cp20=plot(berea_s1(1:100,1), berea_s1(1:100,7),'color',f_cmap3(6,:),'linewidth',2); hold on
berea_cp25=plot(berea_s1(1:100,1), berea_s1(1:100,8),'color',f_cmap3(7,:),'linewidth',2); hold on
berea_cp30=plot(berea_s1(1:100,1), berea_s1(1:100,9),'color',f_cmap3(8,:),'linewidth',2); hold on
berea_cp40=plot(berea_s1(1:100,1), berea_s1(1:100,10),'color',f_cmap3(9,:),'linewidth',2); hold on
berea_cp50=plot(berea_s1(1:100,1), berea_s1(1:100,11),'color',f_cmap3(10,:),'linewidth',2); hold on

xlim([0 berea_s1(100,1)]);
f9_x=xlabel('Axial Strain'); f9_y=ylabel ('Axial Stress (MPa)'); hold on
f9_legend=legend('0 MPa','2 MPa','5 MPa','10 MPa','15 MPa','20 MPa','25 MPa','30 MPa','40 MPa','50 MPa','location','northwest');
set(f9_legend,'FontSize',9);
set(f9_legend,'FontName','Helvetica');
set([f9_x f9_y],'FontName','Helvetica');
set(gcf,'Color','w');
saveas(gcf,'berea_s1','mmat');
export_fig f9 f9_berea_s1 -pdf -png -eps -tiff -q101 -nocrop
%%
%ldb_splots
f10=figure(10);
ldb_cp0=plot(ldb_s1(1:100,1), ldb_s1(1:100,2),'color',f_cmap3(1,:),'linewidth',2); hold on
ldb_cp2=plot(ldb_s1(1:100,1), ldb_s1(1:100,3),'color',f_cmap3(2,:),'linewidth',2); hold on
ldb_cp5=plot(ldb_s1(1:100,1), ldb_s1(1:100,4),'color',f_cmap3(3,:),'linewidth',2); hold on
ldb_cp10=plot(ldb_s1(1:100,1), ldb_s1(1:100,5),'color',f_cmap3(4,:),'linewidth',2); hold on
ldb_cp15=plot(ldb_s1(1:100,1), ldb_s1(1:100,6),'color',f_cmap3(5,:),'linewidth',2); hold on
ldb_cp20=plot(ldb_s1(1:100,1), ldb_s1(1:100,7),'color',f_cmap3(6,:),'linewidth',2); hold on
ldb_cp25=plot(ldb_s1(1:100,1), ldb_s1(1:100,8),'color',f_cmap3(7,:),'linewidth',2); hold on
ldb_cp30=plot(ldb_s1(1:100,1), ldb_s1(1:100,9),'color',f_cmap3(8,:),'linewidth',2); hold on
ldb_cp40=plot(ldb_s1(1:100,1), ldb_s1(1:100,10),'color',f_cmap3(9,:),'linewidth',2); hold on
ldb_cp50=plot(ldb_s1(1:100,1), ldb_s1(1:100,11),'color',f_cmap3(10,:),'linewidth',2); hold on

xlim([0 ldb_s1(100,1)]);
f10_x=xlabel('Axial Strain'); f10_y=ylabel ('Axial Stress (MPa)'); hold on
f10_legend=legend('0 MPa','2 MPa','5 MPa','10 MPa','15 MPa','20 MPa','25 MPa','30 MPa','40 MPa','50 MPa','location','northeast');
set(f10_legend,'FontSize',9);
set(f10_legend,'FontName','Helvetica');
set([f10_x f10_y],'FontName','Helvetica');
set(gcf,'Color','w');
saveas(gcf,'ldb_s1','mmat');
export_fig f10 f10_ldb_s1 -pdf -png -eps -tiff -q101 -nocrop

%%
%unconfined splots
f_cmap5=linspecer(2);
f11=figure(11);
berea_cp=plot(berea_s1(1:100,1), berea_s1(1:100,2),'color',f_cmap5(1,:),'linewidth',2); hold on
ldb_cp=plot(ldb_s1(1:100,1), ldb_s1(1:100,2),'color',f_cmap5(2,:),'linewidth',2); hold on

xlim([0 ldb_s1(100,1)]);
f11_x=xlabel('Axial Strain'); f11_y=ylabel ('Axial Stress (MPa)'); hold on
f11_legend=legend('Berea Sandstone','Lac du Bonnet Granite','location','northeast');
set(f11_legend,'FontSize',9);
set(f11_legend,'FontName','Helvetica');
set([f11_x f11_y],'FontName','Helvetica');
set(gcf,'Color','w');
saveas(gcf,'splot_unconfined','mmat');
export_fig f11 f11_splot_unconfined -pdf -png -eps -tiff -q101 -nocrop
%%
%energy bar plot - Berea
f12=figure(12);
f_cmap3=linspecer(2);
e_bb1=bar(cp_proxy,[berea_energy_s berea_energy_t],0.9,'stack');colormap(f_cmap3); 
set(e_bb1,'EdgeColor','k','linewidth',1);
hold on
set(gca, 'xtickLabel', cp_data); hold on
axis tight; 
ylim([0 9]);
f12_x=xlabel('Confining Pressure (MPa)'); f12_y=ylabel ('Fracture Energy (Joules)'); hold on
f12_legend=legend('Shear', 'Tensile','location','northwest');
set(f12_legend,'FontSize',9);
set(f12_legend,'FontName','Helvetica');
set([f12_x f12_y],'FontName','Helvetica');
set(gcf,'Color','w');
saveas(gcf,'berea_cp_energy','mmat');
export_fig f12 f12_berea_cp_energy -painters -pdf -png -eps -tiff -q101 -nocrop 
%%
%energy bar plots - LDB
f13=figure(13);
e_bb2=bar(cp_proxy,[ldb_energy_s ldb_energy_t],0.9,'stack');colormap(f_cmap3);
set(e_bb2,'EdgeColor','k','linewidth',1);
hold on
set(gca, 'xtickLabel', cp_data); hold on
axis tight; 
ylim([0 35]);
f13_x=xlabel('Confining Pressure (MPa)'); f13_y=ylabel ('Fracture Energy (Joules)'); hold on
f13_legend=legend('Shear', 'Tensile','location','northwest');
set(f13_legend,'FontSize',9);
set(f13_legend,'FontName','Helvetica');
set([f13_x f13_y],'FontName','Helvetica');
set(gcf,'Color','w');
saveas(gcf,'ldb_cp_energy','mmat');
export_fig f13 f13_ldb_cp_energy -painters -pdf -png -eps -tiff -q101 -nocrop
%%
%b-values @ 15 mpa plot
f14=figure(14);
berea_b15=plot(b15_data(:,1), b15_data(:,2),'o','color',f_cmap5(1,:),'linewidth',1.2); hold on
ldb_b15=plot(b15_data(:,3), b15_data(:,4),'o','color',f_cmap5(2,:),'linewidth',1.2); hold on

% xlim([0 ldb_s1(100,1)]);
axis tight
f14_x=xlabel('Event Magnitude (M_{e})'); f14_y=ylabel ('Cumulative Number of Events, Log N (>M_{e})'); hold on
f14_legend=legend('Berea Sandstone','Lac du Bonnet Granite','location','southwest');
set(f14_legend,'FontSize',9);
set(f14_legend,'FontName','Helvetica');
set([f14_x f14_y],'FontName','Helvetica');
set(gcf,'Color','w');
saveas(gcf,'bplot_cp15','mmat');
export_fig f14 f14_bplot_cp15 -pdf -png -eps -tiff -q101 -nocrop
%%
%D-values @ 15 mpa plot
f15=figure(15);
% berea_d15=semilogx(d15_data(:,1), d15_data(:,2),'o','color',f_cmap5(1,:),'linewidth',1.2); hold on
berea_d15=semilogx(d15_data(:,5), d15_data(:,6),'o','color',f_cmap5(1,:),'linewidth',1.2); hold on
ldb_d15=semilogx(d15_data(:,3), d15_data(:,4),'o','color',f_cmap5(2,:),'linewidth',1.2); hold on
% xlim([0 ldb_s1(100,1)]);
axis tight
f15_x=xlabel('log (C(r))'); f15_y=ylabel ('Radius (r)'); hold on
f15_legend=legend('Berea Sandstone','Lac du Bonnet Granite','location','southeast');
set(f15_legend,'FontSize',9);
set(f15_legend,'FontName','Helvetica');
set([f15_x f15_y],'FontName','Helvetica');
set(gcf,'Color','w');
saveas(gcf,'dplot_cp15','mmat');
export_fig f15 f15_dplot_cp15 -pdf -png -eps -tiff -q101 -nocrop
%%
%berea-DIplots
f_cmap3=linspecer(10);
f16=figure(16);
berea_di0=plot(berea_di(:,1), berea_di(:,2),'color',f_cmap3(1,:),'linewidth',2); hold on
berea_di2=plot(berea_di(:,3), berea_di(:,4),'color',f_cmap3(2,:),'linewidth',2); hold on
berea_di5=plot(berea_di(:,5), berea_di(:,6),'color',f_cmap3(3,:),'linewidth',2); hold on
berea_di10=plot(berea_di(:,7), berea_di(:,8),'color',f_cmap3(4,:),'linewidth',2); hold on
berea_di15=plot(berea_di(:,9), berea_di(:,10),'color',f_cmap3(5,:),'linewidth',2); hold on
berea_di20=plot(berea_di(:,11), berea_di(:,12),'color',f_cmap3(6,:),'linewidth',2); hold on
berea_di25=plot(berea_di(:,13), berea_di(:,14),'color',f_cmap3(7,:),'linewidth',2); hold on
berea_di30=plot(berea_di(:,15), berea_di(:,16),'color',f_cmap3(8,:),'linewidth',2); hold on
berea_di40=plot(berea_di(:,17), berea_di(:,18),'color',f_cmap3(9,:),'linewidth',2); hold on
berea_di50=plot(berea_di(:,19), berea_di(:,20),'color',f_cmap3(10,:),'linewidth',2); hold on

% xlim([0 berea_s1(100,1)]);
f16_x=xlabel('\sigma/\sigma_{failure}'); 
f16_y=ylabel('Damage Index (\phi - \phi_{i})/(1-\phi_{i})'); 
hold on
f16_legend=legend('0 MPa','2 MPa','5 MPa','10 MPa','15 MPa','20 MPa','25 MPa','30 MPa','40 MPa','50 MPa','location','southwest');
axis tight;
set(f16_legend,'FontSize',9);
set(f16_legend,'FontName','Helvetica');
set([f16_x f16_y],'FontName','Helvetica');
set(gcf,'Color','w');
saveas(gcf,'berea_DI','mmat');
export_fig f16 f16_berea_DI -pdf -png -eps -tiff -q101 -nocrop
%%
%ldb - DIplots

f17=figure(17);
ldb_di0=plot(ldb_di(:,1), ldb_di(:,2),'color',f_cmap3(1,:),'linewidth',2); hold on
ldb_di2=plot(ldb_di(:,3), ldb_di(:,4),'color',f_cmap3(2,:),'linewidth',2); hold on
ldb_di5=plot(ldb_di(:,5), ldb_di(:,6),'color',f_cmap3(3,:),'linewidth',2); hold on
ldb_di10=plot(ldb_di(:,7), ldb_di(:,8),'color',f_cmap3(4,:),'linewidth',2); hold on
ldb_di15=plot(ldb_di(:,9), ldb_di(:,10),'color',f_cmap3(5,:),'linewidth',2); hold on
ldb_di20=plot(ldb_di(:,11), ldb_di(:,12),'color',f_cmap3(6,:),'linewidth',2); hold on
ldb_di25=plot(ldb_di(:,13), ldb_di(:,14),'color',f_cmap3(7,:),'linewidth',2); hold on
ldb_di30=plot(ldb_di(:,15), ldb_di(:,16),'color',f_cmap3(8,:),'linewidth',2); hold on
ldb_di40=plot(ldb_di(:,17), ldb_di(:,18),'color',f_cmap3(9,:),'linewidth',2); hold on
ldb_di50=plot(ldb_di(:,19), ldb_di(:,20),'color',f_cmap3(10,:),'linewidth',2); hold on

% xlim([0 berea_s1(100,1)]);
f17_x=xlabel('\sigma/\sigma_{failure}'); 
f17_y=ylabel('Damage Index (\phi - \phi_{i})/(1-\phi_{i})'); 
hold on
f17_legend=legend('0 MPa','2 MPa','5 MPa','10 MPa','15 MPa','20 MPa','25 MPa','30 MPa','40 MPa','50 MPa','location','northwest');
axis tight;
set(f17_legend,'FontSize',9);
set(f17_legend,'FontName','Helvetica');
set([f17_x f17_y],'FontName','Helvetica');
set(gcf,'Color','w');
saveas(gcf,'ldb_DI','mmat');
export_fig f17 f17_ldb_DI -pdf -png -eps -tiff -q101 -nocrop
%%
f18=figure(18);
berea_di15_c=plot(berea_di(:,9), berea_di(:,10),'color',f_cmap5(1,:),'linewidth',2); hold on
ldb_di15_c=plot(ldb_di(:,9), ldb_di(:,10),'color',f_cmap5(2,:),'linewidth',2); hold on

% xlim([0 ldb_s1(100,1)]);
f18_x=xlabel('\sigma/\sigma_{failure}'); 
f18_y=ylabel('Damage Index (\phi - \phi_{i})/(1-\phi_{i})'); 
hold on
f18_legend=legend('Berea Sandstone','Lac du Bonnet Granite','location','southwest');
axis tight;
set(f18_legend,'FontSize',9);
set(f18_legend,'FontName','Helvetica');
set([f18_x f18_y],'FontName','Helvetica');
set(gcf,'Color','w');
saveas(gcf,'DI_cp15','mmat');
export_fig f18 f18_DI_cp15 -pdf -png -eps -tiff -q101 -nocrop

