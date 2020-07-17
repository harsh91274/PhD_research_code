clc
clf
clear all
close all

path(path,genpath('C:\matlab_adds'));
path(path,genpath('C:\export_fig\update2'));
path(path,genpath('C:\Program Files\Glyph & Cog\XpdfReader-win64'));
path(path,genpath('C:\Program Files\gs\gs9.22\bin'));
path(path,genpath('C:\Program Files\MATLAB\R2010a\toolbox\linspecer'));
path(path,genpath('C:\xpdf\xpdf-tools-win-4.00'));
% path(path,genpath('C:\Users\Harsh\cmu'));
path(path,genpath('C:\xpdf\xpdf-tools-win-4.00\bin32'));
warning('off','all')
%%

[k_data,txt1,raw1]=xlsread('plots_data_wrr6.xlsx','k_plot');
[t_data,txt2,raw2]=xlsread('plots_data_wrr6.xlsx','t_plot');
[k_ratios,txt3,raw3]=xlsread('plots_data_wrr6.xlsx','k_ratios');
[k_mondol,txt4,raw4]=xlsread('plots_data_wrr6.xlsx','k_mondol');
[k_neuzil,txt5,raw5]=xlsread('plots_data_wrr6.xlsx','k_neuzil');
[kf_ratios,txt6,raw6]=xlsread('plots_data_wrr6.xlsx','k_frac');
[crs799,txt7,raw7]=xlsread('CRS799.xlsx','crs799','m5:r954');
[crs015,txt8,raw8]=xlsread('CRS015.xlsx','crs015_final','m5:r1341');

[nm1, txt9, raw9]=xlsread('lb_core_ref2.xlsx','PLOT_DATA', 'a2:b7');
[nm2, txt10, raw10]=xlsread('lb_core_ref2.xlsx','PLOT_DATA', 'c2:d7');

[nm1_r, txt9_r, raw9_r]=xlsread('lb_core_ref2.xlsx','PLOT_DATA2', 'a2:b7');
[nm2_r, txt10_r, raw10_r]=xlsread('lb_core_ref2.xlsx','PLOT_DATA2', 'c2:d7');

[nm1_tv, ~, ~]=xlsread('lb_core_ref2.xlsx','PLOT_DATA3', 'a2:b7');
[nm2_tv, ~, ~]=xlsread('lb_core_ref2.xlsx','PLOT_DATA3', 'c2:d7');

[fito1,gof1]=fit(k_data(1:9,1),k_data(1:9,2),'poly1');  %compaction
[fito2,gof2]=fit(k_data(1:6,3),k_data(1:6,4),'poly1');  %mf
[fito3,gof3]=fit(k_data(1:6,5),k_data(1:6,6),'poly1');  %fracture
[fito4,gof4]=fit(k_data(1:8,7),k_data(1:8,8),'poly1');  %kaolinite
[fito5,gof5]=fit(k_data(1:7,9),k_data(1:7,10),'poly1'); %smectite
[fito6,gof6]=fit(k_data(1:9,11),k_data(1:9,12),'poly1');  %GM M=50
[fito7,gof7]=fit(k_data(1:10,13),k_data(1:10,14),'poly1'); %GM M=20

[fito8,gof8]=fit(k_data(1:8,15),k_data(1:8,16),'poly1'); %GM_OR
[fito9,gof9]=fit(k_data(1:8,17),k_data(1:8,18),'poly1'); %KAO_OR
[fito10,gof10]=fit(k_data(1:7,19),k_data(1:7,20),'poly1'); %SMEC_OR

[fito11,gof11]=fit(nm1(:,1),nm1(:,2),'poly1');  %NM1
[fito12,gof12]=fit(nm2(:,1),nm2(:,2),'poly1'); %NM2

[k_mondol_comp, txt11, raw11]=xlsread('plots_data_wrr3.xlsx','k_mondol_compare');

smec_rsq2=rsquare(k_mondol_comp(1:10,3), k_mondol_comp(1:10,5),false);
kao_rsq=rsquare(k_mondol_comp(1:13,8),k_mondol_comp(1:13,10));
kao_rsq2=rsquare(k_mondol_comp(1:13,8),k_mondol_comp(1:13,10),false);

nm1_rsq=rsquare(k_mondol_comp(1:950,12),k_mondol_comp(1:950,13));
nm1_rsq2=rsquare(k_mondol_comp(1:950,12),k_mondol_comp(1:950,13),false);
nm2_rsq=rsquare(k_mondol_comp(:,15),k_mondol_comp(:,16));
nm2_rsq2=rsquare(k_mondol_comp(:,15),k_mondol_comp(:,16),false);


[k_frac_comp, txt12, raw12]=xlsread('plots_data_wrr3.xlsx','k_fracture');
%%
f4=figure(4);
n_orig=fill(k_neuzil(3:13,2),k_neuzil(3:13,1),[0.91 0.91 0.91],'edgecolor','none'); %1
hold on
fill(k_neuzil(3:5,4),k_neuzil(3:5,3),[0.91 0.91 0.91],'edgecolor','none'); %2
hold on
fill(k_neuzil(3:13,6),k_neuzil(3:13,5),[0.91 0.91 0.91],'edgecolor','none');   %3
hold on
fill(k_neuzil(3:7,8),k_neuzil(3:7,7),[0.91 0.91 0.91],'edgecolor','none');     %4
hold on
fill(k_neuzil(3:17,10), k_neuzil(3:17,9),[0.91 0.91 0.91],'edgecolor','none'); %5
hold on
fill(k_neuzil(3:14,12), k_neuzil(3:14,11),[0.91 0.91 0.91],'edgecolor','none');   %6
hold on
fill(k_neuzil(3:11,14), k_neuzil(3:11,13),[0.91 0.91 0.91],'edgecolor','none');    %7
hold on
fill(k_neuzil(3:13,16), k_neuzil(3:13,15),[0.91 0.91 0.91],'edgecolor','none');    %8
hold on
fill(k_neuzil(3:7,18), k_neuzil(3:7,17),[0.91 0.91 0.91],'edgecolor','none');      %9
hold on
fill(k_neuzil(3:7,20), k_neuzil(3:7,19),[0.91 0.91 0.91],'edgecolor','none');  %10
hold on
fill(k_neuzil(3:11,22), k_neuzil(3:11,21),[0.91 0.91 0.91],'edgecolor','none');  %11
hold on
fill(k_neuzil(3:20,24), k_neuzil(3:20,23),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
%natural data plot
fill(k_neuzil(3:7,26), k_neuzil(3:7,25),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:10,28), k_neuzil(3:10,27),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:7,30), k_neuzil(3:7,29),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:7,32), k_neuzil(3:7,31),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
plot(k_neuzil(3:8,34), k_neuzil(3:8,33),'color',[0.91 0.91 0.91],'linewidth',2);
hold on

crs799_plot=plot(crs799(:,1), crs799(:,6),'color', [0, 0.75, 0.75],'linewidth',2);
hold on 
crs015_plot=plot(crs015(:,1), crs015(:,6),'color',[0.9290, 0.6940, 0.1250] , 'linewidth',2);
hold on

nm1_fit=plot(nm1(:,1), fito11(nm1(:,1)),':','linewidth',2); hold on;
nm2_fit=plot(nm2(:,1), fito12(nm2(:,1)),'r:','linewidth',2); hold on;
nm1_trend=plot(nm1(:,1),nm1(:,2),'o'); hold on; 
nm2_trend=plot(nm2(:,1),nm2(:,2),'o'); hold on;
set(nm1_trend,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor','b'); hold on;
set(nm2_trend,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor','r'); hold on;

nm_legend=legend([n_orig, crs799_plot, nm1_trend, nm1_fit, crs015_plot, nm2_trend, nm2_fit], 'Datasets from Neuzil [1994]', 'CRS799 data [Long et al., 2008]', 'NM1 Model data', 'k_{v}^{NM1} trend', 'CRS015B data [Long et al., 2008]', 'NM2 Model data', 'k_{v}^{NM2} trend', 'location','southeast');
set(nm_legend,'FontName','Times New Roman'); hold on;
set(nm_legend, 'FontSize',7.5); hold on; 
legend boxon

y4=ylabel('Log Permeability (log( k_{v} ))'); x4=xlabel ('Porosity ( {\phi} )'); hold on
ylim([-23 -14]); xlim([0 0.9]);
set(gcf,'Color','w');
export_fig f4 wrr9_f4_k_het -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
f1=figure(1);
set(gcf,'Color','w');
n_orig=fill(k_neuzil(3:13,2),k_neuzil(3:13,1),[0.91 0.91 0.91],'edgecolor','none'); %1
hold on
fill(k_neuzil(3:5,4),k_neuzil(3:5,3),[0.91 0.91 0.91],'edgecolor','none'); %2
hold on
fill(k_neuzil(3:13,6),k_neuzil(3:13,5),[0.91 0.91 0.91],'edgecolor','none');   %3
hold on
fill(k_neuzil(3:7,8),k_neuzil(3:7,7),[0.91 0.91 0.91],'edgecolor','none');     %4
hold on
fill(k_neuzil(3:17,10), k_neuzil(3:17,9),[0.91 0.91 0.91],'edgecolor','none'); %5
hold on
fill(k_neuzil(3:14,12), k_neuzil(3:14,11),[0.91 0.91 0.91],'edgecolor','none');   %6
hold on
fill(k_neuzil(3:11,14), k_neuzil(3:11,13),[0.91 0.91 0.91],'edgecolor','none');    %7
hold on
fill(k_neuzil(3:13,16), k_neuzil(3:13,15),[0.91 0.91 0.91],'edgecolor','none');    %8
hold on
fill(k_neuzil(3:7,18), k_neuzil(3:7,17),[0.91 0.91 0.91],'edgecolor','none');      %9
hold on
fill(k_neuzil(3:7,20), k_neuzil(3:7,19),[0.91 0.91 0.91],'edgecolor','none');  %10
hold on
fill(k_neuzil(3:11,22), k_neuzil(3:11,21),[0.91 0.91 0.91],'edgecolor','none');  %11
hold on
fill(k_neuzil(3:20,24), k_neuzil(3:20,23),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
%natural data plot
fill(k_neuzil(3:7,26), k_neuzil(3:7,25),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:10,28), k_neuzil(3:10,27),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:7,30), k_neuzil(3:7,29),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:7,32), k_neuzil(3:7,31),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
plot(k_neuzil(3:8,34), k_neuzil(3:8,33),'color',[0.91 0.91 0.91],'linewidth',2);
hold on

%mondol plot
mondol_s=plot(k_mondol(:,1), k_mondol(:,2), 'r:','linewidth',1.5);
hold on
mondol_k=plot(k_mondol(:,3), k_mondol(:,4), 'r-.','linewidth',1.5);
hold on
%plot model curves
c_fit=plot (k_data(:,15), fito8(k_data(:,15)),'k', 'linewidth',2); hold on;
kaolinite_fit=plot(k_data(:,17), fito9(k_data(:,17)),'-.','color',[0.1 0.45 0.82],'linewidth',2); hold on;
smectite_fit=plot(k_data(:,19), fito10(k_data(:,19)),':','color',[0.1 0.45 0.82],'linewidth',2); hold on;
hold on

c_trend=plot(k_data(1:end,15),k_data(1:end,16), 'o');
hold on
kaolinite_trend=plot(k_data(1:8,17),k_data(1:8,18),'v');
hold on
smectite_trend=plot(k_data(1:7,19),k_data(1:7,20),'^');
hold on
set(c_trend,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7]);
set(kaolinite_trend,'Marker','^','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.1 0.45 0.82]);
set(smectite_trend,'Marker','v','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.1 0.45 0.82]);
hold on

yl3=ylabel('Log Permeability (log( k_{v} ))'); xl3=xlabel ('Porosity ( {\phi} )'); hold on
ylim([-23 -14]); xlim([0 0.9]);

neuzil_legend=legend([n_orig, mondol_k, kaolinite_trend, kaolinite_fit, mondol_s, smectite_trend, smectite_fit, c_trend, c_fit], 'Datasets from Neuzil [1994]', 'Kaolinite [Mondol et al., 2008]', 'Kaolinite Model Data', 'k_{v}^{kaolinite} trend', 'Smectite [Mondol et al., 2008]', 'Smectite Model Data', 'k_{v}^{smectite} trend', 'Intermediate Mudstone Model Data', 'k_{v}^{int} trend', 'location','southeast');
set(neuzil_legend,'FontName','Times New Roman'); hold on;
set(neuzil_legend, 'FontSize',7); hold on; 
legend boxon
% resize_legend(neuzil_legend,1);
set([xl3, yl3],'FontName','Times New Roman');
set(gcf,'Color','w');
box on
export_fig f1 wrr9_f2_kcomp_neuzil -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
f5=figure(5);
set(gcf,'Color','w');
n_orig=fill(k_neuzil(3:13,2),k_neuzil(3:13,1),[0.91 0.91 0.91],'edgecolor','none'); %1
hold on
fill(k_neuzil(3:5,4),k_neuzil(3:5,3),[0.91 0.91 0.91],'edgecolor','none'); %2
hold on
fill(k_neuzil(3:13,6),k_neuzil(3:13,5),[0.91 0.91 0.91],'edgecolor','none');   %3
hold on
fill(k_neuzil(3:7,8),k_neuzil(3:7,7),[0.91 0.91 0.91],'edgecolor','none');     %4
hold on
fill(k_neuzil(3:17,10), k_neuzil(3:17,9),[0.91 0.91 0.91],'edgecolor','none'); %5
hold on
fill(k_neuzil(3:14,12), k_neuzil(3:14,11),[0.91 0.91 0.91],'edgecolor','none');   %6
hold on
fill(k_neuzil(3:11,14), k_neuzil(3:11,13),[0.91 0.91 0.91],'edgecolor','none');    %7
hold on
fill(k_neuzil(3:13,16), k_neuzil(3:13,15),[0.91 0.91 0.91],'edgecolor','none');    %8
hold on
fill(k_neuzil(3:7,18), k_neuzil(3:7,17),[0.91 0.91 0.91],'edgecolor','none');      %9
hold on
fill(k_neuzil(3:7,20), k_neuzil(3:7,19),[0.91 0.91 0.91],'edgecolor','none');  %10
hold on
fill(k_neuzil(3:11,22), k_neuzil(3:11,21),[0.91 0.91 0.91],'edgecolor','none');  %11
hold on
fill(k_neuzil(3:20,24), k_neuzil(3:20,23),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
%natural data plot
fill(k_neuzil(3:7,26), k_neuzil(3:7,25),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:10,28), k_neuzil(3:10,27),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:7,30), k_neuzil(3:7,29),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:7,32), k_neuzil(3:7,31),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
plot(k_neuzil(3:8,34), k_neuzil(3:8,33),'color',[0.91 0.91 0.91],'linewidth',2);
hold on

c_fit=plot (k_data(:,15), fito8(k_data(:,15)),'k', 'linewidth',2); hold on;
mf_fit=plot (k_data(:,3), fito2(k_data(:,3)),'k--', 'linewidth',2); hold on;
f_fit=plot (k_data(:,5), fito3(k_data(:,5)),'k:', 'linewidth',2); hold on;

c_trend=plot(k_data(1:end,15),k_data(1:end,16), 'o');
hold on
mf_trend=plot(k_data(1:end,3),k_data(1:end,4), 'd');
hold on
f_trend=plot(k_data(1:end,5),k_data(1:end,6), 's');
hold on

set(c_trend,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7]);
set(mf_trend,'Marker','d','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9]);
set(f_trend,'Marker','s','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5]);

yl3=ylabel('Log Permeability (log( k_{v} ))'); xl3=xlabel ('Porosity ( {\phi} )'); hold on
ylim([-23 -14]); xlim([0 0.9]);

neuzil_legend2=legend([n_orig, c_trend, c_fit, mf_trend, mf_fit, f_trend, f_fit], 'Datasets from Neuzil [1994]', 'Intermediate Mudstone Model Data', 'k_{v}^{int} trend', 'Microfracture Network Data', 'k_{v}^{mf} trend', 'Macrofracture Propagation Data','k_{v}^{frac} trend','location','southeast');
set(neuzil_legend2,'FontSize',7.5);
set(neuzil_legend2,'FontName','Times New Roman');
legend boxon
set([xl3, yl3],'FontName','Times New Roman');
set(gcf,'Color','w');
box on
export_fig f2 wrr9_f6_k_fracture -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
% f9=figure(9);
% emf_trend=plot (k_frac_comp(:,3), k_frac_comp(:,2),'k:','linewidth',1.5); hold on;
% ef_trend=plot (k_frac_comp(:,6), k_frac_comp(:,5),'k--','linewidth',1.5); hold on;
% 
% set(emf_trend,'Marker','d','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9]);
% set(ef_trend,'Marker','s','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5]);
% 
% y9=ylabel('Log Permeability (log( k_{v} ))'); x9=xlabel ('Effective Fracture Width ( {\epsilon}_{eff} )'); hold on
% ylim([-23 -12]); 
% f9_l=legend('Microfracture Network Data ( {\epsilon}^{mf}_{eff} )', 'Fracture Propagation Data ( {\epsilon}^{frac}_{eff} )','location','southeast');
% % set(f9_l,'interpreter','tex');
% set([xl3, yl3],'FontName','Times New Roman');
% set(gcf,'Color','w');
% box on
% export_fig f9 wrr3b_f9_k_fracture -q101 -painters -nocrop -pdf -png -tiff -eps 
% % xlim([0 0.9]);
%%
f2=figure(2);

[fit_t1,t1_gof]=fit(t_data(1:9,1),t_data(1:9,2),'poly1');
[fit_t4,t4_gof]=fit(t_data(1:8,7),t_data(1:8,8),'poly1');
[fit_t5,t5_gof]=fit(t_data(1:7,9),t_data(1:7,10),'poly1');

% [fit_t2,t2_gof]=fit(t_data(1:8,3),t_data(1:8,4),'power2');
% [fit_t3,t3_gof]=fit(t_data(1:7,5),t_data(1:7,6),'exp2'); 

t_mf=plot(t_data(1:end,3),t_data(1:end,4),'k--','linewidth',1);
hold on
t_f=plot(t_data(1:end,5),t_data(1:end,6),'k:','linewidth',1);
hold on
t_comp=plot(t_data(1:end,15),t_data(1:end,16),'k-','linewidth',1);
hold on
% t_kaolinite=plot(t_data(1:8,7),t_data(1:8,8),'k-.','linewidth',1);
% hold on
% t_smectite=plot(t_data(1:7,9),t_data(1:7,10),'k:','linewidth',1);
% hold on

% t_gm50=plot(t_data(1:end,11),t_data(1:end,12),'k--','linewidth',1);
hold on
% t_gm20=plot(t_data(1:end,13),t_data(1:end,14),'k:','linewidth',1);
hold on

set(gcf,'Color','w');
% set(t_gm50,'Marker','p','Markersize',6, 'MarkerEdgeColor','k');
% set(t_gm20,'Marker','*','Markersize',6, 'MarkerEdgeColor','k');
set(t_mf,'Marker','d','Markersize',6, 'MarkerFaceColor',[0.9 0.9 0.9], 'MarkerEdgeColor','k');
set(t_f,'Marker','s','Markersize',6,'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerEdgeColor','k');
set(t_comp,'Marker','o','Markersize',6,'MarkerFaceColor',[0.7 0.7 0.7], 'MarkerEdgeColor','k');
% set(t_kaolinite,'Marker','^','Markersize',6,'MarkerFaceColor',[0.55 0.57 0.67], 'MarkerEdgeColor','k');
% set(t_smectite,'Marker','v','Markersize',6,'MarkerFaceColor',[0.55 0.57 0.67], 'MarkerEdgeColor','k');

% t_f_fit=plot (t_data(:,5), fit_t3(t_data(:,5)),'k:', 'linewidth',2); hold on;
% t_mf_fit=plot (t_data(:,3), fit_t2(t_data(:,3)),'k--', 'linewidth',2); hold on;
% t_c_fit=plot (t_data(:,1), fit_t1(t_data(:,1)),'k', 'linewidth',2); hold on;
% t_kaolinite_fit=plot (t_data(1:8,7), fit_t4(t_data(1:8,7)),'color',[0 0.19 0.56] ,'linewidth',2); hold on;
% t_smectite_fit=plot (t_data(1:7,9), fit_t5(t_data(1:7,9)),':','color',[0 0.19 0.56] ,'linewidth',2); hold on;

xl2=xlabel('Porosity ( {\phi} )'); yl2=ylabel ('Vertical Tortuosity ( {\tau}_{v} )'); hold on
xlim([0 0.9]); 
% ylim([0.8 2.6]);

t_legend=legend([t_comp, t_mf, t_f], 'Intermediate Mudstone Model', 'Microfracture Network Model', 'Macrofracture Model','location','northeast'); 
set(t_legend,'FontSize',8);
set(t_legend,'FontName','Times New Roman');
% set(t_legend,'position',[ 0.6256 0.3 0.2696 0.22]); 
set([xl2, yl2],'FontName','Times New Roman');
set(gcf,'Color','w');
export_fig f2 wrr6b_fS2c_t_frac -pdf -png -eps -tiff -q101 -nocrop

%%
f7=figure(7);

% [fit_t1,t1_gof]=fit(t_data(1:9,1),t_data(1:9,2),'poly1');
% [fit_t4,t4_gof]=fit(t_data(1:8,7),t_data(1:8,8),'poly1');
% [fit_t5,t5_gof]=fit(t_data(1:7,9),t_data(1:7,10),'poly1');
[fit_t6,t6_gof]=fit(t_data(1:9,11),t_data(1:9,12),'poly1');
[fit_t7,t7_gof]=fit(t_data(1:10,13),t_data(1:10,14),'poly1');

% [fit_t2,t2_gof]=fit(t_data(1:8,3),t_data(1:8,4),'power2');
% [fit_t3,t3_gof]=fit(t_data(1:7,5),t_data(1:7,6),'exp2'); 

% t_mf=plot(t_data(1:end,3),t_data(1:end,4),'k--','linewidth',1);
% hold on
% t_f=plot(t_data(1:end,5),t_data(1:end,6),'k:','linewidth',1);
% hold on
t_comp=plot(t_data(1:end,15),t_data(1:end,16),'k-','linewidth',1);
hold on
t_kaolinite=plot(t_data(1:8,17),t_data(1:8,18),'k-.','linewidth',1);
hold on
t_smectite=plot(t_data(1:7,19),t_data(1:7,20),'k:','linewidth',1);
hold on
% t_mf=plot(t_data(1:end,3),t_data(1:end,4),'k--','linewidth',1);
% hold on
% t_f=plot(t_data(1:end,5),t_data(1:end,6),'k:','linewidth',1);
% hold on

%plot rotated
% t_comp_or=plot(t_data(1:8,15),t_data(1:8,16),'color',[0.6350, 0.0780, 0.1840],'linewidth',1);
% hold on
% t_kao_or=plot(t_data(1:8,17),t_data(1:8,18),'-.','color',[0.6350, 0.0780, 0.1840],'linewidth',1);
% hold on
% t_smec_or=plot(t_data(1:7,19),t_data(1:7,20),':','color',[0.6350, 0.0780, 0.1840],'linewidth',1);
% hold on
%linestyles
set(t_comp,'Marker','o','Markersize',6,'MarkerFaceColor',[0.7 0.7 0.7], 'MarkerEdgeColor','k'); hold on; 
set(t_kaolinite,'Marker','^','Markersize',6,'MarkerFaceColor',[0.55 0.57 0.67], 'MarkerEdgeColor','k'); hold on; 
set(t_smectite,'Marker','v','Markersize',6,'MarkerFaceColor',[0.55 0.57 0.67], 'MarkerEdgeColor','k'); hold on;
% set(t_mf,'Marker','d','Markersize',6, 'MarkerFaceColor',[0.9 0.9 0.9], 'MarkerEdgeColor','k');
% set(t_f,'Marker','s','Markersize',6,'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerEdgeColor','k');

nm1_tv_plot=plot(nm1_tv(:,1), nm1_tv(:,2),'b:','linewidth',2); hold on;
nm2_tv_plot=plot(nm2_tv(:,1), nm2_tv(:,2),'r:','linewidth',2); hold on;
set(nm1_tv_plot,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor','b'); hold on;
set(nm2_tv_plot,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor','r'); hold on;

%rotated linestyles
% set(t_comp_or,'Marker','o','Markersize',6,'MarkerFaceColor',[0.6350, 0.0780, 0.1840], 'MarkerEdgeColor','k'); hold on;
% set(t_kao_or,'Marker','^','Markersize',6,'MarkerFaceColor',[0.6350, 0.0780, 0.1840], 'MarkerEdgeColor','k'); hold on; 
% set(t_smec_or,'Marker','v','Markersize',6,'MarkerFaceColor',[0.6350, 0.0780, 0.1840], 'MarkerEdgeColor','k'); hold on;

xl7=xlabel('Porosity ( {\phi} )'); yl7=ylabel ('Vertical Tortuosity ( {\tau}_{v} )'); hold on
xlim([0 0.9]); 
% ylim([0.8 2.6]);

t7_legend2=legend('Kaolinite Model', 'Intermediate Mudstone Model', 'Smectite Model','NM1 Model','NM2 Model', 'location','northeast'); 
set(t7_legend2,'FontSize',8);
set(t7_legend2,'FontName','Times New Roman');
legend boxon
% set(t_legend,'position',[ 0.6256 0.3 0.2696 0.22]); 
set([xl7, yl7],'FontName','Times New Roman');
set(gcf,'Color','w');
export_fig f7 wrr9_fS2_torutosity -pdf -png -eps -tiff -q101 -nocrop
keyboard
% print(gcf,'f4_tortuosity.png','-dpng','-r300'); 
% export_fig f1 f4_tortuosity -pdf -png -tiff -q101
%%
f31=figure(31);

% mf_ratio=semilogy(k_ratios(1:end,3),k_ratios(1:end,4),'k--','linewidth',1);
% hold on
% f_ratio=semilogy(k_ratios(1:end,5),k_ratios(1:end,6),'k:','linewidth',1);
% hold on
% c_ratio=plot(k_ratios(1:end,1),k_ratios(1:end,2),'k-','linewidth',1);
% hold on
% kaolinite_ratio=plot(k_ratios(1:8,7),k_ratios(1:8,8),'k-.','linewidth',1);
% hold on
% smectite_ratio=plot(k_ratios(1:7,9),k_ratios(1:7,10),'k:','linewidth',1);
% hold on
%plot_rotated
c_or_ratio=semilogy(k_ratios(1:end,11),k_ratios(1:end,12),'k-', 'linewidth',1);
hold on
kao_or_ratio=semilogy(k_ratios(1:8,13),k_ratios(1:8,14),'k-.','linewidth',1);
hold on
smec_or_ratio=semilogy(k_ratios(1:7,15),k_ratios(1:7,16),'k:','linewidth',1);
hold on

%
set(gcf,'Color','w');
% set(mf_ratio,'Marker','d','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9]); hold on;
% set(f_ratio,'Marker','s','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5]); hold on; 
% set(c_ratio,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7]);hold on;
% set(kaolinite_ratio,'Marker','^','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.55 0.57 0.67]); hold on;
% set(smectite_ratio,'Marker','v','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.55 0.57 0.67]); hold on;


set(c_or_ratio,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',	[0.7 0.7 0.7]);hold on;
set(kao_or_ratio,'Marker','^','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.55 0.57 0.67]); hold on;
set(smec_or_ratio,'Marker','v','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.55 0.57 0.67]); hold on;

hold on
xl3=xlabel('Porosity ( {\phi} )'); yl3=ylabel ('Permeability Anisotropy ( k_{h}/k_{v} )'); hold on
xlim([0 0.9]); 
% ylim([0.00001 100]);

rat_legend=legend([kao_or_ratio c_or_ratio smec_or_ratio],'Kaolinite Model', 'Intermediate Mudstone Model','Smectite Model','location','Northeast');
% set(rat_legend,'position',[0.493137254901961 0.214228617106315 0.376470588235294 0.235011990407674]); 
set(rat_legend,'FontSize',7);
set(rat_legend,'FontName','Times New Roman');
set([xl3, yl3],'FontName','Times New Roman');
set(gcf,'Color','w');
% saveas(gcf,'d6_f5_ratios','eps');
export_fig f31 wrr6b_f3_krat_homo -pdf -png -eps -tiff -q101 -nocrop
%%
f32=figure(32);

% mf_ratio=semilogy(k_ratios(1:end,3),k_ratios(1:end,4),'k--','linewidth',1);
% hold on
% f_ratio=semilogy(k_ratios(1:end,5),k_ratios(1:end,6),'k:','linewidth',1);
% hold on
% c_ratio=plot(k_ratios(1:end,1),k_ratios(1:end,2),'k-','linewidth',1);
% hold on
% kaolinite_ratio=plot(k_ratios(1:8,7),k_ratios(1:8,8),'k-.','linewidth',1);
% hold on
% smectite_ratio=plot(k_ratios(1:7,9),k_ratios(1:7,10),'k:','linewidth',1);
% hold on
%plot_rotated
c_or_ratio=semilogy(k_ratios(1:end,11),k_ratios(1:end,12),'k-', 'linewidth',1);
hold on
kao_or_ratio=semilogy(k_ratios(1:8,13),k_ratios(1:8,14),'k-.','linewidth',1);
hold on
smec_or_ratio=semilogy(k_ratios(1:7,15),k_ratios(1:7,16),'k:','linewidth',1);
hold on


nm1_r_plot=semilogy(nm1_r(:,1), nm1_r(:,2),'b:','linewidth',2); hold on;
nm2_r_plot=semilogy(nm2_r(:,1), nm2_r(:,2),'r:','linewidth',2); hold on;
% nm1_r_trend=plot(nm1_r(:,1),nm1_r(:,2),':'); hold on; 
% nm2_r_trend=plot(nm2_r(:,1),nm2_r(:,2),'r:'); hold on;
set(nm1_r_plot,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor','b'); hold on;
set(nm2_r_plot,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor','r'); hold on;

%
set(gcf,'Color','w');
% set(mf_ratio,'Marker','d','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9]); hold on;
% set(f_ratio,'Marker','s','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5]); hold on; 
% set(c_ratio,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7]);hold on;
% set(kaolinite_ratio,'Marker','^','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.55 0.57 0.67]); hold on;
% set(smectite_ratio,'Marker','v','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.55 0.57 0.67]); hold on;


set(c_or_ratio,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',	[0.7 0.7 0.7]);hold on;
set(kao_or_ratio,'Marker','^','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.55 0.57 0.67]); hold on;
set(smec_or_ratio,'Marker','v','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.55 0.57 0.67]); hold on;

hold on
xl4=xlabel('Porosity ( {\phi} )'); yl4=ylabel ('Permeability Anisotropy ( k_{h}/k_{v} )'); hold on
xlim([0 0.9]); 
% ylim([0.00001 100]);

rat2_legend=legend('Kaolinite Model', 'Intermediate Mudstone Model','Smectite Model','NM1 Model', 'NM2 Model','location','Northeast');
% set(rat_legend,'position',[0.493137254901961 0.214228617106315 0.376470588235294 0.235011990407674]); 
set(rat2_legend,'FontSize',7);
set(rat2_legend,'FontName','Times New Roman');
set([xl4, yl4],'FontName','Times New Roman');
set(gcf,'Color','w');
% saveas(gcf,'d6_f5_ratios','eps');
export_fig f32 wrr6b_f3_krat_all -pdf -png -eps -tiff -q101 -nocrop
%%
f33=figure(33);

xl5=xlabel('Porosity ( {\phi} )'); yl5=ylabel ('Vertical Tortuosity ( {\tau}_{v} )'); hold on
xlim([0 0.9]); 
ylim([0 60]);

t2_legend=legend('NM1 Model', 'NM2 Model','location','northeast'); 
set(t2_legend,'FontSize',8);
set(t2_legend,'FontName','Times New Roman');
% set(t_legend,'position',[ 0.6256 0.3 0.2696 0.22]); 
set([xl5, yl5],'FontName','Times New Roman');
set(gcf,'Color','w');
export_fig f33 wrr6b_fS2b_t_het -pdf -png -eps -tiff -q101 -nocrop
%%
f39=figure(39);
emf_trend=plot (k_frac_comp(1:6,3), k_frac_comp(1:6,2),'k:','linewidth',1.5); hold on;
ef_trend=plot (k_frac_comp(1:6,6), k_frac_comp(1:6,5),'k--','linewidth',1.5); hold on;

set(emf_trend,'Marker','d','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9]);
set(ef_trend,'Marker','s','Markersize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5]);

y9=ylabel('Log Permeability (log( k_{v} ))'); x9=xlabel ('Effective Fracture Width ( {\epsilon}_{eff} )'); hold on
ylim([-23 -12]); 
f9_l=legend('Microfracture Network Data ( {\epsilon}^{mf}_{eff} )', 'Fracture Propagation Data ( {\epsilon}^{frac}_{eff} )','location','southeast');
% set(f9_l,'interpreter','tex');
set([xl3, yl3],'FontName','Times New Roman');
set(gcf,'Color','w');
box on
export_fig f39 wrr6b_f39_eff_kfrac -q101 -painters -nocrop -pdf -png -tiff -eps 
% xlim([0 0.9]);

%%
f11=figure(11);
set(gcf,'Color','w');
n_orig=fill(k_neuzil(3:13,2),k_neuzil(3:13,1),[0.91 0.91 0.91],'edgecolor','none'); %1
hold on
fill(k_neuzil(3:5,4),k_neuzil(3:5,3),[0.91 0.91 0.91],'edgecolor','none'); %2
hold on
fill(k_neuzil(3:13,6),k_neuzil(3:13,5),[0.91 0.91 0.91],'edgecolor','none');   %3
hold on
fill(k_neuzil(3:7,8),k_neuzil(3:7,7),[0.91 0.91 0.91],'edgecolor','none');     %4
hold on
fill(k_neuzil(3:17,10), k_neuzil(3:17,9),[0.91 0.91 0.91],'edgecolor','none'); %5
hold on
fill(k_neuzil(3:14,12), k_neuzil(3:14,11),[0.91 0.91 0.91],'edgecolor','none');   %6
hold on
fill(k_neuzil(3:11,14), k_neuzil(3:11,13),[0.91 0.91 0.91],'edgecolor','none');    %7
hold on
fill(k_neuzil(3:13,16), k_neuzil(3:13,15),[0.91 0.91 0.91],'edgecolor','none');    %8
hold on
fill(k_neuzil(3:7,18), k_neuzil(3:7,17),[0.91 0.91 0.91],'edgecolor','none');      %9
hold on
fill(k_neuzil(3:7,20), k_neuzil(3:7,19),[0.91 0.91 0.91],'edgecolor','none');  %10
hold on
fill(k_neuzil(3:11,22), k_neuzil(3:11,21),[0.91 0.91 0.91],'edgecolor','none');  %11
hold on
fill(k_neuzil(3:20,24), k_neuzil(3:20,23),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
%natural data plot
fill(k_neuzil(3:7,26), k_neuzil(3:7,25),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:10,28), k_neuzil(3:10,27),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:7,30), k_neuzil(3:7,29),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:7,32), k_neuzil(3:7,31),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
plot(k_neuzil(3:8,34), k_neuzil(3:8,33),'color',[0.91 0.91 0.91],'linewidth',2);
hold on

%mondol plot
mondol_s=plot(k_mondol(:,1), k_mondol(:,2), 'r:','linewidth',1.5);
hold on
mondol_k=plot(k_mondol(:,3), k_mondol(:,4), 'r-.','linewidth',1.5);
hold on

yl3=ylabel('Log Permeability (log( k_{v} ))'); xl3=xlabel ('Porosity ( {\phi} )'); hold on
ylim([-23 -12]); xlim([0 0.9]);

neuzil_legend=legend([n_orig, mondol_k, mondol_s], 'Datasets from Neuzil, 1994', 'Kaolinite (Mondol et al., 2008)', 'Smectite (Mondol et al., 2008)', 'location','southeast');
set(neuzil_legend,'FontName','Times New Roman'); hold on;
set(neuzil_legend, 'FontSize',6); hold on; 
legend boxon
% resize_legend(neuzil_legend,1);
set([xl3, yl3],'FontName','Times New Roman');
set(gcf,'Color','w');
box on
export_fig f11 wrr6b_f11_k_exponly -q101 -painters -nocrop -pdf -png -tiff -eps 

%%
f41=figure(41);
n_orig=fill(k_neuzil(3:13,2),k_neuzil(3:13,1),[0.91 0.91 0.91],'edgecolor','none'); %1
hold on
fill(k_neuzil(3:5,4),k_neuzil(3:5,3),[0.91 0.91 0.91],'edgecolor','none'); %2
hold on
fill(k_neuzil(3:13,6),k_neuzil(3:13,5),[0.91 0.91 0.91],'edgecolor','none');   %3
hold on
fill(k_neuzil(3:7,8),k_neuzil(3:7,7),[0.91 0.91 0.91],'edgecolor','none');     %4
hold on
fill(k_neuzil(3:17,10), k_neuzil(3:17,9),[0.91 0.91 0.91],'edgecolor','none'); %5
hold on
fill(k_neuzil(3:14,12), k_neuzil(3:14,11),[0.91 0.91 0.91],'edgecolor','none');   %6
hold on
fill(k_neuzil(3:11,14), k_neuzil(3:11,13),[0.91 0.91 0.91],'edgecolor','none');    %7
hold on
fill(k_neuzil(3:13,16), k_neuzil(3:13,15),[0.91 0.91 0.91],'edgecolor','none');    %8
hold on
fill(k_neuzil(3:7,18), k_neuzil(3:7,17),[0.91 0.91 0.91],'edgecolor','none');      %9
hold on
fill(k_neuzil(3:7,20), k_neuzil(3:7,19),[0.91 0.91 0.91],'edgecolor','none');  %10
hold on
fill(k_neuzil(3:11,22), k_neuzil(3:11,21),[0.91 0.91 0.91],'edgecolor','none');  %11
hold on
fill(k_neuzil(3:20,24), k_neuzil(3:20,23),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
%natural data plot
fill(k_neuzil(3:7,26), k_neuzil(3:7,25),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:10,28), k_neuzil(3:10,27),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:7,30), k_neuzil(3:7,29),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:7,32), k_neuzil(3:7,31),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
plot(k_neuzil(3:8,34), k_neuzil(3:8,33),'color',[0.91 0.91 0.91],'linewidth',2);
hold on

crs799_plot=plot(crs799(:,1), crs799(:,6),'color', [0, 0.75, 0.75],'linewidth',2);
hold on 
crs015_plot=plot(crs015(:,1), crs015(:,6),'color',[0.9290, 0.6940, 0.1250] , 'linewidth',2);
hold on

nm_legend=legend([n_orig, crs799_plot, crs015_plot], 'Datasets from Neuzil, 1994', 'CRS799 data (Long et al., 2008)', 'CRS015B data (Long et al., 2008)', 'location','southeast');
set(nm_legend,'FontName','Times New Roman'); hold on;
set(nm_legend, 'FontSize',6); hold on; 
legend boxon

y4=ylabel('Log Permeability (log( k_{v} ))'); x4=xlabel ('Porosity ( {\phi} )'); hold on
ylim([-23 -12]); xlim([0 0.9]);
set(gcf,'Color','w');
export_fig f41 wrr6b_f41_k_het_exponly -q101 -painters -nocrop -pdf -png -tiff -eps 

%%
f12=figure(12);
set(gcf,'Color','w');
n_orig=fill(k_neuzil(3:13,2),k_neuzil(3:13,1),[0.91 0.91 0.91],'edgecolor','none'); %1
hold on
fill(k_neuzil(3:5,4),k_neuzil(3:5,3),[0.91 0.91 0.91],'edgecolor','none'); %2
hold on
fill(k_neuzil(3:13,6),k_neuzil(3:13,5),[0.91 0.91 0.91],'edgecolor','none');   %3
hold on
fill(k_neuzil(3:7,8),k_neuzil(3:7,7),[0.91 0.91 0.91],'edgecolor','none');     %4
hold on
fill(k_neuzil(3:17,10), k_neuzil(3:17,9),[0.91 0.91 0.91],'edgecolor','none'); %5
hold on
fill(k_neuzil(3:14,12), k_neuzil(3:14,11),[0.91 0.91 0.91],'edgecolor','none');   %6
hold on
fill(k_neuzil(3:11,14), k_neuzil(3:11,13),[0.91 0.91 0.91],'edgecolor','none');    %7
hold on
fill(k_neuzil(3:13,16), k_neuzil(3:13,15),[0.91 0.91 0.91],'edgecolor','none');    %8
hold on
fill(k_neuzil(3:7,18), k_neuzil(3:7,17),[0.91 0.91 0.91],'edgecolor','none');      %9
hold on
fill(k_neuzil(3:7,20), k_neuzil(3:7,19),[0.91 0.91 0.91],'edgecolor','none');  %10
hold on
fill(k_neuzil(3:11,22), k_neuzil(3:11,21),[0.91 0.91 0.91],'edgecolor','none');  %11
hold on
fill(k_neuzil(3:20,24), k_neuzil(3:20,23),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
%natural data plot
fill(k_neuzil(3:7,26), k_neuzil(3:7,25),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:10,28), k_neuzil(3:10,27),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:7,30), k_neuzil(3:7,29),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
fill(k_neuzil(3:7,32), k_neuzil(3:7,31),[0.91 0.91 0.91],'edgecolor','none');   %12
hold on
plot(k_neuzil(3:8,34), k_neuzil(3:8,33),'color',[0.91 0.91 0.91],'linewidth',2);
hold on


%plot model curves
c_fit=plot (k_data(:,15), fito8(k_data(:,15)),'k', 'linewidth',2); hold on;
c_trend=plot(k_data(1:end,15),k_data(1:end,16), 'o');hold on
set(c_trend,'Marker','o','Markersize',6,'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7]); hold on

yl3=ylabel('Log Permeability (log( k_{v} ))'); xl3=xlabel ('Porosity ( {\phi} )'); hold on
ylim([-23 -12]); xlim([0 0.9]);

neuzil_legend=legend([n_orig, c_trend, c_fit], 'Datasets from Neuzil, 1994', 'Intermediate Mudstone Model Data', 'k_{v}^{int} trend', 'location','southeast');
set(neuzil_legend,'FontName','Times New Roman'); hold on;
set(neuzil_legend, 'FontSize',6); hold on; 
legend boxon
% resize_legend(neuzil_legend,1);
set([xl3, yl3],'FontName','Times New Roman');
set(gcf,'Color','w');
box on
export_fig f12 wrr6b_f12_kcomp_neuzil_int -q101 -painters -nocrop -pdf -png -tiff -eps 