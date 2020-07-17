%%
%10 mpa
e_kur10=data10(:,26);
mc_var10=data10(:,21);
e_rate10=data10(:,27);

f41=figure(41);
%microcrack rate and MC variance 
s1=subplot(2,1,1);
yyaxis(s1,'left');
plot(s1,strain10,mc_rate10,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain10,mc_var10,'linewidth',2); hold on; 
box on
% axis tight;
xlim([0 max(strain10)]);
ylabel ('Microcracking Variance', 'fontsize', 8);
% legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%shear fraction and D-value
s2=subplot(2,1,2);
yyaxis(s2,'left');
plot(s2,strain10,shear_fraction10,'linewidth',2); hold on; 
% axis tight; 
ylabel('Shear Microcrack Fraction', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain10,d_value10,'linewidth',2); hold on; 
% axis tight;
box on
xlim([0 max(strain10)]);
ylabel ('D-value', 'fontsize', 8);
% legend(s2,'Microcrack Energy Variance', 'D-value','fontsize',6,'Location','East');

set(gcf,'Color','w'); 
export_fig f41 f41_10mpa_composite1 -q101 -painters -nocrop -pdf -png -tiff -eps 

f42=figure(42);
%microcrack rate and MC variance 
s1=subplot(2,1,1);
yyaxis(s1,'left');
plot(s1,strain10,e_rate10,'linewidth',2); hold on; 
% axis tight; 
ylabel('Energy Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain10,e_var10,'linewidth',2); hold on; 
% axis tight;
box on
xlim([0 max(strain10)]);
ylabel ('Energy Variance', 'fontsize', 8);
% legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%moment difference and b-value
s2=subplot(2,1,2);
yyaxis(s2,'left');
plot(s2,strain10,moment_diff10,'linewidth',2); hold on; 
% axis tight; 
ylabel('Moment Range', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain10,b_value10,'linewidth',2); hold on; 
% axis tight;
box on
xlim([0 max(strain10)]);
ylabel ('b-value', 'fontsize', 8);
% legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f42 f42_10mpa_composite2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%35 mpa composites
e_kur35=data35(:,26);
mc_var35=data35(:,21);
e_rate35=data35(:,27);

f43=figure(43);
%microcrack rate and MC variance 
s1=subplot(2,1,1);
yyaxis(s1,'left');
plot(s1,strain35,mc_rate35,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain35,mc_var35,'linewidth',2); hold on; 
% axis tight;
box on
xlim([0 max(strain35)]);
ylabel ('Microcracking Variance', 'fontsize', 8);
% legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%shear fraction and D-value
s2=subplot(2,1,2);
yyaxis(s2,'left');
plot(s2,strain35,shear_fraction35,'linewidth',2); hold on; 
% axis tight; 
ylabel('Shear Microcrack Fraction', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain35,d_value35,'linewidth',2); hold on; 
% axis tight;
box on
xlim([0 max(strain35)]);
ylabel ('D-value', 'fontsize', 8);
% legend(s2,'Microcrack Energy Variance', 'D-value','fontsize',6,'Location','East');

set(gcf,'Color','w'); 
export_fig f43 f43_35mpa_composite1 -q101 -painters -nocrop -pdf -png -tiff -eps 

f44=figure(44);
%microcrack rate and MC variance 
s1=subplot(2,1,1);
yyaxis(s1,'left');
plot(s1,strain35,e_rate35,'linewidth',2); hold on; 
% axis tight; 
ylabel('Energy Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain35,e_var35,'linewidth',2); hold on; 
% axis tight;
box on
xlim([0 max(strain35)]);
ylabel ('Energy Variance', 'fontsize', 8);
% legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%moment difference and b-value
s2=subplot(2,1,2);
yyaxis(s2,'left');
plot(s2,strain35,moment_diff35,'linewidth',2); hold on; 
% axis tight; 
ylabel('Moment Range', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain35,b_value35,'linewidth',2); hold on; 
% axis tight;
box on
xlim([0 max(strain35)]);
ylabel ('b-value', 'fontsize', 8);
% legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f44 f44_35mpa_composite2 -q101 -painters -nocrop -pdf -png -tiff -eps 



%%
%10 mpa b, D plots
data10b=xlsread('p3_data_binned3.xlsx','10b');
data10d=xlsread('p3_data_binned3.xlsx','10d');
f_cmap3=linspecer(4);

init_bx=data10b(:,1);
init_by=data10b(:,2);
nu_bx=data10b(:,3);
nu_by=data10b(:,4);
rup_bx=data10b(:,5);
rup_by=data10b(:,6);
sld_bx=data10b(:,7);
sld_by=data10b(:,8);

f45=figure(45);
plot(init_bx, init_by, 'color',f_cmap3(1,:),'linewidth',2); hold on;
plot(nu_bx, nu_by, 'color',f_cmap3(2,:),'linewidth',2); hold on;
plot(rup_bx, rup_by, 'color',f_cmap3(3,:),'linewidth',2); hold on;
plot(sld_bx, sld_by, 'color',f_cmap3(4,:),'linewidth',2); hold on;
axis tight; 
xlabel('Moment Magnitude'); ylabel('Frequency');
legend('Axial Strain = 0.0216', 'Axial Strain = 0.0329', 'Axial Strain = 0.0567', 'Axial Strain = 0.0724','location','northeast');
box on
set(gcf,'Color','w'); 
export_fig f45 f45_10mpa_b_value -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%10 mpa d-value
init_dx=data10d(:,1);
init_dy=data10d(:,2);
nu_dx=data10d(:,3);
nu_dy=data10d(:,4);
rup_dx=data10d(:,5);
rup_dy=data10d(:,6);
sld_dx=data10d(:,7);
sld_dy=data10d(:,8);

f46=figure(46);
semilogx(init_dx, init_dy, 'color',f_cmap3(1,:),'linewidth',2); hold on;
semilogx(nu_dx, nu_dy, 'color',f_cmap3(2,:),'linewidth',2); hold on;
semilogx(rup_dx, rup_dy, 'color',f_cmap3(3,:),'linewidth',2); hold on;
semilogx(sld_dx, sld_dy, 'color',f_cmap3(4,:),'linewidth',2); hold on;
axis tight; 
xlabel('Radius (m)'); ylabel('Correlation Integral (C(R))');
legend('Axial Strain = 0.0216', 'Axial Strain = 0.0329', 'Axial Strain = 0.0567', 'Axial Strain = 0.0724','location','northwest');
box on
set(gcf,'Color','w'); 
export_fig f46 f46_10mpa_d_value -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%0 MPa stress plot
f47=figure(47); 
strain0=data0(:,29);
stress0=data0(:,28);
plot(strain0, stress0, 'color',f_cmap3(1,:),'linewidth',2); hold on;
xlabel('Axial Strain'); ylabel('Axial Stress (MPa)');
box on
set(gcf,'Color','w'); 
export_fig f47 f47_0Mpa_stress_strain -q101 -painters -nocrop -pdf -png -tiff -eps 

%%
%25 mpa composites
e_kur25=data25(:,26);
mc_var25=data25(:,21);
e_rate25=data25(:,27);

f51=figure(51);
%microcrack rate and MC variance 
s1=subplot(2,1,1);
yyaxis(s1,'left');
plot(s1,strain25,mc_rate25,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain25,mc_var25,'linewidth',2); hold on; 
% axis tight;
box on
xlim([0 max(strain25)]);
ylabel ('Microcracking Variance', 'fontsize', 8);
% legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%shear fraction and D-value
s2=subplot(2,1,2);
yyaxis(s2,'left');
plot(s2,strain25,shear_fraction25,'linewidth',2); hold on; 
% axis tight; 
ylabel('Shear Microcrack Fraction', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain25,d_value25,'linewidth',2); hold on; 
% axis tight;
box on
xlim([0 max(strain25)]);
ylabel ('D-value', 'fontsize', 8);
% legend(s2,'Microcrack Energy Variance', 'D-value','fontsize',6,'Location','East');

set(gcf,'Color','w'); 
export_fig f51 f51_25mpa_composite1 -q101 -painters -nocrop -pdf -png -tiff -eps 

f52=figure(52);
%microcrack rate and MC variance 
s1=subplot(2,1,1);
yyaxis(s1,'left');
plot(s1,strain25,e_rate25,'linewidth',2); hold on; 
% axis tight; 
ylabel('Energy Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain25,e_var25,'linewidth',2); hold on; 
% axis tight;
box on
xlim([0 max(strain25)]);
ylabel ('Energy Variance', 'fontsize', 8);
% legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%moment difference and b-value
s2=subplot(2,1,2);
yyaxis(s2,'left');
plot(s2,strain25,moment_diff25,'linewidth',2); hold on; 
% axis tight; 
ylabel('Moment Range', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain25,b_value25,'linewidth',2); hold on; 
% axis tight;
box on
xlim([0 max(strain25)]);
ylabel ('b-value', 'fontsize', 8);
% legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f52 f52_25mpa_composite2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
% %individual plots
% f91=figure(91);
% %microcrack rate and shear fraction
% e_kur15=data15(:,26);
% 
% %microcrack energy variance and D-value
% subplot(2,1,1);
% % yyaxis(s2,'left');
% plot(strain15,e_var15,'linewidth',2); hold on; 
% % axis tight; 
% ylabel('Microcrack Energy Variance', 'fontsize', 8); 
% xlabel('Axial Strain');
% % yyaxis(s1,'right');
% % plot(s1,strain15,e_kur15,'linewidth',2); hold on; 
% % axis tight;
% xlim([0 max(strain15)]);
% % ylabel ('Microcrack Energy Kurtosis', 'fontsize', 8);
% 
% subplot(2,1,2);
% % yyaxis(s3,'left');
% % plot(s3,strain15,moment_diff15,'linewidth',2); hold on; 
% % % axis tight; 
% % ylabel('Moment Range', 'fontsize', 8); 
% % xlabel('Axial Strain');
% % yyaxis(s3,'right');
% plot(strain15,e_kur15,'linewidth',2); hold on; 
% % axis tight;
% xlim([0 max(strain15)]);
% ylabel ('Microcrack Energy Kurtosis', 'fontsize', 8);
% %legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);
% 
% set(gcf,'Color','w'); 
% export_fig f91 f91_15mpa_e_stats -q101 -painters -nocrop -pdf -png -tiff -eps
% 
% f92=figure(92);
% 
% %microcrack energy variance and D-value
% subplot(2,1,1);
% % yyaxis(s2,'left');
% plot(strain15,d_value15,'linewidth',2); hold on; 
% % axis tight; 
% ylabel('Fractal Dimension (D-value)', 'fontsize', 8); 
% xlabel('Axial Strain');
% % yyaxis(s1,'right');
% % plot(s1,strain15,e_kur15,'linewidth',2); hold on; 
% % axis tight;
% xlim([0 max(strain15)]);
% % ylabel ('Microcrack Energy Kurtosis', 'fontsize', 8);
% 
% set(gcf,'Color','w'); 
% export_fig f92 f91_15mpa_dvalue -q101 -painters -nocrop -pdf -png -tiff -eps
% 
% f93=figure(93);
% %microcrack rate and shear fraction
% subplot(2,1,1);
% % yyaxis(s2,'left');
% plot(strain15,mc_rate15,'linewidth',2); hold on; 
% ylabel('Microcracking Rate', 'fontsize', 8); 
% xlabel('Axial Strain');
% % yyaxis(s1,'right');
% % plot(s1,strain15,e_kur15,'linewidth',2); hold on; 
% % axis tight;
% xlim([0 max(strain15)]);
% % ylabel ('Microcrack Energy Kurtosis', 'fontsize', 8);
% 
% subplot(2,1,2);
% % yyaxis(s3,'left');
% % plot(s3,strain15,moment_diff15,'linewidth',2); hold on; 
% % % axis tight; 
% % ylabel('Moment Range', 'fontsize', 8); 
% % xlabel('Axial Strain');
% % yyaxis(s3,'right');
% plot(strain15,shear_fraction15,'linewidth',2); hold on; 
% % axis tight;
% xlim([0 max(strain15)]);
% ylabel ('Microcrack Shear Fraction', 'fontsize', 8);
% %legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);
% 
% set(gcf,'Color','w'); 
% export_fig f91 f91_15mpa_mc_stats -q101 -painters -nocrop -pdf -png -tiff -eps
% 
% f94=figure(94);
% %microcrack rate and shear fraction
% subplot(2,1,1);
% % yyaxis(s2,'left');
% plot(strain15,b_value15,'linewidth',2); hold on; 
% ylabel('b-value', 'fontsize', 8); 
% xlabel('Axial Strain');
% % yyaxis(s1,'right');
% % plot(s1,strain15,e_kur15,'linewidth',2); hold on; 
% % axis tight;
% xlim([0 max(strain15)]);
% % ylabel ('Microcrack Energy Kurtosis', 'fontsize', 8);
% 
% subplot(2,1,2);
% % yyaxis(s3,'left');
% % plot(s3,strain15,moment_diff15,'linewidth',2); hold on; 
% % % axis tight; 
% % ylabel('Moment Range', 'fontsize', 8); 
% % xlabel('Axial Strain');
% % yyaxis(s3,'right');
% plot(strain15,moment_diff15,'linewidth',2); hold on; 
% % axis tight;
% xlim([0 max(strain15)]);
% ylabel ('Moment Range', 'fontsize', 8);
% %legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);
% 
% set(gcf,'Color','w'); 
% export_fig f94 f94_15mpa_moment_stats -q101 -painters -nocrop -pdf -png -tiff -eps
% 


