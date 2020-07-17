close all
f_cmap=linspecer(4);

time_tf0=data0(:,46);
time_tf5=data5(:,46);
time_tf10=data10(:,46);
time_tf15=data15(:,46);
time_tf20=data20(:,46);
time_tf25=data25(:,46);
time_tf30=data30(:,46);
time_tf35=data35(:,46);
time_tf40=data40(:,46);
time_tf45=data45(:,46);
time_tf50=data50(:,46);


stress_tf0=data0(:,42);
strain_tf0=data0(:,43);
stress_tf5=data5(:,42);
strain_tf5=data5(:,43);
stress_tf10=data10(:,42);
strain_tf10=data10(:,43);
stress_tf15=data15(:,42);
strain_tf15=data15(:,43);
stress_tf20=data20(:,42);
strain_tf20=data20(:,43);
stress_tf25=data25(:,42);
strain_tf25=data25(:,43);
stress_tf30=data30(:,42);
strain_tf30=data30(:,43);
stress_tf35=data35(:,42);
strain_tf35=data35(:,43);
stress_tf40=data40(:,42);
strain_tf40=data40(:,43);
stress_tf45=data45(:,42);
strain_tf45=data45(:,43);
stress_tf50=data50(:,42);
strain_tf50=data50(:,43);
%modify stress to failure
%negative prior to failure
%positive after failure
ind0=find(stress_tf0==0);
stress_tf0(ind0+1:end)=-1*stress_tf0(ind0+1:end);
ind5=find(stress_tf5==0);
stress_tf5(ind5+1:end)=-1*stress_tf5(ind5+1:end);
ind10=find(stress_tf10==0);
stress_tf10(ind10+1:end)=-1*stress_tf10(ind10+1:end);
ind15=find(stress_tf15==0);
stress_tf15(ind15+1:end)=-1*stress_tf15(ind15+1:end);
ind20=find(stress_tf20==0);
stress_tf20(ind20+1:end)=-1*stress_tf20(ind20+1:end);
ind25=find(stress_tf25==0);
stress_tf25(ind25+1:end)=-1*stress_tf25(ind25+1:end);
ind30=find(stress_tf30==0);
stress_tf30(ind30+1:end)=-1*stress_tf30(ind30+1:end);
ind35=find(stress_tf35==0);
stress_tf35(ind35+1:end)=-1*stress_tf35(ind35+1:end);
ind40=find(stress_tf40==0);
stress_tf40(ind40+1:end)=-1*stress_tf40(ind40+1:end);
ind45=find(stress_tf45==0);
stress_tf45(ind45+1:end)=-1*stress_tf45(ind45+1:end);
ind50=find(stress_tf50==0);
stress_tf50(ind50+1:end)=-1*stress_tf50(ind50+1:end);
%%
mc_var0=data0(:,21);
mc_var5=data5(:,21);
mc_var15=data15(:,21);
mc_var20=data20(:,21);
mc_var25=data25(:,21);
mc_var30=data30(:,21);
mc_var40=data40(:,21);
mc_var45=data45(:,21);
mc_var50=data50(:,21);
%%
%MC rate vs strain
f12=figure(12);
plot(strain0,mc_rate0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain15,mc_rate15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain20,mc_rate20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain45,mc_rate45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([0 max(strain45)]);
ylabel('Microcracking Rate'); 
xlabel('Axial Strain');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f12 f12a_MCrate2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%MC rate vs strain-to-failure 
f13=figure(13);
plot(strain_tf0,mc_rate0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain_tf15,mc_rate15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain_tf20,mc_rate20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain_tf45,mc_rate45,'linewidth',2,'color',f_cmap(4,:)); hold on;
% xlim([0 max(strain_tf45)]);
xlim([min(strain_tf45) max(strain_tf0)]);
ylabel('Microcracking Rate'); 
xlabel('Axial Strain to failure');
legend('0 MPa','20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f13 f12b_MCrate2_straintf -q101 -painters -nocrop -pdf -png -tiff -eps 

%MC rate vs stress-to-failure
f14=figure(14);
plot(stress_tf0,mc_rate0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(stress_tf15,mc_rate15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(stress_tf20,mc_rate20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(stress_tf45,mc_rate45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(stress_tf45) max(stress_tf45)]);
ylabel('Microcracking Rate'); 
xlabel('S1 to failure (MPa)');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f14 f12c_MCrate2_stresstf -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%MC var vs strain
f32=figure(32);
plot(strain0,mc_var0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain15,mc_rate15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain20,mc_var20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain45,mc_var45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([0 max(strain45)]);
ylabel('Microcracking Variance'); 
xlabel('Axial Strain');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);
set(gcf,'Color','w'); 
export_fig f32 f32a_MCvar2 -q101 -painters -nocrop -pdf -png -tiff -eps 

%MC var vs strain-to-failure 
f33=figure(33);
plot(strain_tf0,mc_var0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain_tf15,mc_rate15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain_tf20,mc_var20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain_tf45,mc_var45,'linewidth',2,'color',f_cmap(4,:)); hold on;
% xlim([0 max(strain_tf45)]);
xlim([min(strain_tf45) max(strain_tf0)]);
ylabel('Microcracking Variance'); 
xlabel('Axial Strain to failure');
legend('0 MPa','20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f33 f32b_MCvar2_straintf -q101 -painters -nocrop -pdf -png -tiff -eps 

%MC rate vs stress-to-failure
f34=figure(34);
plot(stress_tf0,mc_var0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(stress_tf15,mc_rate15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(stress_tf20,mc_var20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(stress_tf45,mc_var45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(stress_tf45) max(stress_tf45)]);
ylabel('Microcracking Variance'); 
xlabel('S1 to failure (MPa)');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);
set(gcf,'Color','w'); 
export_fig f34 f32c_MCvar2_stresstf -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
% Shear fraction vs strain
f15=figure(15);
plot(strain0,shear_fraction0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain15,shear_fraction15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain20,shear_fraction20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain45,shear_fraction45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([0 max(strain45)]); ylim([0 1.2]);
ylabel('Shear Microcrack Fraction'); 
xlabel('Axial Strain');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f15 f13a_SFraction2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%SF vs strain-to-failure 
f16=figure(16);
plot(strain_tf0,shear_fraction0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain_tf15,shear_fraction15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain_tf20,shear_fraction20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain_tf45,shear_fraction45,'linewidth',2,'color',f_cmap(4,:)); hold on;
% xlim([0 max(strain_tf45)]);
xlim([min(strain_tf45) max(strain_tf0)]); ylim([0 1.2]);
ylabel('Shear Microcrack Fraction'); 
xlabel('Axial Strain to failure');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f16 f13b_SFraction2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%SF vs stress-to-failure
f17=figure(17);
plot(stress_tf0,shear_fraction0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(stress_tf15,shear_fraction15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(stress_tf20,shear_fraction20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(stress_tf45,shear_fraction45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(stress_tf45) max(stress_tf45)]); ylim([0 1.2]);
ylabel('Shear Microcrack Fraction'); 
xlabel('S1 to failure (MPa)');
legend('0 MPa','20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f17 f13c_SFraction2 -q101 -painters -nocrop -pdf -png -tiff -eps 

%%
% Energy variance vs strain
f18=figure(18);
plot(strain0,e_var0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain15,e_var15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain20,e_var20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain45,e_var45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([0 max(strain45)]);
ylabel('Microcrack Energy Variance'); 
xlabel('Axial Strain');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f18 f14a_Evar2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%Energy Variance vs strain-to-failure 
f19=figure(19);
plot(strain_tf0,e_var0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain_tf15,e_var15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain_tf20,e_var20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain_tf45,e_var45,'linewidth',2,'color',f_cmap(4,:)); hold on;
% xlim([0 max(strain_tf45)]);
xlim([min(strain_tf45) max(strain_tf0)]);
ylabel('Microcrack Energy Variance'); 
xlabel('Axial Strain to failure');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f19 f14b_Evar2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%Energy Variance vs stress-to-failure
f20=figure(20);
plot(stress_tf0,e_var0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(stress_tf15,e_var15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(stress_tf20,e_var20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(stress_tf45,e_var45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(stress_tf45) max(stress_tf45)]);
ylabel('Microcrack Energy Variance'); 
xlabel('S1 to failure (MPa)');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f20 f14c_Evar2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
% Energy kurtosis vs strain

e_kur0=data0(:,26);
e_kur15=data15(:,26);
e_kur20=data20(:,26);
e_kur45=data45(:,26);
e_kur5=data5(:,26);
e_kur30=data30(:,26);
e_kur40=data40(:,26);
e_kur50=data50(:,26);
e_kur25=data25(:,26);

e_rate0=data0(:,27);
e_rate15=data15(:,27);
e_rate20=data20(:,27);
e_rate45=data45(:,27);
e_rate5=data5(:,27);
e_rate30=data30(:,27);
e_rate40=data40(:,27);
e_rate50=data50(:,27);
e_rate25=data25(:,27);

f32=figure(32);
plot(strain0,e_rate0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain15,e_kur15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain20,e_rate20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain45,e_rate45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([0 max(strain45)]);
ylabel('Microcrack Energy Rate'); 
xlabel('Axial Strain');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f32 f32a_Erate2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%Energy Variance vs strain-to-failure 
f33=figure(33);
plot(strain_tf0,e_rate0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain_tf15,e_kur15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain_tf20,e_rate20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain_tf45,e_rate45,'linewidth',2,'color',f_cmap(4,:)); hold on;
% xlim([0 max(strain_tf45)]);
xlim([min(strain_tf45) max(strain_tf0)]);
ylabel('Microcrack Energy Rate'); 
xlabel('Axial Strain to failure');
legend('0 MPa', '20 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f33 f33b_Erate2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%Energy Variance vs stress-to-failure
f34=figure(34);
plot(stress_tf0,e_rate0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(stress_tf15,e_kur15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(stress_tf20,e_rate20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(stress_tf45,e_rate50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(stress_tf45) max(stress_tf45)]);
ylabel('Microcrack Energy Rate'); 
xlabel('S1 to failure (MPa)');
legend('0 MPa', '20 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f34 f34c_Erate2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
% Moment Range vs strain
f21=figure(21);
plot(strain0,moment_diff0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain15,moment_diff15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain20,moment_diff20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain45,moment_diff45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([0 max(strain45)]);
ylabel('Moment Range'); 
xlabel('Axial Strain');
legend('0 MPa','20 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f21 f15a_MRange2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%Moment Range vs strain-to-failure 
f22=figure(22);
plot(strain_tf0,moment_diff0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain_tf15,moment_diff15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain_tf20,moment_diff20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain_tf45,moment_diff45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(strain_tf45) max(strain_tf0)]);
ylabel('Moment Range'); 
xlabel('Axial Strain to failure');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f22 f15b_MRange2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%Energy Variance vs stress-to-failure
f23=figure(23);
plot(stress_tf0,moment_diff0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(stress_tf15,moment_diff15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(stress_tf20,moment_diff20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(stress_tf45,moment_diff45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(stress_tf45) max(stress_tf45)]);
ylabel('Moment Range'); 
xlabel('S1 to failure (MPa)');
legend('0 MPa','20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f23 f15c_MRange2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
% b_value vs strain
f24=figure(24);
plot(strain0,b_value0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain15,b_value15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain20,b_value20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain45,b_value45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([0 max(strain45)]);
ylabel('b-value'); 
xlabel('Axial Strain');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f24 f16a_bvalue2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%b_value vs strain-to-failure 
f25=figure(25);
plot(strain_tf0,b_value0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain_tf15,b_value15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain_tf20,b_value20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain_tf45,b_value45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(strain_tf45) max(strain_tf0)]);
ylabel('b-value'); 
xlabel('Axial Strain to failure');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f25 f16b_bvalue2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%b_value vs stress-to-failure
f26=figure(26);
plot(stress_tf0,b_value0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(stress_tf15,b_value15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(stress_tf20,b_value20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(stress_tf45,b_value45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(stress_tf45) max(stress_tf45)]);
ylabel('b-value'); 
xlabel('S1 to failure (MPa)');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f26 f16c_bvalue2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
% D_value vs strain
f27=figure(27);
plot(strain0,d_value0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain15,d_value15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain20,d_value20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain45,d_value45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([0 max(strain45)]);
ylabel('D-value'); 
xlabel('Axial Strain');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);
set(gcf,'Color','w'); 
export_fig f27 f17a_Dvalue2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%D_value vs strain-to-failure 
f28=figure(28);
plot(strain_tf0,d_value0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain_tf15,d_value15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain_tf20,d_value20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain_tf45,d_value45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(strain_tf45) max(strain_tf0)]);
ylabel('D-value'); 
xlabel('Axial Strain to failure');
legend('0 MPa','20 MPa','45 MPa','location','best','fontsize',9);
set(gcf,'Color','w'); 
export_fig f28 f17b_dvalue2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%D_value vs stress-to-failure
f29=figure(29);
plot(stress_tf0,d_value0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(stress_tf15,d_value15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(stress_tf20,d_value20,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(stress_tf45,d_value45,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(stress_tf45) max(stress_tf45)]);
ylabel('D-value'); 
xlabel('S1 to failure (MPa)');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9);
set(gcf,'Color','w'); 
export_fig f29 f17c_Dvalue2 -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%b-value_all
f_cmap2=linspecer(11);
f30=figure(30);
plot(strain0,b_value0,'linewidth',2,'color',f_cmap2(1,:)); hold on;
plot(strain5,b_value5,'linewidth',2,'color',f_cmap2(2,:)); hold on;
plot(strain10,b_value10,'linewidth',2,'color',f_cmap2(3,:)); hold on;
plot(strain15,b_value15,'linewidth',2,'color',f_cmap2(4,:)); hold on;
plot(strain20,b_value20,'linewidth',2,'color',f_cmap2(5,:)); hold on;
plot(strain25,b_value25,'linewidth',2,'color',f_cmap2(6,:)); hold on;
plot(strain30,b_value30,'linewidth',2,'color',f_cmap2(7,:)); hold on;
plot(strain35,b_value35,'linewidth',2,'color',f_cmap2(8,:)); hold on;
plot(strain40,b_value40,'linewidth',2,'color',f_cmap2(9,:)); hold on;
plot(strain45,b_value45,'linewidth',2,'color',f_cmap2(10,:)); hold on;
plot(strain50,b_value50,'linewidth',2,'color',f_cmap2(11,:)); hold on;
xlim([0 max(strain50)]);
ylabel('b-value'); 
xlabel('Axial Strain');
legend('0 MPa', '5 MPa','10 MPa','15 MPa','20 MPa','25 MPa','30 MPa','35 MPa','40 MPa','45 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f30 f30_bvalue_all -q101 -painters -nocrop -pdf -png -tiff -eps 

f31=figure(31);
plot(strain0,d_value0,'linewidth',2,'color',f_cmap2(1,:)); hold on;
plot(strain5,d_value5,'linewidth',2,'color',f_cmap2(2,:)); hold on;
plot(strain10,d_value10,'linewidth',2,'color',f_cmap2(3,:)); hold on;
plot(strain15,d_value15,'linewidth',2,'color',f_cmap2(4,:)); hold on;
plot(strain20,d_value20,'linewidth',2,'color',f_cmap2(5,:)); hold on;
plot(strain25,d_value25,'linewidth',2,'color',f_cmap2(6,:)); hold on;
plot(strain30,d_value30,'linewidth',2,'color',f_cmap2(7,:)); hold on;
plot(strain35,d_value35,'linewidth',2,'color',f_cmap2(8,:)); hold on;
plot(strain40,d_value40,'linewidth',2,'color',f_cmap2(9,:)); hold on;
plot(strain45,d_value45,'linewidth',2,'color',f_cmap2(10,:)); hold on;
plot(strain50,d_value50,'linewidth',2,'color',f_cmap2(11,:)); hold on;
xlim([0 max(strain50)]);
ylabel('D-value'); 
xlabel('Axial Strain');
legend('0 MPa', '5 MPa','10 MPa','15 MPa','20 MPa','25 MPa','30 MPa','35 MPa','40 MPa','45 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f31 f31_Dvalue_all -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%b value vs shear fraction
b_comp=[b_value0;b_value20; b_value45];
sf_comp=[shear_fraction0;shear_fraction20; shear_fraction45];
[fito1, gof1]=fit(sf_comp, b_comp,'poly1');
eq=[fito1.p1 fito1.p2]; 
sf_test=0:0.1:1; 
line_pts=polyval(eq, sf_test);

f39=figure(39);
scatter(shear_fraction0, b_value0, [],f_cmap(1,:)); hold on;
% scatter(shear_fraction5, b_value5, [],f_cmap2(2,:)); hold on;
% scatter(shear_fraction10, b_value10, [],f_cmap2(3,:)); hold on;
% scatter(shear_fraction15, b_value15, [],f_cmap2(4,:)); hold on;
scatter(shear_fraction20, b_value20, [],f_cmap(2,:)); hold on;
% scatter(shear_fraction25, b_value25, [],f_cmap2(6,:)); hold on;
% scatter(shear_fraction30, b_value30, [],f_cmap2(7,:)); hold on;
% scatter(shear_fraction35, b_value35, [],f_cmap2(8,:)); hold on;
% scatter(shear_fraction40, b_value40, [],f_cmap2(9,:)); hold on;
scatter(shear_fraction45, b_value45, [],f_cmap(3,:)); hold on;


plot(sf_test,line_pts,'k--','linewidth',2); hold on;
axis tight
box on
ylabel('b-value'); 
xlabel('Shear Fraction (SF)');
legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9); hold on;
set(gcf,'Color','w'); 
export_fig f39 f39_sf_vs_b -q101 -painters -nocrop -pdf -png -tiff -eps 
%%

b_comp=[b_value0;b_value20; b_value45];
d_comp=[d_value0;d_value20; d_value45];
[fito2, gof2]=fit(d_comp, b_comp,'poly1');
eq=[fito2.p1 fito2.p2]; 
% d_test=0:0.1:1; 
% line_pts=polyval(eq, sf_test);

f40=figure(40);
scatter(d_value0, b_value0, [],f_cmap(1,:)); hold on;
% scatter(d_value5, b_value5, [],f_cmap2(2,:)); hold on;
% scatter(d_value10, b_value10, [],f_cmap2(3,:)); hold on;
% scatter(d_value15, b_value15, [],f_cmap2(4,:)); hold on;
scatter(d_value20, b_value20, [],f_cmap(2,:)); hold on;
% scatter(d_value25, b_value25, [],f_cmap2(6,:)); hold on;
% scatter(d_value30, b_value30, [],f_cmap2(7,:)); hold on;
% scatter(d_value35, b_value35, [],f_cmap2(8,:)); hold on;
% scatter(d_value40, b_value40, [],f_cmap2(9,:)); hold on;
scatter(d_value45, b_value45, [],f_cmap(3,:)); hold on;

% plot(sf_test,line_pts,'k--','linewidth',2); hold on;
axis tight
box on
ylabel('b-value'); 
xlabel('D-value');
% legend('0 MPa', '20 MPa','45 MPa','location','best','fontsize',9); hold on;
set(gcf,'Color','w'); 
% export_fig f39 f39_sf_vs_b -q101 -painters -nocrop -pdf -png -tiff -eps 

%%
stress=xlsread('p3_data_binned4.xlsx','initial_S1');
strain_all=stress(:,1);
time_all=stress(:,2);
b_proxy=nan(10,1);d_proxy=nan(10,1);sf_proxy=nan(10,1);mcvar_proxy=nan(10,1);

b_value0_t=[b_proxy; b_value0];
d_value0_t=[d_proxy; d_value0];
shear_fraction0_t=[sf_proxy; shear_fraction0];
mc_var0_t=[mcvar_proxy; mc_var0];

b_value10_t=[b_proxy; b_value10];
d_value10_t=[d_proxy; d_value10];
shear_fraction10_t=[sf_proxy; shear_fraction10];
mc_var10_t=[mcvar_proxy; mc_var10];

b_value20_t=[b_proxy; b_value20];
d_value20_t=[d_proxy; d_value20];
shear_fraction20_t=[sf_proxy; shear_fraction20];
mc_var20_t=[mcvar_proxy; mc_var20];

b_value30_t=[b_proxy; b_value30];
d_value30_t=[d_proxy; d_value30];
shear_fraction30_t=[sf_proxy; shear_fraction30];
mc_var30_t=[mcvar_proxy; mc_var30];

b_value40_t=[b_proxy; b_value40];
d_value40_t=[d_proxy; d_value40];
shear_fraction40_t=[sf_proxy; shear_fraction40];
mc_var40_t=[mcvar_proxy; mc_var40];

b_value50_t=[b_proxy; b_value50];
d_value50_t=[d_proxy; d_value50];
shear_fraction50_t=[sf_proxy; shear_fraction50];
mc_var50_t=[mcvar_proxy; mc_var50];
%%
f510=figure(510);
%mcvar and sf 
s1=subplot(2,1,1);
yyaxis(s1,'left');
plot(s1,time_all,mc_var0_t,'linewidth',2); hold on; 
ylabel ('MC_{var}', 'fontsize', 8); hold on; 
yyaxis(s1,'right');
plot(s1,time_all,shear_fraction0_t,'linewidth',2); hold on; 
ylabel('Shear Fraction (SF)','fontsize',8);
box on
xlim([0 max(time_all)]);
xlabel('Time (s)'); hold on;

s2=subplot(2,1,2);
yyaxis(s2,'left');
plot(s2,time_all,b_value0_t,'linewidth',2); hold on; 
ylabel ('b-value', 'fontsize', 8);
yyaxis(s2,'right');
plot(s2,time_all,d_value0_t,'linewidth',2); hold on; 
ylabel('Fractal Dimension (D_{2})','fontsize',8);
box on
xlim([0 max(time_all)]);
xlabel('Time (s)'); hold on;

set(gcf,'Color','w'); 
export_fig f510 f510_0mpa_time -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
f511=figure(511);
%mcvar and sf 
s1=subplot(2,1,1);
yyaxis(s1,'left');
plot(s1,time_all,mc_var10_t,'linewidth',2); hold on; 
ylabel ('MC_{var}', 'fontsize', 8); hold on;
yyaxis(s1,'right');
plot(s1,time_all,shear_fraction10_t,'linewidth',2); hold on; 
ylabel('Shear Fraction (SF)','fontsize',8);
box on
xlim([0 max(time_all)]);
xlabel('Time (s)'); hold on;

s2=subplot(2,1,2);
yyaxis(s2,'left');
plot(s2,time_all,b_value10_t,'linewidth',2); hold on; 
ylabel ('b-value', 'fontsize', 8);
yyaxis(s2,'right');
plot(s2,time_all,d_value10_t,'linewidth',2); hold on; 
ylabel('Fractal Dimension (D_{2})','fontsize',8);
box on
xlim([0 max(time_all)]);
xlabel('Time (s)'); hold on;

set(gcf,'Color','w'); 
export_fig f511 f511_10mpa_time -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
f512=figure(512);
%mcvar and sf 
s1=subplot(2,1,1);
yyaxis(s1,'left');
plot(s1,time_all,mc_var20_t,'linewidth',2); hold on; 
ylabel ('MC_{var}', 'fontsize', 8); hold on;
yyaxis(s1,'right');
plot(s1,time_all,shear_fraction20_t,'linewidth',2); hold on; 
ylabel('Shear Fraction (SF)','fontsize',8);
box on
xlim([0 max(time_all)]);
xlabel('Time (s)'); hold on;

s2=subplot(2,1,2);
yyaxis(s2,'left');
plot(s2,time_all,b_value20_t,'linewidth',2); hold on; 
ylabel ('b-value', 'fontsize', 8);
yyaxis(s2,'right');
plot(s2,time_all,d_value20_t,'linewidth',2); hold on; 
ylabel('Fractal Dimension (D_{2})','fontsize',8);
box on
xlim([0 max(time_all)]);
xlabel('Time (s)'); hold on;

set(gcf,'Color','w'); 
export_fig f512 f512_20mpa_time -q101 -painters -nocrop -pdf -png -tiff -eps 

%%
f513=figure(513);
%mcvar and sf 
s1=subplot(2,1,1);
yyaxis(s1,'left');
plot(s1,time_all,mc_var30_t,'linewidth',2); hold on; 
ylabel ('MC_{var}', 'fontsize', 8); hold on;
yyaxis(s1,'right');
plot(s1,time_all,shear_fraction30_t,'linewidth',2); hold on; 
ylabel('Shear Fraction (SF)','fontsize',8);
box on
xlim([0 max(time_all)]);
xlabel('Time (s)'); hold on;

s2=subplot(2,1,2);
yyaxis(s2,'left');
plot(s2,time_all,b_value30_t,'linewidth',2); hold on; 
ylabel ('b-value', 'fontsize', 8);
yyaxis(s2,'right');
plot(s2,time_all,d_value30_t,'linewidth',2); hold on; 
ylabel('Fractal Dimension (D_{2})','fontsize',8);
box on
xlim([0 max(time_all)]);
xlabel('Time (s)'); hold on;

set(gcf,'Color','w'); 
export_fig f513 f513_30mpa_time -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
f514=figure(514);
%mcvar and sf 
s1=subplot(2,1,1);
yyaxis(s1,'left');
plot(s1,time_all,mc_var40_t,'linewidth',2); hold on; 
ylabel ('MC_{var}', 'fontsize', 8); hold on;
yyaxis(s1,'right');
plot(s1,time_all,shear_fraction40_t,'linewidth',2); hold on; 
ylabel('Shear Fraction (SF)','fontsize',8);
box on
xlim([0 max(time_all)]);
xlabel('Time (s)'); hold on;

s2=subplot(2,1,2);
yyaxis(s2,'left');
plot(s2,time_all,b_value40_t,'linewidth',2); hold on; 
ylabel ('b-value', 'fontsize', 8);
yyaxis(s2,'right');
plot(s2,time_all,d_value40_t,'linewidth',2); hold on; 
ylabel('Fractal Dimension (D_{2})','fontsize',8);
box on
xlim([0 max(time_all)]);
xlabel('Time (s)'); hold on;

set(gcf,'Color','w'); 
export_fig f514 f514_40mpa_time -q101 -painters -nocrop -pdf -png -tiff -eps 

%%
f515=figure(515);
%mcvar and sf 
s1=subplot(2,1,1);
yyaxis(s1,'left');
plot(s1,time_all,mc_var50_t,'linewidth',2); hold on; 
ylabel ('MC_{var}', 'fontsize', 8); hold on;
yyaxis(s1,'right');
plot(s1,time_all,shear_fraction50_t,'linewidth',2); hold on; 
ylabel('Shear Fraction (SF)','fontsize',8);
box on
xlim([0 max(time_all)]);
xlabel('Time (s)'); hold on;

s2=subplot(2,1,2);
yyaxis(s2,'left');
plot(s2,time_all,b_value50_t,'linewidth',2); hold on; 
ylabel ('b-value', 'fontsize', 8);
yyaxis(s2,'right');
plot(s2,time_all,d_value50_t,'linewidth',2); hold on; 
ylabel('Fractal Dimension (D_{2})','fontsize',8);
box on
xlim([0 max(time_all)]);
xlabel('Time (s)'); hold on;

set(gcf,'Color','w'); 
export_fig f515 f515_50mpa_time -q101 -painters -nocrop -pdf -png -tiff -eps 