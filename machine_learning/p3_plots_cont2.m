close all
f_cmap=linspecer(4);

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
stress_tf0(1:ind0)=-1*stress_tf0(1:ind0);
ind5=find(stress_tf5==0);
stress_tf5(1:ind5)=-1*stress_tf5(1:ind5);
ind10=find(stress_tf10==0);
stress_tf10(1:ind10)=-1*stress_tf10(1:ind10);
ind15=find(stress_tf15==0);
stress_tf15(1:ind15)=-1*stress_tf15(1:ind15);
ind20=find(stress_tf20==0);
stress_tf20(1:ind20)=-1*stress_tf20(1:ind20);
ind25=find(stress_tf25==0);
stress_tf25(1:ind25)=-1*stress_tf25(1:ind25);
ind30=find(stress_tf30==0);
stress_tf30(1:ind30)=-1*stress_tf30(1:ind30);
ind35=find(stress_tf35==0);
stress_tf35(1:ind35)=-1*stress_tf35(1:ind35);
ind40=find(stress_tf40==0);
stress_tf40(1:ind40)=-1*stress_tf40(1:ind40);
ind45=find(stress_tf45==0);
stress_tf45(1:ind45)=-1*stress_tf45(1:ind45);
ind50=find(stress_tf50==0);
stress_tf50(1:ind50)=-1*stress_tf50(1:ind50);
%%
%MC rate vs strain
f12=figure(12);
plot(strain0,mc_rate0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain15,mc_rate15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain30,mc_rate30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain50,mc_rate50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([0 max(strain50)]);
ylabel('Microcracking Rate'); 
xlabel('Axial Strain');
legend('0 MPa', '30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f12 f12a_MCrate -q101 -painters -nocrop -pdf -png -tiff -eps 
%MC rate vs strain-to-failure 
f13=figure(13);
plot(strain_tf0,mc_rate0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain_tf15,mc_rate15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain_tf30,mc_rate30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain_tf50,mc_rate50,'linewidth',2,'color',f_cmap(4,:)); hold on;
% xlim([0 max(strain_tf50)]);
xlim([min(strain_tf50) max(strain_tf0)]);
ylabel('Microcracking Rate'); 
xlabel('Axial Strain to failure');
legend('0 MPa','30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f13 f12b_MCrate_straintf -q101 -painters -nocrop -pdf -png -tiff -eps 
%MC rate vs stress-to-failure
f14=figure(14);
plot(stress_tf0,mc_rate0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(stress_tf15,mc_rate15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(stress_tf30,mc_rate30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(stress_tf50,mc_rate50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(stress_tf50) max(stress_tf50)]);
ylabel('Microcracking Rate'); 
xlabel('S1 to failure (MPa)');
legend('0 MPa', '30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f14 f12c_MCrate_stresstf -q101 -painters -nocrop -pdf -png -tiff -eps 

%%
% Shear fraction vs strain
f15=figure(15);
plot(strain0,shear_fraction0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain15,shear_fraction15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain30,shear_fraction30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain50,shear_fraction50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([0 max(strain50)]);
ylabel('Shear Microcrack Fraction'); 
xlabel('Axial Strain');
legend('0 MPa', '30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f15 f13a_SFraction -q101 -painters -nocrop -pdf -png -tiff -eps 
%SF vs strain-to-failure 
f16=figure(16);
plot(strain_tf0,shear_fraction0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain_tf15,shear_fraction15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain_tf30,shear_fraction30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain_tf50,shear_fraction50,'linewidth',2,'color',f_cmap(4,:)); hold on;
% xlim([0 max(strain_tf50)]);
xlim([min(strain_tf50) max(strain_tf0)]);
ylabel('Shear Microcrack Fraction'); 
xlabel('Axial Strain to failure');
legend('0 MPa', '30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f16 f13b_SFraction -q101 -painters -nocrop -pdf -png -tiff -eps 
%SF vs stress-to-failure
f17=figure(17);
plot(stress_tf0,shear_fraction0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(stress_tf15,shear_fraction15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(stress_tf30,shear_fraction30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(stress_tf50,shear_fraction50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(stress_tf50) max(stress_tf50)]);
ylabel('Shear Microcrack Fraction'); 
xlabel('S1 to failure (MPa)');
legend('0 MPa','30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f17 f13c_SFraction -q101 -painters -nocrop -pdf -png -tiff -eps 

%%
% Energy variance vs strain
f18=figure(18);
plot(strain0,e_var0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain15,e_var15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain30,e_var30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain50,e_var50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([0 max(strain50)]);
ylabel('Microcrack Energy Variance'); 
xlabel('Axial Strain');
legend('0 MPa', '30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f18 f14a_Evar -q101 -painters -nocrop -pdf -png -tiff -eps 
%Energy Variance vs strain-to-failure 
f19=figure(19);
plot(strain_tf0,e_var0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain_tf15,e_var15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain_tf30,e_var30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain_tf50,e_var50,'linewidth',2,'color',f_cmap(4,:)); hold on;
% xlim([0 max(strain_tf50)]);
xlim([min(strain_tf50) max(strain_tf0)]);
ylabel('Microcrack Energy Variance'); 
xlabel('Axial Strain to failure');
legend('0 MPa', '30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f19 f14b_Evar -q101 -painters -nocrop -pdf -png -tiff -eps 
%Energy Variance vs stress-to-failure
f20=figure(20);
plot(stress_tf0,e_var0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(stress_tf15,e_var15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(stress_tf30,e_var30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(stress_tf50,e_var50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(stress_tf50) max(stress_tf50)]);
ylabel('Microcrack Energy Variance'); 
xlabel('S1 to failure (MPa)');
legend('0 MPa', '30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f20 f14c_Evar -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
% Energy kurtosis vs strain

e_kur0=data0(:,26);
e_kur15=data15(:,26);
e_kur30=data30(:,26);
e_kur50=data50(:,26);

f32=figure(32);
plot(strain0,e_kur0,'linewidth',2,'color',f_cmap(1,:)); hold on;
% plot(strain15,e_kur15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain30,e_kur30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain50,e_kur50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([0 max(strain50)]);
ylabel('Microcrack Energy Kurtosis'); 
xlabel('Axial Strain');
legend('0 MPa', '30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f32 f32a_Ekur -q101 -painters -nocrop -pdf -png -tiff -eps 
%Energy Variance vs strain-to-failure 
f33=figure(33);
plot(strain_tf0,e_kur0,'linewidth',2,'color',f_cmap(1,:)); hold on;
plot(strain_tf15,e_kur15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain_tf30,e_kur30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain_tf50,e_kur50,'linewidth',2,'color',f_cmap(4,:)); hold on;
% xlim([0 max(strain_tf50)]);
xlim([min(strain_tf50) max(strain_tf0)]);
ylabel('Microcrack Energy Kurtosis'); 
xlabel('Axial Strain to failure');
legend('0 MPa', '15 MPa','30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f33 f33b_Ekur -q101 -painters -nocrop -pdf -png -tiff -eps 
%Energy Variance vs stress-to-failure
f34=figure(34);
plot(stress_tf0,e_kur0,'linewidth',2,'color',f_cmap(1,:)); hold on;
plot(stress_tf15,e_kur15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(stress_tf30,e_kur30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(stress_tf50,e_kur50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(stress_tf50) max(stress_tf50)]);
ylabel('Microcrack Energy Kurtosis'); 
xlabel('S1 to failure (MPa)');
legend('0 MPa', '15 MPa','30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f34 f34c_Evar -q101 -painters -nocrop -pdf -png -tiff -eps 

%%
% Moment Range vs strain
f21=figure(21);
plot(strain0,moment_diff0,'linewidth',2,'color',f_cmap(1,:)); hold on;
plot(strain15,moment_diff15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain30,moment_diff30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain50,moment_diff50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([0 max(strain50)]);
ylabel('Moment Range'); 
xlabel('Axial Strain');
legend('0 MPa', '15 MPa','30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f21 f15a_MRange -q101 -painters -nocrop -pdf -png -tiff -eps 
%Moment Range vs strain-to-failure 
f22=figure(22);
plot(strain_tf0,moment_diff0,'linewidth',2,'color',f_cmap(1,:)); hold on;
plot(strain_tf15,moment_diff15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain_tf30,moment_diff30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain_tf50,moment_diff50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(strain_tf50) max(strain_tf0)]);
ylabel('Moment Range'); 
xlabel('Axial Strain to failure');
legend('0 MPa', '15 MPa','30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f22 f15b_MRange -q101 -painters -nocrop -pdf -png -tiff -eps 
%Energy Variance vs stress-to-failure
f23=figure(23);
plot(stress_tf0,moment_diff0,'linewidth',2,'color',f_cmap(1,:)); hold on;
plot(stress_tf15,moment_diff15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(stress_tf30,moment_diff30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(stress_tf50,moment_diff50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(stress_tf50) max(stress_tf50)]);
ylabel('Moment Range'); 
xlabel('S1 to failure (MPa)');
legend('0 MPa', '15 MPa','30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f23 f15c_MRange -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
% b_value vs strain
f24=figure(24);
plot(strain0,b_value0,'linewidth',2,'color',f_cmap(1,:)); hold on;
plot(strain15,b_value15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain30,b_value30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain50,b_value50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([0 max(strain50)]);
ylabel('b-value'); 
xlabel('Axial Strain');
legend('0 MPa', '15 MPa','30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f24 f16a_bvalue -q101 -painters -nocrop -pdf -png -tiff -eps 
%b_value vs strain-to-failure 
f25=figure(25);
plot(strain_tf0,b_value0,'linewidth',2,'color',f_cmap(1,:)); hold on;
plot(strain_tf15,b_value15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain_tf30,b_value30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain_tf50,b_value50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(strain_tf50) max(strain_tf0)]);
ylabel('b-value'); 
xlabel('Axial Strain to failure');
legend('0 MPa', '15 MPa','30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f25 f16b_bvalue -q101 -painters -nocrop -pdf -png -tiff -eps 
%b_value vs stress-to-failure
f26=figure(26);
plot(stress_tf0,b_value0,'linewidth',2,'color',f_cmap(1,:)); hold on;
plot(stress_tf15,b_value15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(stress_tf30,b_value30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(stress_tf50,b_value50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(stress_tf50) max(stress_tf50)]);
ylabel('b-value'); 
xlabel('S1 to failure (MPa)');
legend('0 MPa', '15 MPa','30 MPa','50 MPa','location','best','fontsize',9);

set(gcf,'Color','w'); 
export_fig f26 f16c_bvalue -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
% D_value vs strain
f27=figure(27);
plot(strain0,d_value0,'linewidth',2,'color',f_cmap(1,:)); hold on;
plot(strain15,d_value15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain30,d_value30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain50,d_value50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([0 max(strain50)]);
ylabel('D-value'); 
xlabel('Axial Strain');
legend('0 MPa', '15 MPa','30 MPa','50 MPa','location','best','fontsize',9);
set(gcf,'Color','w'); 
export_fig f27 f17a_Dvalue -q101 -painters -nocrop -pdf -png -tiff -eps 
%D_value vs strain-to-failure 
f28=figure(28);
plot(strain_tf0,d_value0,'linewidth',2,'color',f_cmap(1,:)); hold on;
plot(strain_tf15,d_value15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(strain_tf30,d_value30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(strain_tf50,d_value50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(strain_tf50) max(strain_tf0)]);
ylabel('b-value'); 
xlabel('Axial Strain to failure');
legend('0 MPa', '15 MPa','30 MPa','50 MPa','location','best','fontsize',9);
set(gcf,'Color','w'); 
export_fig f28 f17b_dvalue -q101 -painters -nocrop -pdf -png -tiff -eps 
%D_value vs stress-to-failure
f29=figure(29);
plot(stress_tf0,d_value0,'linewidth',2,'color',f_cmap(1,:)); hold on;
plot(stress_tf15,d_value15,'linewidth',2,'color',f_cmap(2,:)); hold on;
plot(stress_tf30,d_value30,'linewidth',2,'color',f_cmap(3,:)); hold on;
plot(stress_tf50,d_value50,'linewidth',2,'color',f_cmap(4,:)); hold on;
xlim([min(stress_tf50) max(stress_tf50)]);
ylabel('D-value'); 
xlabel('S1 to failure (MPa)');
legend('0 MPa', '15 MPa','30 MPa','50 MPa','location','best','fontsize',9);
set(gcf,'Color','w'); 
export_fig f29 f17c_Dvalue -q101 -painters -nocrop -pdf -png -tiff -eps 
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
run p3_ML_io
