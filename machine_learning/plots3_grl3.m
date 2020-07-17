clc
clf
clear all
close all
%paths
path(path,genpath('C:\export_fig\update2'));
path(path,genpath('C:\Program Files\Glyph & Cog\XpdfReader-win64'));
path(path,genpath('C:\kakearney-legendflex-pkg-f29cb4e'));
path(path,genpath('C:\kakearney-legendflex-pkg-f29cb4e\legendflex'));
path(path,genpath('C:\Program Files\gs\gs9.22\bin'));
path(path,genpath('C:\Program Files\MATLAB\R2010a\toolbox\linspecer'));
path(path,genpath('C:\xpdf\xpdf-tools-win-4.00'));
path(path,genpath('C:\xpdf\xpdf-tools-win-4.00\bin32'));
path(path,genpath('C:\matlab_adds'));
warning('off','all')
%%
%10 mpa
data10=xlsread('p3_data_binned4.xlsx','10mpa');
strain10=data10(:,29);
mc_rate10=data10(:,23);
shear_fraction10=data10(:,12);
e_var10=data10(:,25);
d_value10=data10(:,38);
moment_diff10=data10(:,33);
b_value10=data10(:,37);
blocs10=find(b_value10(:)>1);
b_value10(blocs10)=b_value10(blocs10)-0.3;

f3=figure(3);
%microcrack rate and shear fraction
s1=subplot(3,1,1);
yyaxis(s1,'left');
plot(s1,strain10,mc_rate10,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain10,shear_fraction10,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain10)]);
ylabel ('Shear Microcrack Fraction', 'fontsize', 8);
% legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%microcrack energy variance and D-value
s2=subplot(3,1,2);
yyaxis(s2,'left');
plot(s2,strain10,e_var10,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcrack Energy Variance', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain10,d_value10,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain10)]);
ylabel ('D-value', 'fontsize', 8);
% legend(s2,'Microcrack Energy Variance', 'D-value','fontsize',6,'Location','East');

%moment difference and b-value
s3=subplot(3,1,3);
yyaxis(s3,'left');
plot(s3,strain10,moment_diff10,'linewidth',2); hold on; 
% axis tight; 
ylabel('Moment Range', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s3,'right');
plot(s3,strain10,b_value10,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain10)]);
ylabel ('b-value', 'fontsize', 8);
% legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f3 f3_10mpa_stats -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%5 mpa
data5=xlsread('p3_data_binned4.xlsx','5mpa');
strain5=data5(:,29);
mc_rate5=data5(:,23);
shear_fraction5=data5(:,12);
e_var5=data5(:,25);
d_value5=data5(:,38);
moment_diff5=data5(:,33);
b_value5=data5(:,37);
blocs5=find(b_value5(:)>1);
b_value5(blocs5)=b_value5(blocs5)-0.3;

f2=figure(2);
%microcrack rate and shear fraction
s1=subplot(3,1,1);
yyaxis(s1,'left');
plot(s1,strain5,mc_rate5,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain5,shear_fraction5,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain5)]);
ylabel ('Shear Microcrack Fraction', 'fontsize', 8);
%legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%microcrack energy variance and D-value
s2=subplot(3,1,2);
yyaxis(s2,'left');
plot(s2,strain5,e_var5,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcrack Energy Variance', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain5,d_value5,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain5)]);
ylabel ('D-value', 'fontsize', 8);
%legend(s2,'Microcrack Energy Variance', 'D-value','fontsize',6,'Location','East');

%moment difference and b-value
s3=subplot(3,1,3);
yyaxis(s3,'left');
plot(s3,strain5,moment_diff5,'linewidth',2); hold on; 
% axis tight; 
ylabel('Moment Range', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s3,'right');
plot(s3,strain5,b_value5,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain5)]);
ylabel ('b-value', 'fontsize', 8);
%legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f2 f2_5mpa_stats -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%0 mpa
data0=xlsread('p3_data_binned4.xlsx','0mpa');
strain0=data0(:,29);
mc_rate0=data0(:,23);
shear_fraction0=data0(:,12);
e_var0=data0(:,25);
d_value0=data0(:,38);
moment_diff0=data0(:,33);
b_value0=data0(:,37);
blocs0=find(b_value0(:)>1);
b_value0(blocs0)=b_value0(blocs0)-0.3;

f1=figure(1);
%microcrack rate and shear fraction
s1=subplot(3,1,1);
yyaxis(s1,'left');
plot(s1,strain0,mc_rate0,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain0,shear_fraction0,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain0)]);
ylabel ('Shear Microcrack Fraction', 'fontsize', 8);
%legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%microcrack energy variance and D-value
s2=subplot(3,1,2);
yyaxis(s2,'left');
plot(s2,strain0,e_var0,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcrack Energy Variance', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain0,d_value0,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain0)]);
ylabel ('D-value', 'fontsize', 8);
%legend(s2,'Microcrack Energy Variance', 'D-value','fontsize',6,'Location','East');

%moment difference and b-value
s3=subplot(3,1,3);
yyaxis(s3,'left');
plot(s3,strain0,moment_diff0,'linewidth',2); hold on; 
% axis tight; 
ylabel('Moment Range', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s3,'right');
plot(s3,strain0,b_value0,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain0)]);
ylabel ('b-value', 'fontsize', 8);
%legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f1 f1_0mpa_stats -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%15 mpa
data15=xlsread('p3_data_binned4.xlsx','15mpa');
strain15=data15(:,29);
mc_rate15=data15(:,23);
shear_fraction15=data15(:,12);
e_var15=data15(:,25);
d_value15=data15(:,38);
moment_diff15=data15(:,33);
b_value15=data15(:,37);
blocs15=find(b_value15(:)>1);
b_value15(blocs15)=b_value15(blocs15)-0.3;

f4=figure(4);
%microcrack rate and shear fraction
s1=subplot(3,1,1);
yyaxis(s1,'left');
plot(s1,strain15,mc_rate15,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain15,shear_fraction15,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain15)]);
ylabel ('Shear Microcrack Fraction', 'fontsize', 8);
%legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%microcrack energy variance and D-value
s2=subplot(3,1,2);
yyaxis(s2,'left');
plot(s2,strain15,e_var15,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcrack Energy Variance', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain15,d_value15,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain15)]);
ylabel ('D-value', 'fontsize', 8);
%legend(s2,'Microcrack Energy Variance', 'D-value','fontsize',6,'Location','East');

%moment difference and b-value
s3=subplot(3,1,3);
yyaxis(s3,'left');
plot(s3,strain15,moment_diff15,'linewidth',2); hold on; 
% axis tight; 
ylabel('Moment Range', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s3,'right');
plot(s3,strain15,b_value15,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain15)]);
ylabel ('b-value', 'fontsize', 8);
%legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f4 f4_15mpa_stats -q101 -painters -nocrop -pdf -png -tiff -eps
% keyboard
%%
%20 mpa
data20=xlsread('p3_data_binned4.xlsx','20mpa');
strain20=data20(:,29);
mc_rate20=data20(:,23);
shear_fraction20=data20(:,12);
e_var20=data20(:,25);
d_value20=data20(:,38);
moment_diff20=data20(:,33);
b_value20=data20(:,37);
blocs20=find(b_value20(:)>1);
b_value20(blocs20)=b_value20(blocs20)-0.3;

f5=figure(5);
%microcrack rate and shear fraction
s1=subplot(3,1,1);
yyaxis(s1,'left');
plot(s1,strain20,mc_rate20,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain20,shear_fraction20,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain20)]);
ylabel ('Shear Microcrack Fraction', 'fontsize', 8);
%legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%microcrack energy variance and D-value
s2=subplot(3,1,2);
yyaxis(s2,'left');
plot(s2,strain20,e_var20,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcrack Energy Variance', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain20,d_value20,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain20)]);
ylabel ('D-value', 'fontsize', 8);
%legend(s2,'Microcrack Energy Variance', 'D-value','fontsize',6,'Location','East');

%moment difference and b-value
s3=subplot(3,1,3);
yyaxis(s3,'left');
plot(s3,strain20,moment_diff20,'linewidth',2); hold on; 
% axis tight; 
ylabel('Moment Range', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s3,'right');
plot(s3,strain20,b_value20,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain20)]);
ylabel ('b-value', 'fontsize', 8);
%legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f5 f5_20mpa_stats -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%25 mpa
data25=xlsread('p3_data_binned4.xlsx','25mpa');
strain25=data25(:,29);
mc_rate25=data25(:,23);
shear_fraction25=data25(:,12);
e_var25=data25(:,25);
d_value25=data25(:,38);
moment_diff25=data25(:,33);
b_value25=data25(:,37)-0.5;
blocs25=find(b_value25(:)>1);
b_value25(blocs25)=b_value25(blocs25)-0.3;

f6=figure(6);
%microcrack rate and shear fraction
s1=subplot(3,1,1);
yyaxis(s1,'left');
plot(s1,strain25,mc_rate25,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain25,shear_fraction25,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain25)]);
ylabel ('Shear Microcrack Fraction', 'fontsize', 8);
%legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%microcrack energy variance and D-value
s2=subplot(3,1,2);
yyaxis(s2,'left');
plot(s2,strain25,e_var25,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcrack Energy Variance', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain25,d_value25,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain25)]);
ylabel ('D-value', 'fontsize', 8);
%legend(s2,'Microcrack Energy Variance', 'D-value','fontsize',6,'Location','East');

%moment difference and b-value
s3=subplot(3,1,3);
yyaxis(s3,'left');
plot(s3,strain25,moment_diff25,'linewidth',2); hold on; 
% axis tight; 
ylabel('Moment Range', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s3,'right');
plot(s3,strain25,b_value25,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain25)]);
ylabel ('b-value', 'fontsize', 8);
%legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f6 f6_25mpa_stats -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%30 mpa
data30=xlsread('p3_data_binned4.xlsx','30mpa');
strain30=data30(:,29);
mc_rate30=data30(:,23);
shear_fraction30=data30(:,12);
e_var30=data30(:,25);
d_value30=data30(:,38);
moment_diff30=data30(:,33);
b_value30=data30(:,37);
blocs30=find(b_value30(:)>1);
b_value30(blocs30)=b_value30(blocs30)-0.3;

f7=figure(7);
%microcrack rate and shear fraction
s1=subplot(3,1,1);
yyaxis(s1,'left');
plot(s1,strain30,mc_rate30,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain30,shear_fraction30,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain30)]);
ylabel ('Shear Microcrack Fraction', 'fontsize', 8);
%legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%microcrack energy variance and D-value
s2=subplot(3,1,2);
yyaxis(s2,'left');
plot(s2,strain30,e_var30,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcrack Energy Variance', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain30,d_value30,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain30)]);
ylabel ('D-value', 'fontsize', 8);
%legend(s2,'Microcrack Energy Variance', 'D-value','fontsize',6,'Location','East');

%moment difference and b-value
s3=subplot(3,1,3);
yyaxis(s3,'left');
plot(s3,strain30,moment_diff30,'linewidth',2); hold on; 
% axis tight; 
ylabel('Moment Range', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s3,'right');
plot(s3,strain30,b_value30,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain30)]);
ylabel ('b-value', 'fontsize', 8);
%legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f7 f7_30mpa_stats -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%35 mpa
data35=xlsread('p3_data_binned4.xlsx','35mpa');
strain35=data35(:,29);
mc_rate35=data35(:,23);
shear_fraction35=data35(:,12);
e_var35=data35(:,25);
d_value35=data35(:,38);
moment_diff35=data35(:,33);
b_value35=data35(:,37);
blocs35=find(b_value35(:)>1);
b_value35(blocs35)=b_value35(blocs35)-0.3;

f8=figure(8);
%microcrack rate and shear fraction
s1=subplot(3,1,1);
yyaxis(s1,'left');
plot(s1,strain35,mc_rate35,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain35,shear_fraction35,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain35)]);
ylabel ('Shear Microcrack Fraction', 'fontsize', 8);
%legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%microcrack energy variance and D-value
s2=subplot(3,1,2);
yyaxis(s2,'left');
plot(s2,strain35,e_var35,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcrack Energy Variance', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain35,d_value35,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain35)]);
ylabel ('D-value', 'fontsize', 8);
%legend(s2,'Microcrack Energy Variance', 'D-value','fontsize',6,'Location','East');

%moment difference and b-value
s3=subplot(3,1,3);
yyaxis(s3,'left');
plot(s3,strain35,moment_diff35,'linewidth',2); hold on; 
% axis tight; 
ylabel('Moment Range', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s3,'right');
plot(s3,strain35,b_value35,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain35)]);
ylabel ('b-value', 'fontsize', 8);
%legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f8 f8_35mpa_stats -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%40 mpa
data40=xlsread('p3_data_binned4.xlsx','40mpa');
strain40=data40(:,29);
mc_rate40=data40(:,23);
shear_fraction40=data40(:,12);
e_var40=data40(:,25);
d_value40=data40(:,38);
moment_diff40=data40(:,33);
b_value40=data40(:,37);
blocs40=find(b_value40(:)>1);
b_value40(blocs40)=b_value40(blocs40)-0.3;

f9=figure(9);
%microcrack rate and shear fraction
s1=subplot(3,1,1);
yyaxis(s1,'left');
plot(s1,strain40,mc_rate40,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain40,shear_fraction40,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain40)]);
ylabel ('Shear Microcrack Fraction', 'fontsize', 8);
%legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%microcrack energy variance and D-value
s2=subplot(3,1,2);
yyaxis(s2,'left');
plot(s2,strain40,e_var40,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcrack Energy Variance', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain40,d_value40,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain40)]);
ylabel ('D-value', 'fontsize', 8);
%legend(s2,'Microcrack Energy Variance', 'D-value','fontsize',6,'Location','East');

%moment difference and b-value
s3=subplot(3,1,3);
yyaxis(s3,'left');
plot(s3,strain40,moment_diff40,'linewidth',2); hold on; 
% axis tight; 
ylabel('Moment Range', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s3,'right');
plot(s3,strain40,b_value40,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain40)]);
ylabel ('b-value', 'fontsize', 8);
%legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f9 f9_40mpa_stats -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
%45 mpa
data45=xlsread('p3_data_binned4.xlsx','45mpa');
strain45=data45(:,29);
mc_rate45=data45(:,23);
shear_fraction45=data45(:,12);
e_var45=data45(:,25);
d_value45=data45(:,38);
moment_diff45=data45(:,33);
b_value45=data45(:,37);
blocs45=find(b_value45(:)>1);
b_value45(blocs45)=b_value45(blocs45)-0.3;

f10=figure(10);
%microcrack rate and shear fraction
s1=subplot(3,1,1);
yyaxis(s1,'left');
plot(s1,strain45,mc_rate45,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain45,shear_fraction45,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain45)]);
ylabel ('Shear Microcrack Fraction', 'fontsize', 8);
%legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%microcrack energy variance and D-value
s2=subplot(3,1,2);
yyaxis(s2,'left');
plot(s2,strain45,e_var45,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcrack Energy Variance', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain45,d_value45,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain45)]);
ylabel ('D-value', 'fontsize', 8);
%legend(s2,'Microcrack Energy Variance', 'D-value','fontsize',6,'Location','East');

%moment difference and b-value
s3=subplot(3,1,3);
yyaxis(s3,'left');
plot(s3,strain45,moment_diff45,'linewidth',2); hold on; 
% axis tight; 
ylabel('Moment Range', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s3,'right');
plot(s3,strain45,b_value45,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain45)]);
ylabel ('b-value', 'fontsize', 8);
%legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f10 f10_45mpa_stats -q101 -painters -nocrop -pdf -png -tiff -eps 

%%
%50 mpa
data50=xlsread('p3_data_binned4.xlsx','50mpa');
strain50=data50(:,29);
mc_rate50=data50(:,23);
shear_fraction50=data50(:,12);
e_var50=data50(:,25);
d_value50=data50(:,38);
moment_diff50=data50(:,33);
b_value50=data50(:,37);
blocs50=find(b_value50(:)>1);
b_value50(blocs50)=b_value50(blocs50)-0.3;

f11=figure(11);
%microcrack rate and shear fraction
s1=subplot(3,1,1);
yyaxis(s1,'left');
plot(s1,strain50,mc_rate50,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s1,'right');
plot(s1,strain50,shear_fraction50,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain50)]);
ylabel ('Shear Microcrack Fraction', 'fontsize', 8);
%legend(s1,'Microcrack Rate', 'Microcrack Shear Fraction','Location','NorthEast','fontsize',6);

%microcrack energy variance and D-value
s2=subplot(3,1,2);
yyaxis(s2,'left');
plot(s2,strain50,e_var50,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcrack Energy Variance', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s2,'right');
plot(s2,strain50,d_value50,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain50)]);
ylabel ('D-value', 'fontsize', 8);
%legend(s2,'Microcrack Energy Variance', 'D-value','fontsize',6,'Location','East');

%moment difference and b-value
s3=subplot(3,1,3);
yyaxis(s3,'left');
plot(s3,strain50,moment_diff50,'linewidth',2); hold on; 
% axis tight; 
ylabel('Moment Range', 'fontsize', 8); 
xlabel('Axial Strain');
yyaxis(s3,'right');
plot(s3,strain50,b_value50,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain50)]);
ylabel ('b-value', 'fontsize', 8);
%legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f11 f11_50mpa_stats -q101 -painters -nocrop -pdf -png -tiff -eps 
%%
% run p3_plots2_c10_35m_new
run plots3_c10_35m
% run plots3_all
run plots3_all2
run plots3_ML_IO6
run p3_prediction_RF5
run grl3_f4
% 
% run p3_plots_cont3
% run p3_plots_cont4
% run p3_ML_io2