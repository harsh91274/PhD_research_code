f91=figure(91);
%microcrack rate and shear fraction
e_kur15=data15(:,26);

%microcrack energy variance and D-value
subplot(2,1,1);
% yyaxis(s2,'left');
plot(strain15,e_var15,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcrack Energy Variance', 'fontsize', 8); 
xlabel('Axial Strain');
% yyaxis(s1,'right');
% plot(s1,strain15,e_kur15,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain15)]);
% ylabel ('Microcrack Energy Kurtosis', 'fontsize', 8);

subplot(2,1,2);
% yyaxis(s3,'left');
% plot(s3,strain15,moment_diff15,'linewidth',2); hold on; 
% % axis tight; 
% ylabel('Moment Range', 'fontsize', 8); 
% xlabel('Axial Strain');
% yyaxis(s3,'right');
plot(strain15,e_kur15,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain15)]);
ylabel ('Microcrack Energy Kurtosis', 'fontsize', 8);
%legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f91 f91_15mpa_e_stats -q101 -painters -nocrop -pdf -png -tiff -eps
%%
f92=figure(92);

%microcrack energy variance and D-value
subplot(2,1,1);
% yyaxis(s2,'left');
plot(strain15,d_value15,'linewidth',2); hold on; 
% axis tight; 
ylabel('Fractal Dimension (D-value)', 'fontsize', 8); 
xlabel('Axial Strain');
% yyaxis(s1,'right');
% plot(s1,strain15,e_kur15,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain15)]);
% ylabel ('Microcrack Energy Kurtosis', 'fontsize', 8);

set(gcf,'Color','w'); 
export_fig f92 f91_15mpa_dvalue -q101 -painters -nocrop -pdf -png -tiff -eps
%%
f93=figure(93);
%microcrack rate and shear fraction
subplot(2,1,1);
% yyaxis(s2,'left');
plot(strain15,mc_rate15,'linewidth',2); hold on; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
% yyaxis(s1,'right');
% plot(s1,strain15,e_kur15,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain15)]);
% ylabel ('Microcrack Energy Kurtosis', 'fontsize', 8);

subplot(2,1,2);
% yyaxis(s3,'left');
% plot(s3,strain15,moment_diff15,'linewidth',2); hold on; 
% % axis tight; 
% ylabel('Moment Range', 'fontsize', 8); 
% xlabel('Axial Strain');
% yyaxis(s3,'right');
plot(strain15,shear_fraction15,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain15)]);
ylabel ('Microcrack Shear Fraction', 'fontsize', 8);
%legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f91 f91_15mpa_mc_stats -q101 -painters -nocrop -pdf -png -tiff -eps
%%
f94=figure(94);
%microcrack rate and shear fraction
subplot(2,1,1);
% yyaxis(s2,'left');
plot(strain15,b_value15,'linewidth',2); hold on; 
ylabel('b-value', 'fontsize', 8); 
xlabel('Axial Strain');
% yyaxis(s1,'right');
% plot(s1,strain15,e_kur15,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain15)]);
% ylabel ('Microcrack Energy Kurtosis', 'fontsize', 8);

subplot(2,1,2);
% yyaxis(s3,'left');
% plot(s3,strain15,moment_diff15,'linewidth',2); hold on; 
% % axis tight; 
% ylabel('Moment Range', 'fontsize', 8); 
% xlabel('Axial Strain');
% yyaxis(s3,'right');
plot(strain15,moment_diff15,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain15)]);
ylabel ('Moment Range', 'fontsize', 8);
%legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f94 f94_15mpa_moment_stats -q101 -painters -nocrop -pdf -png -tiff -eps



