f91=figure(91);
%microcrack rate and shear fraction
e_kur10=data10(:,26);

%microcrack energy variance and D-value
subplot(2,1,1);
% yyaxis(s2,'left');
plot(strain10,e_var10,'linewidth',2); hold on; 
% axis tight; 
ylabel('Microcrack Energy Variance', 'fontsize', 8); 
xlabel('Axial Strain');
% yyaxis(s1,'right');
% plot(s1,strain10,e_kur10,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain10)]);
% ylabel ('Microcrack Energy Kurtosis', 'fontsize', 8);

subplot(2,1,2);
% yyaxis(s3,'left');
% plot(s3,strain10,moment_diff10,'linewidth',2); hold on; 
% % axis tight; 
% ylabel('Moment Range', 'fontsize', 8); 
% xlabel('Axial Strain');
% yyaxis(s3,'right');
plot(strain10e_kur10,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain10)]);
ylabel ('Microcrack Energy Kurtosis', 'fontsize', 8);
%legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f91 f91_10mpa_e_stats -q101 -painters -nocrop -pdf -png -tiff -eps
%%
f92=figure(92);

%microcrack energy variance and D-value
subplot(2,1,1);
% yyaxis(s2,'left');
plot(strain10,d_value10,'linewidth',2); hold on; 
% axis tight; 
ylabel('Fractal Dimension (D-value)', 'fontsize', 8); 
xlabel('Axial Strain');
% yyaxis(s1,'right');
% plot(s1,strain10,e_kur10,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain10)]);
% ylabel ('Microcrack Energy Kurtosis', 'fontsize', 8);

set(gcf,'Color','w'); 
export_fig f92 f91_10mpa_dvalue -q101 -painters -nocrop -pdf -png -tiff -eps



%%
f93=figure(93);
%microcrack rate and shear fraction
subplot(2,1,1);
% yyaxis(s2,'left');
plot(strain10,mc_rate10,'linewidth',2); hold on; 
ylabel('Microcracking Rate', 'fontsize', 8); 
xlabel('Axial Strain');
% yyaxis(s1,'right');
% plot(s1,strain10,e_kur10,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain10)]);
% ylabel ('Microcrack Energy Kurtosis', 'fontsize', 8);

subplot(2,1,2);
% yyaxis(s3,'left');
% plot(s3,strain10,moment_diff10,'linewidth',2); hold on; 
% % axis tight; 
% ylabel('Moment Range', 'fontsize', 8); 
% xlabel('Axial Strain');
% yyaxis(s3,'right');
plot(strain10,shear_fraction10,'linewidth',2); hold on; 
% axis tight;
xlim([0 max(strain10)]);
ylabel ('Microcrack Shear Fraction', 'fontsize', 8);
%legend(s3,'Moment Range', 'b-value','Location','SouthEast','fontsize',6);

set(gcf,'Color','w'); 
export_fig f91 f91_10mpa_mc_stats -q101 -painters -nocrop -pdf -png -tiff -eps


