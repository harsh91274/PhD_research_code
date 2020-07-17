idata1=input_master1';
idata2=input_master1';
targets1=o_time_tf';
targets2=o_stress_tf';
strain_sc=0;
stress_sc=0;
time_sc=0;
path(path,genpath('C:\matlab_adds'));
%%
%Neural Network Build
while time_sc~=1
close all
hiddenlayersize1=5;
net1=fitnet(hiddenlayersize1);

%setup division of data
net1.divideParam.trainRatio = 70/100;
net1.divideParam.valRatio = 15/100;
net1.divideParam.testRatio = 15/100;
net1.trainFcn = 'trainlm';
net1.performfcn='sse';
% net1.layers{1}.transferfcn='poslin';

[net1,tr1] = train(net1,idata1,targets1);

%training performance
outputs1 = net1(idata1);
errors1 = gsubtract(outputs1,targets1);
performance1 = perform(net1,targets1,outputs1);

trOut1 = outputs1(tr1.trainInd);
tsOut1 = outputs1(tr1.testInd);
trVal1 = outputs1(tr1.valInd);
tsVal1 = outputs1(tr1.valInd);
trTarg1 = targets1(tr1.trainInd);
tsTarg1 = targets1(tr1.testInd);

% view(net1);
f1_r=figure;
plotregression(trTarg1, trOut1, 'Train', tsTarg1, tsOut1, 'Testing',targets1,outputs1,'All');
set(gcf,'color','w');
export_fig f1_new ANN_time_regression -pdf -png -eps -tiff -dsvg -q101 -nocrop
%%
%locate indexes for plotting
strain_all=[time0; time5; time10; time15; time20; time30; time35; time40; time45; time50];

strain_all_train=strain_all(tr1.trainInd);
strain_all_val=strain_all(tr1.valInd);
% strain_all_train=[strain_all_train; strain_all_val];
strain_all_test=strain_all(tr1.testInd);
strain_all_test=[strain_all_test; strain_all_val];
% strain_all_test=[strain_all_test; strain_all_val];

strain_tf_train=targets1(tr1.trainInd)';
strain_tf_val=targets1(tr1.valInd)';
% strain_tf_train=[strain_tf_train; strain_tf_val];
strain_tf_test=targets1(tr1.testInd)';
strain_tf_test=[strain_tf_test; strain_tf_val];
% strain_tf_test=[strain_tf_test; strain_tf_val];

cp_all=input_master1(:,1);
cp_train=cp_all(tr1.trainInd);
cp_val=cp_all(tr1.valInd);
% cp_train=[cp_train; cp_val];
cp_test=cp_all(tr1.testInd);
cp_test=[cp_test; cp_val];
% cp_test=[cp_test; cp_val];

cp_train0=find(cp_train==0);
cp_test0=find(cp_test==0);
cp_train5=find(cp_train==5);
cp_test5=find(cp_test==5);
cp_train10=find(cp_train==10);
cp_test10=find(cp_test==10);
cp_train15=find(cp_train==15);
cp_test15=find(cp_test==15);
cp_train20=find(cp_train==20);
cp_test20=find(cp_test==20);
cp_train30=find(cp_train==30);
cp_test30=find(cp_test==30);
cp_train35=find(cp_train==35);
cp_test35=find(cp_test==35);
cp_train40=find(cp_train==40);
cp_test40=find(cp_test==40);
cp_train45=find(cp_train==45);
cp_test45=find(cp_test==45);
cp_train50=find(cp_train==50);
cp_test50=find(cp_test==50);

f_cmap11=linspecer(11);
%%
%training plot
outputs1_train=outputs1(tr1.trainInd)';
outputs1_val=outputs1(tr1.valInd)';
% outputs1_train=[outputs1_train; outputs1_val];
outputs1_test=outputs1(tr1.testInd)';
outputs1_test=[outputs1_test; outputs1_val];
%%
f_cmap12=linspecer(3);
f2=figure(2);
subplot(3,3,1);
title('0 MPa'); hold on;
p0=plot(strain_all_train(cp_train0), strain_tf_train (cp_train0), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test(cp_test0), strain_tf_test (cp_test0), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p0=plot(strain_all_train(cp_train0), outputs1_train (cp_train0), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p0=plot(strain_all_test(cp_test0), outputs1_test (cp_test0), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); 
ylabel ('Time to Failure (s)'); 
box on

subplot(3,3,2);
title('5 MPa'); hold on;
p5=plot(strain_all_train(cp_train5), strain_tf_train (cp_train5), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test(cp_test5), strain_tf_test (cp_test5), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p5=plot(strain_all_train(cp_train5), outputs1_train (cp_train5), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p5=plot(strain_all_test(cp_test5), outputs1_test (cp_test5), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); 
% ylabel ('Strain to Failure'); 
box on

subplot(3,3,3);
title('10 MPa'); hold on;
p10=plot(strain_all_train(cp_train10), strain_tf_train (cp_train10), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test(cp_test10), strain_tf_test (cp_test10), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p10=plot(strain_all_train(cp_train10), outputs1_train (cp_train10), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p10=plot(strain_all_test(cp_test10), outputs1_test (cp_test10), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); 
% ylabel ('Strain to Failure'); 
box on

subplot(3,3,4);
title('15 MPa'); hold on;
p15=plot(strain_all_train(cp_train15), strain_tf_train (cp_train15), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test(cp_test15), strain_tf_test (cp_test15), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p15=plot(strain_all_train(cp_train15), outputs1_train (cp_train15), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p15=plot(strain_all_test(cp_test15), outputs1_test (cp_test15), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); 
% ylabel ('Strain to Failure'); 
box on

subplot(3,3,5);
title('20 MPa'); hold on;
p20=plot(strain_all_train(cp_train20), strain_tf_train (cp_train20), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test(cp_test20), strain_tf_test (cp_test20), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p20=plot(strain_all_train(cp_train20), outputs1_train (cp_train20), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p20=plot(strain_all_test(cp_test20), outputs1_test (cp_test20), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); ylabel ('Time to Failure (s)'); 
box on

subplot(3,3,6);
title('30 MPa'); hold on;
p30=plot(strain_all_train(cp_train30), strain_tf_train (cp_train30), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test(cp_test30), strain_tf_test (cp_test30), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p30=plot(strain_all_train(cp_train30), outputs1_train (cp_train30), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p30=plot(strain_all_test(cp_test30), outputs1_test (cp_test30), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); 
% ylabel ('Strain to Failure'); 
box on

subplot(3,3,7);
title('35 MPa'); hold on;
p35=plot(strain_all_train(cp_train35), strain_tf_train (cp_train35), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test(cp_test35), strain_tf_test (cp_test35), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p35=plot(strain_all_train(cp_train35), outputs1_train (cp_train35), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p35=plot(strain_all_test(cp_test35), outputs1_test (cp_test35), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); 
% ylabel ('Strain to Failure'); 
box on

subplot(3,3,8);
title('40 MPa'); hold on;
p40=plot(strain_all_train(cp_train40), strain_tf_train (cp_train40), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test(cp_test40), strain_tf_test (cp_test40), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p40=plot(strain_all_train(cp_train40), outputs1_train (cp_train40), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p40=plot(strain_all_test(cp_test40), outputs1_test (cp_test40), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Axial Strain'); 
% ylabel ('Strain to Failure'); 
box on

subplot(3,3,9);
title('45 MPa'); hold on;
p45=plot(strain_all_train(cp_train45), strain_tf_train (cp_train45), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test(cp_test45), strain_tf_test (cp_test45), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p45=plot(strain_all_train(cp_train45), outputs1_train (cp_train45), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p45=plot(strain_all_test(cp_test45), outputs1_test (cp_test45), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); 
ylabel ('Time to Failure'); 
box on

set(gcf,'color','w');
saveas(gcf,'strain_tf','fig');
export_fig f2 ANN_time_tvt -pdf -png -eps -tiff -dsvg -q101 -nocrop
%%
%25 mpa testing
targets1_25=o_25t_time_tf;
outputs1_25 = sim(net1,input_25t_master1');
errors1_25 = gsubtract(outputs1_25,targets1_25);
performance1_25 = perform(net1,targets1_25,outputs1_25);

f3=figure(3);
title('25 MPa'); hold on;
p25=plot(time25, targets1_25, 'o', 'markersize',4,'color',f_cmap12(1,:)); hold on;
test_p25=plot(time25, outputs1_25, 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',2.5); hold on;
axis tight; xlabel('Axial Strain'); ylabel ('Output: Time to Failure'); 
yyaxis right
test_s25=plot(time25, data25(:,39)./10^6,'linewidth',2); hold on;
ylabel('Axial Stress (MPa)'); hold on;
box on
legend('Testing Set', 'Predictions','location','southeast');
set(gcf,'color','w');
export_fig f3 ANN_Time_test25 -pdf -png -eps -tiff -dsvg -q101 -nocrop

strain_rsq1=rsquare(targets1_25, outputs1_25');
strain_rsq2=rsquare(targets1_25, outputs1_25',false);
disp(strcat('25 MPa Strain-to-Failure Prediction R-squared: ',num2str(strain_rsq1)));
disp(strcat('25 MPa Strain-to-Failure Prediction Adjusted R-squared: ',num2str(strain_rsq2)));
%%
%50 mpa testing
targets1_50=o_50t_time_tf;
outputs1_50 = sim(net1,input_50t_master1');
errors1_50 = gsubtract(outputs1_50,targets1_50);
performance1_50 = perform(net1,targets1_50,outputs1_50);

f31=figure(31);
title('50 MPa'); hold on;
p50=plot(time50, targets1_50, 'o', 'markersize',4,'color',f_cmap12(1,:),'linewidth',1.5); hold on;
test_p50=plot(time50, outputs1_50, 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',2.5); hold on;
axis tight; xlabel('Axial Strain'); ylabel ('Output: Time to Failure'); 
yyaxis right
% set(gca,'XAxisLocation','top');
test_s50=plot(time50, data50(:,39)./10^6,'linewidth',2); hold on;
ylabel('Axial Stress (MPa)'); hold on;
box on

legend('Testing Set', 'Predictions','location','southeast');
set(gcf,'color','w');
export_fig f31 ANN_Time_test50 -pdf -png -eps -tiff -dsvg -q101 -nocrop

strain_rsq11=rsquare(targets1_50, outputs1_50');
strain_rsq22=rsquare(targets1_50, outputs1_50',false);
disp(strcat('50 MPa Strain-to-Failure Prediction R-squared: ',num2str(strain_rsq11)));
disp(strcat('50 MPa Strain-to-Failure Prediction Adjusted R-squared: ',num2str(strain_rsq22)));
%%
% cmap1=linspecer(3);
% cmap2=linspecer(2);
% 
% f33=figure(33);
% title('Training and Testing (0-45 MPa)'); hold on;
% p0=plot(strain_all_train(cp_train0), strain_tf_train (cp_train0), 'o', 'markersize',4,'color',cmap1(1,:)); hold on;
% plot(strain_all_test(cp_test0), strain_tf_test (cp_test0), 'o', 'markersize',4,'color',cmap1(1,:)); hold on;
% train_p0=plot(strain_all_train(cp_train0), outputs1_train (cp_train0), 'x', 'markersize',4,'color',cmap1(1,:),'linewidth',1.5); hold on;
% test_p0=plot(strain_all_test(cp_test0), outputs1_test (cp_test0), 's', 'markersize',4,'color',cmap1(1,:),'linewidth',1.5); hold on;
% 
% p20=plot(strain_all_train(cp_train20), strain_tf_train (cp_train20), 'o', 'markersize',4,'color',cmap1(2,:)); hold on;
% plot(strain_all_test(cp_test20), strain_tf_test (cp_test20), 'o', 'markersize',4,'color',cmap1(2,:)); hold on;
% train_p20=plot(strain_all_train(cp_train20), outputs1_train (cp_train20), 'x', 'markersize',4,'color',cmap1(2,:),'linewidth',1.5); hold on;
% test_p20=plot(strain_all_test(cp_test20), outputs1_test (cp_test20), 's', 'markersize',4,'color',cmap1(2,:),'linewidth',1.5); hold on;
% 
% p45=plot(strain_all_train(cp_train45), strain_tf_train (cp_train45), 'o', 'markersize',4,'color',cmap1(3,:)); hold on;
% plot(strain_all_test(cp_test45), strain_tf_test (cp_test45), 'o', 'markersize',4,'color',cmap1(3,:)); hold on;
% train_p45=plot(strain_all_train(cp_train45), outputs1_train (cp_train45), 'x', 'markersize',4,'color',cmap1(3,:),'linewidth',1.5); hold on;
% test_p45=plot(strain_all_test(cp_test45), outputs1_test (cp_test45), 's', 'markersize',4,'color',cmap1(3,:),'linewidth',1.5); hold on;
% 
% legend([p0, p20, p45], '0 MPa', '20 MPa', '45 MPa','location','southeast'); hold on;
% 
% axis tight; xlabel('Axial Strain'); 
% box on
% ylabel ('Strain to Failure'); 
% set(gcf,'color','w');
% saveas(gcf,'strain_tf','fig');
% 
% 
% 
% export_fig f33 strain_test3 -pdf -png -eps -tiff -dsvg -q101 -nocrop
%%
strain_w1=net1.IW{1};
strain_w2=net1.LW{2};
strain_b1=net1.b{1};
strain_b2=net1.b{2};

[no_neurons, no_inputs]=size(strain_w1);       %changed format based on DTSpaper
ri_mat1=zeros(size(strain_w1));

for i=1:no_inputs
    for j=1:no_neurons
        ri_mat1(j,i)=abs(strain_w1(j,i).*strain_w2(j));
    end
end
ri_divide=sum(ri_mat1,2);

ri_mat2=ri_mat1;
for i=1:no_neurons
    ri_mat2(i,:)=ri_mat2(i,:)./ri_divide(i);
end

% for i=1:no_neurons
%     for j=1:no_inputs
%         ri_mat3(i,j)= ri_mat1(i,j)./ ri_divide(j);
%     end
% end

ri_strain=sum(ri_mat2); 
ri_strain_cumsum=sum(ri_strain); 
ri_strain2=(ri_strain./ri_strain_cumsum).*100;

f84=figure(84);
subplot(2,1,1)
% bb1=bar(bond_master2(:,2),[bond_master2(:,4) bond_master2(:,6) bond_master2(:,8)],'stack');colormap(f1_color);
bb1=bar(ri_strain2);
hold on
set(bb1,'EdgeColor','k','linewidth',0.25);
hold on
ylabel('Relative Importance'); 
title('Strain-to-Failure Prediction');
xticklabels({'CP', 'MC_{var}', 'SF', 'b', 'D'});
hold on
%%
time_sc=input('Satisfied with Strain-to-Failure Predictions? (1=Yes): ');
if time_sc==1
    save net1;
else
    close all;
end

end
%%
% view(net1);
%plot residuals and biases
dlmwrite('strain_net_w1.dat',strain_w1);
dlmwrite('strain_net_w2.dat',strain_w2);

f11=figure(11); 
plotwb(net1);
set(gcf,'color','w');
export_fig f11 ANN_Strain_Hilton -pdf -png -eps -tiff -dsvg -q101 -nocrop
%stepwise regression fit
[b1,se1,pval1,inmodel1,stats1,nextstep1,history1] = stepwisefit(idata1',targets1', 'penter',0.05,'premove',0.10,'scale','on');
strain_swf=[b1, se1, pval1, inmodel1'];
dlmwrite('strain_swf.dat', strain_swf);
%C1 = coefficients 
%C2 = std error
%C3 = pvalue
%C4 = Inmodel (1=yes)

%5 neurons
%6 inputs 

dlmwrite('ri_strain_matrix1.dat', ri_mat1);
dlmwrite('ri_strain_matrix2.dat', ri_mat2);
dlmwrite('ri_strain.dat', ri_strain);
dlmwrite('ri_strain2.dat', ri_strain2);

%%
while stress_sc~=1
close all
hiddenlayersize2=5;
net2=fitnet(hiddenlayersize2);

%setup division of data
net2.divideParam.trainRatio = 70/100;
net2.divideParam.valRatio = 15/100;
net2.divideParam.testRatio = 15/100;
net2.trainFcn = 'trainlm';
net2.performfcn='sse';
% net2.layers{1}.transferfcn='tansig';

[net2,tr2] = train(net2,idata2,targets2);

%training performance
outputs2 = net2(idata2);
errors2 = gsubtract(outputs2,targets2);
performance2 = perform(net2,targets2,outputs2);

trOut2 = outputs2(tr2.trainInd);
tsOut2 = outputs2(tr2.testInd);
trVal2 = outputs2(tr2.valInd);
tsVal2 = outputs2(tr2.valInd);
trTarg2 = targets2(tr2.trainInd);
tsTarg2 = targets2(tr2.testInd);

% view(net2);

f2_r=figure(4);
plotregression(trTarg2, trOut2, 'Train', tsTarg2, tsOut2, 'Testing',targets2,outputs2,'All');
set(gcf,'color','w');
export_fig f2_r ANN_stress_regression -pdf -png -eps -tiff -dsvg -q101 -nocrop
%%
%locate indexes for plotting
strain_all_train2=strain_all(tr2.trainInd);
strain_all_val2=strain_all(tr2.valInd);
% strain_all_train2=[strain_all_train2; strain_all_val2];
strain_all_test2=strain_all(tr2.testInd);
strain_all_test=[strain_all_test; strain_all_val];

strain_tf_train2=targets2(tr2.trainInd)';
strain_tf_val2=targets2(tr2.valInd)';
% strain_tf_train2=[strain_tf_train2; strain_tf_val2];
strain_tf_test2=targets2(tr2.testInd)';
strain_tf_test=[strain_tf_test; strain_tf_val];

cp2_train=cp_all(tr2.trainInd);
cp2_val=cp_all(tr2.valInd);
% cp2_train=[cp2_train; cp2_val];
cp2_test=cp_all(tr2.testInd);
cp_test=[cp_test; cp_val];

cp2_train0=find(cp2_train==0);
cp2_test0=find(cp2_test==0);
cp2_train5=find(cp2_train==5);
cp2_test5=find(cp2_test==5);
cp2_train10=find(cp2_train==10);
cp2_test10=find(cp2_test==10);
cp2_train15=find(cp2_train==15);
cp2_test15=find(cp2_test==15);
cp2_train20=find(cp2_train==20);
cp2_test20=find(cp2_test==20);
cp2_train30=find(cp2_train==30);
cp2_test30=find(cp2_test==30);
cp2_train35=find(cp2_train==35);
cp2_test35=find(cp2_test==35);
cp2_train40=find(cp2_train==40);
cp2_test40=find(cp2_test==40);
cp2_train45=find(cp2_train==45);
cp2_test45=find(cp2_test==45);
cp2_train50=find(cp2_train==50);
cp2_test50=find(cp2_test==50);
%%
%training plot
outputs2_train=outputs2(tr2.trainInd)';
outputs2_val=outputs2(tr2.valInd)';
% outputs2_train=[outputs2_train; outputs2_val];
outputs2_test=outputs2(tr2.testInd)';
outputs1_test=[outputs1_test; outputs1_val];
%%
f5=figure(5);
subplot(3,3,1);
title('0 MPa'); hold on;
p0=plot(strain_all_train2(cp2_train0), strain_tf_train2 (cp2_train0), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test2(cp2_test0), strain_tf_test2 (cp2_test0), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p0=plot(strain_all_train2(cp2_train0), outputs2_train (cp2_train0), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p0=plot(strain_all_test2(cp2_test0), outputs2_test (cp2_test0), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); 
ylabel ('Stress to Failure (MPa)'); 
box on

subplot(3,3,2);
title('5 MPa'); hold on;
p5=plot(strain_all_train2(cp2_train5), strain_tf_train2 (cp2_train5), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test2(cp2_test5), strain_tf_test2 (cp2_test5), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p5=plot(strain_all_train2(cp2_train5), outputs2_train (cp2_train5), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p5=plot(strain_all_test2(cp2_test5), outputs2_test (cp2_test5), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); 
% ylabel ('Stress to Failure (MPa)'); 
box on

subplot(3,3,3);
title('10 MPa'); hold on;
p10=plot(strain_all_train2(cp2_train10), strain_tf_train2 (cp2_train10), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test2(cp2_test10), strain_tf_test2 (cp2_test10), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p10=plot(strain_all_train2(cp2_train10), outputs2_train (cp2_train10), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p10=plot(strain_all_test2(cp2_test10), outputs2_test (cp2_test10), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); 
% ylabel ('Stress to Failure (MPa)'); 
box on

subplot(3,3,4);
title('15 MPa'); hold on;
p15=plot(strain_all_train2(cp2_train15), strain_tf_train2 (cp2_train15), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test2(cp2_test15), strain_tf_test2 (cp2_test15), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p15=plot(strain_all_train2(cp2_train15), outputs2_train (cp2_train15), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p15=plot(strain_all_test2(cp2_test15), outputs2_test (cp2_test15), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); 
% ylabel ('Stress to Failure (MPa)'); 
box on

subplot(3,3,5);
title('20 MPa'); hold on;
p20=plot(strain_all_train2(cp2_train20), strain_tf_train2 (cp2_train20), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test2(cp2_test20), strain_tf_test2 (cp2_test20), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p20=plot(strain_all_train2(cp2_train20), outputs2_train (cp2_train20), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p20=plot(strain_all_test2(cp2_test20), outputs2_test (cp2_test20), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); ylabel ('Stress to Failure (MPa)'); 
box on

subplot(3,3,6);
title('30 MPa'); hold on;
p30=plot(strain_all_train2(cp2_train30), strain_tf_train2 (cp2_train30), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test2(cp2_test30), strain_tf_test2 (cp2_test30), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p30=plot(strain_all_train2(cp2_train30), outputs2_train (cp2_train30), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p30=plot(strain_all_test2(cp2_test30), outputs2_test (cp2_test30), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); 
% ylabel ('Stress to Failure (MPa)'); 
box on

subplot(3,3,7);
title('35 MPa'); hold on;
p35=plot(strain_all_train2(cp2_train35), strain_tf_train2 (cp2_train35), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test2(cp2_test35), strain_tf_test2 (cp2_test35), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p35=plot(strain_all_train2(cp2_train35), outputs2_train (cp2_train35), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p35=plot(strain_all_test2(cp2_test35), outputs2_test (cp2_test35), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); 
% ylabel ('Stress to Failure (MPa)'); 
box on

subplot(3,3,8);
title('40 MPa'); hold on;
p40=plot(strain_all_train2(cp2_train40), strain_tf_train2 (cp2_train40), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test2(cp2_test40), strain_tf_test2 (cp2_test40), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p40=plot(strain_all_train2(cp2_train40), outputs2_train (cp2_train40), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p40=plot(strain_all_test2(cp2_test40), outputs2_test (cp2_test40), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); 
% ylabel ('Stress to Failure (MPa)'); 
box on

subplot(3,3,9);
title('45 MPa'); hold on;
p45=plot(strain_all_train2(cp2_train45), strain_tf_train2 (cp2_train45), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
plot(strain_all_test2(cp2_test45), strain_tf_test2 (cp2_test45), 'o', 'markersize',2,'color',f_cmap12(1,:)); hold on;
train_p45=plot(strain_all_train2(cp2_train45), outputs2_train (cp2_train45), 'x', 'markersize',4,'color',f_cmap12(2,:),'linewidth',1.5); hold on;
test_p45=plot(strain_all_test2(cp2_test45), outputs2_test (cp2_test45), 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',1.5); hold on;
axis tight; xlabel('Time (s)'); 
ylabel ('Stress to Failure (MPa)'); 
box on

set(gcf,'color','w');
saveas(gcf,'strain_tf','fig');
export_fig f5 ANN_stress_time_tvt -pdf -png -eps -tiff -dsvg -q101 -nocrop
%%
%25 mpa testing
targets2_25=o_25t_stress_tf;
outputs2_25 = sim(net2,input_25t_master1');
errors2_25 = gsubtract(outputs2_25,targets2_25);
performance2_25 = perform(net2,targets2_25,outputs2_25);

f6=figure(6);
title('25 MPa'); hold on;
p25=plot(time25, targets2_25, 'o', 'markersize',4,'color',f_cmap12(1,:),'color','b'); hold on;
test_p25=plot(time25, outputs2_25, 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',2.5); hold on;
axis tight; xlabel('Time (s)'); ylabel ('Output: Axial Stress to Failure (MPa)'); 
yyaxis right
test_s25=plot(time25, data25(:,39)./10^6,'linewidth',2); hold on;
ylabel('Axial Stress (MPa)'); hold on;
box on
legend('Testing Set', 'Predictions','location','southeast');
set(gcf,'color','w');
export_fig f6 ANN_Stress_time_test25 -pdf -png -eps -tiff -dsvg -q101 -nocrop

stress_rsq1=rsquare(targets2_25, outputs2_25');
stress_rsq2=rsquare(targets2_25, outputs2_25',false);
disp(strcat('25 MPa Stress-to-Failure Prediction R-squared: ',num2str(stress_rsq1)));
disp(strcat('25 MPa Stress-to-Failure Prediction Adjusted R-squared: ',num2str(stress_rsq2)));

%%
%50 mpa testing
targets2_50=o_50t_stress_tf;
outputs2_50 = sim(net2,input_50t_master1');
errors2_50 = gsubtract(outputs2_50,targets2_50);
performance2_50 = perform(net2,targets2_50,outputs2_50);

f61=figure(61);
title('50 MPa'); hold on;
p50=plot(time50, targets2_50, 'o', 'markersize',4,'color',f_cmap12(1,:),'color','b'); hold on;
test_p50=plot(time50, outputs2_50, 's', 'markersize',4,'color',f_cmap12(3,:),'linewidth',2.5); hold on;
axis tight; xlabel('Time (s)'); ylabel ('Output: Axial Stress to Failure (MPa)'); 
yyaxis right
test_s50=plot(time50, data50(:,39)./10^6,'linewidth',2); hold on;
ylabel('Axial Stress (MPa)'); hold on;
box on
legend('Testing Set', 'Predictions','location','southeast');
set(gcf,'color','w');
export_fig f61 ANN_Stress_time_test50 -pdf -png -eps -tiff -dsvg -q101 -nocrop

stress_rsq1=rsquare(targets2_50, outputs2_50');
stress_rsq2=rsquare(targets2_50, outputs2_50',false);
disp(strcat('50 MPa Stress-to-Failure Prediction R-squared: ',num2str(stress_rsq1)));
disp(strcat('50 MPa Stress-to-Failure Prediction Adjusted R-squared: ',num2str(stress_rsq2)));

%%
stress_w1=net2.IW{1};
stress_w2=net2.LW{2};
stress_b1=net2.b{1};
stress_b2=net2.b{2};

ri_mat3=zeros(size(stress_w1));

for i=1:no_inputs
    for j=1:no_neurons
        ri_mat3(j,i)=abs(stress_w1(j,i).*stress_w2(j));
    end
end
ri_divide2=sum(ri_mat3,2);

ri_mat4=ri_mat3;
for i=1:no_neurons
    ri_mat4(i,:)=ri_mat4(i,:)./ri_divide2(i);
end

% for i=1:no_neurons
%     for j=1:no_inputs
%         ri_mat3(i,j)= ri_mat1(i,j)./ ri_divide(j);
%     end
% end
ri_stress=sum(ri_mat4);
ri_stress_cumsum=sum(ri_stress);
ri_stress2=(ri_stress./ri_strain_cumsum).*100;

f94=figure(94);
subplot(2,1,1)
% bb1=bar(bond_master2(:,2),[bond_master2(:,4) bond_master2(:,6) bond_master2(:,8)],'stack');colormap(f1_color);
bb1=bar(ri_stress2);
hold on
set(bb1,'EdgeColor','k','linewidth',0.25);
hold on
ylabel('Relative Importance'); 
title('Stress-to-Failure Prediction');
xticklabels({'CP', 'MC_{var}', 'SF', 'b', 'D'});
hold on

stress_sc=input('Satisfied with Stress-to-Failure Predictions? (1=Yes): ');

%%
if stress_sc==1
    save net2;
end

end

% view(net2);
dlmwrite('stress_net_w1.dat',stress_w1);
dlmwrite('stress_net_w2.dat',stress_w2);

f22=figure(22); 
plotwb(net2);
set(gcf,'color','w');
export_fig f22 ANN_Stress_Hilton -pdf -png -eps -tiff -dsvg -q101 -nocrop

[b2,se2,pval2,inmodel2,stats2,nextstep2,history2] = stepwisefit(idata2',targets2', 'penter',0.05,'premove',0.10,'scale','on');
stress_swf=[b2, se2, pval2, inmodel2'];
dlmwrite('strain_swf.dat', stress_swf);
%C1 = coefficients 
%C2 = std error
%C3 = pvalue
%C4 = Inmodel (1=yes)


dlmwrite('ri_stress_matrix1.dat', ri_mat3);
dlmwrite('ri_stress_matrix2.dat', ri_mat4);
dlmwrite('ri_stress.dat', ri_stress);
dlmwrite('ri_stress2.dat', ri_stress2);

%%
f24=figure(24);
subplot(2,1,1)
% bb1=bar(bond_master2(:,2),[bond_master2(:,4) bond_master2(:,6) bond_master2(:,8)],'stack');colormap(f1_color);
bb1=bar(ri_stress2);
hold on
set(bb1,'EdgeColor','k','linewidth',0.25);
hold on
ylabel('Relative Importance (%)'); 
title('Stress-to-Failure Prediction');
xticklabels({'CP', 'MC_{var}', 'SF', 'b', 'D'});
hold on

subplot(2,1,2)
% bb1=bar(bond_master2(:,2),[bond_master2(:,4) bond_master2(:,6) bond_master2(:,8)],'stack');colormap(f1_color);
bar(ri_strain2);
hold on
set(bb1,'EdgeColor','k','linewidth',0.25);
hold on
ylabel('Relative Importance (%)'); 
title('Strain-to-Failure Prediction');
xticklabels({'CP', 'MC_{var}', 'SF', 'b', 'D'});
hold on

set(gcf,'Color','w'); 
export_fig f24 f24_RI -q101 -painters -nocrop -pdf -png -tiff -eps 

