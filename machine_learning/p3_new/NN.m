close all
clear all
clc

%paths
path(path,genpath('C:\export_fig\update2'));
path(path,genpath('C:\Program Files\Glyph & Cog\XpdfReader-win64'));
path(path,genpath('C:\kakearney-legendflex-pkg-f29cb4e'));
path(path,genpath('C:\kakearney-legendflex-pkg-f29cb4e\legendflex'));
path(path,genpath('C:\Program Files\gs\gs9.22\bin'));
path(path,genpath('C:\Program Files\MATLAB\R2010a\toolbox\linspecer'));
path(path,genpath('C:\xpdf\xpdf-tools-win-4.00'));
path(path,genpath('C:\xpdf\xpdf-tools-win-4.00\bin32'));
warning('off','all')

inputdata=importdata('test1_in.txt');   %import inputs
outputdata=importdata('test1_out.txt'); %import outputs

idata=inputdata.data';   %data only
odata=outputdata.data;  %data only
f_cmap=linspecer(11);
%%
%predict confining pressure
targets1=odata(:,1)';    %confining pressure
targets2=odata(:,2)';    %stress to failure
targets3=odata(:,3)';    %strain to failure
    
% create a fitting network
% hiddenLayerSize1 = [5 3];
% net1 = fitnet(hiddenLayerSize1);
% 
% %setup division of data
% net1.divideParam.trainRatio = 70/100;
% net1.divideParam.valRatio = 0/100;
% net1.divideParam.testRatio = 15/100;
% net1.trainFcn = 'trainbr';
% % net1.numlayers=3;
% 
% 
% %train network using bayesian regularization
% [net1,tr1] = train(net1,idata,targets1);
% 
% % Test the Network
% outputs1 = net1(idata);
% errors1 = gsubtract(outputs1,targets1);
% performance1 = perform(net1,targets1,outputs1);
% 
% %view the network
% view(net1)
% 
% % Plots
% % figure, plotperform(tr1)
% % figure, plottrainstate(tr1)
% % figure, plotfit(targets1,outputs1)
% figure, ploterrhist(errors1)
% 
% trOut1 = outputs1(tr1.trainInd);
% tsOut1 = outputs1(tr1.testInd);
% trTarg1 = targets1(tr1.trainInd);
% tsTarg1 = targets1(tr1.testInd);
% 
% 
% 
% f1_r=figure;
% plotregression(trTarg1, trOut1, 'Train', tsTarg1, tsOut1, 'Testing',targets1,outputs1,'All');
% set(gcf,'color','w');
% export_fig f1_r NN_CP_regression -pdf -png -eps -tiff -dsvg -q101 -nocrop
% 
% cp_mod0=odata(:,2)';
% cp_mod_train0=odata(tr1.trainInd);
% cp_mod_test0=odata(tr1.testInd);
% 
% f1_fit=figure;
% plot(targets1,'o', 'markersize',6,'linewidth',1.5); hold on;
% plot(tr1.trainInd, trOut1,'x','markersize',6,'linewidth',1.5); hold on;
% plot(tr1.testInd, tsOut1,'*','markersize',6,'linewidth',1.5); hold on;
% legend ('Confining Pressure', 'Training','Testing','location','southeast');
% set(gcf,'color','w');
% xlabel('Index'), ylabel('confining pressure');
% export_fig f1_fit NN_CP_fit -pdf -png -eps -tiff -dsvg -q101 -nocrop
% 
% error_ratio1=errors1./targets1;
% all1=length(error_ratio1);
% positive1=length(find(error_ratio1<=0.1));
% success1=(positive1/all1)*100;
% true_error1=norm(targets1-outputs1,2)/norm(targets1,2);
% 
% disp(strcat('Confining Pressure Prediction - Percentage of Success (<10% Error): ',num2str(success1),' %'));
% disp(strcat('Confining Pressure Prediction - True Error): ',num2str(true_error1),' %'));
% 
% cp_net=net1;
% save cp_net;

%%
%predict peak stress on rock
%create a fitting network
% hiddenLayerSize2 = [13 2];
% net2 = fitnet(hiddenLayerSize2);
% 
% %setup division of data
% net2.divideParam.trainRatio = 70/100;
% net2.divideParam.valRatio = 0/100;
% net2.divideParam.testRatio = 30/100;
% net2.trainFcn = 'trainbr';
% 
% %train network using bayesian regularization
% [net2,tr2] = train(net2,idata,targets2);
% 
% % Test the Network
% outputs2 = net2(idata);
% errors2 = gsubtract(outputs2,targets2);
% performance2 = perform(net2,targets2,outputs2);
% 
% %view the network
% view(net2)
% 
% % Plots
% % figure, plotperform(tr1)
% % figure, plottrainstate(tr1)
% % figure, plotfit(targets1,outputs1)
% figure, ploterrhist(errors2)
% 
% trOut2 = outputs2(tr2.trainInd);
% tsOut2 = outputs2(tr2.testInd);
% trTarg2 = targets2(tr2.trainInd);
% tsTarg2 = targets2(tr2.testInd);
% 
% f2_r=figure;
% plotregression(trTarg2, trOut2, 'Train', tsTarg2, tsOut2, 'Testing',targets2,outputs2,'All');
% 
% set(gcf,'color','w');
% export_fig f2_r NN_Stress_regression -pdf -png -eps -tiff -dsvg -q101 -nocrop
% 
% cp_mod=odata(:,1)';
% cp_mod_train=odata(tr2.trainInd);
% cp_mod_test=odata(tr2.testInd);
% 
% f2_fit=figure;
% plot(cp_mod, targets2,'o', 'markersize',6,'linewidth',1.5); hold on;
% plot(cp_mod_train,trOut2,'x','markersize',6,'linewidth',1.5); hold on;
% plot(cp_mod_test,tsOut2,'*','markersize',6,'linewidth',1.5); hold on;
% xlabel('confining pressure (MPa)'); ylabel ('Axial Stress (MPa)');
% legend ('Stress to Failure', 'Train','Test','location','best');
% 
% set(gcf,'color','w');
% export_fig f2_fit NN_Stress_fit -pdf -png -eps -tiff -dsvg -q101 -nocrop
% 
% error_ratio2=errors2./targets2;
% all2=length(error_ratio2);
% positive2=length(find(error_ratio2<=0.1));
% success2=(positive2/all2)*100;
% true_error2=norm(targets2-outputs2,2)/norm(targets2,2);
% 
% disp(strcat('Stress Prediction - Percentage of Success (<10% Error): ',num2str(success2),' %'));
% disp(strcat('Stress Prediction - True Error): ',num2str(true_error2),' %'))
% 
% stress_net=net2;
% save net2;
%%
%predict strain to failure
%create a fitting network
hiddenLayerSize3 = 500;
net3 = fitnet(hiddenLayerSize3);

%setup division of data
net3.divideParam.trainRatio = 70/100;
net3.divideParam.valRatio = 0/100;
net3.divideParam.testRatio = 30/100;
net3.trainFcn = 'trainbr';
% net3.numlayers=3;

% targets32=targets3+1;
%train network using bayesian regularization
[net3,tr3] = train(net3,idata,targets3);

% Test the Network
outputs3 = net3(idata);
% outputs3=outputs32-1;
errors3 = gsubtract(outputs3,targets3);
% performance3 = perform(net3,targets3,outputs3);

%view the network
%view(net3)

% Plots
% figure, plotperform(tr1)
% figure, plottrainstate(tr1)
% figure, plotfit(targets1,outputs1)
figure, ploterrhist(errors3)

trOut3 = outputs3(tr3.trainInd);
tsOut3 = outputs3(tr3.testInd);
trTarg3 = targets3(tr3.trainInd);
tsTarg3 = targets3(tr3.testInd);

f3_r=figure;
plotregression(trTarg3, trOut3, 'Train', tsTarg3, tsOut3, 'Testing',targets3,outputs3,'All');
set(gcf,'color','w');
export_fig f3_r NN_Strain_regression -pdf -png -eps -tiff -dsvg -q101 -nocrop

cp_mod=odata(:,1)';
cp_mod_train3=odata(tr3.trainInd);
cp_mod_test3=odata(tr3.testInd);

% f2_fit=figure;
% plot(cp_mod, targets2,'o-', 'markersize',6,'linewidth',1.5); hold on;
% plot(cp_mod_train,trOut2,'x','markersize',6,'linewidth',1.5); hold on;
% plot(cp_mod_test,tsOut2,'*','markersize',6,'linewidth',1.5); hold on;
% xlabel('confining pressure (MPa)'); ylabel ('Axial Stress (MPa)');


f3_fit=figure;
plot(cp_mod, targets3,'o','markersize',6,'linewidth',1.5); hold on;
plot(cp_mod_train3,trOut3,'x','markersize',6,'linewidth',1.5); hold on;
plot(cp_mod_test3,tsOut3,'*','markersize',6,'linewidth',1.5); hold on;
xlabel('confining pressure (MPa)'); ylabel ('Axial Strain to Failure (MPa)');
% 
% plot(cp_mod_train3, outputs3,'x','markersize',6,'linewidth',1.5); hold on;
legend ('Strain to Failure', 'Train','Test','location','best');
set(gcf,'color','w');
export_fig f3_fit NN_Strain_fit -pdf -png -eps -tiff -dsvg -q101 -nocrop

error_ratio3=errors3./targets3;
all3=length(error_ratio3);
positive3=length(find(error_ratio3)<=0.1);
success3=(positive3/all3)*100;
true_error3=norm(targets3-outputs3,2)/norm(targets3,2);

disp(strcat('Strain Prediction - Percentage of Success (<10% Error): ',num2str(success3),' %'));
disp(strcat('Strain Prediction - True Error): ',num2str(true_error3),' %'))

strain_net=net3;
save net3;
%%
save('NN_failure_prediction');