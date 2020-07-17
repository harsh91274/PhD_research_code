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

eqdata=importdata('test2.dat'); 
%1 = cycle
%2 - 7 - seismometer displacement 
%8 = cycles from mainshock
f_cmap=linspecer(6);
%%
f1=figure; 
for i=1:6
    plot(eqdata(:,1),eqdata(:,2*i),'color',f_cmap(i,:)); hold on;
end
legend ('S1x','S2x','S3x','S4x','S5x','S6x','location','best');

f2=figure; 
for i=1:6
    plot(eqdata(:,1),eqdata(:,1+2*i),'color',f_cmap(i,:)); hold on;
end
legend ('S1y','S2y','S3y','S4y','S5y','S6y','location','best');


% ev=[43 829 1162];                         %cycle for onset of EQ
% cyc_data=eqdata(:,1);
% targets4_mod=zeros(size(cyc_data));
% counter=0;
% 
% for i=1:length(ev)
%     for j=counter+1:ev(i)
%         if i==1
%             targets4_mod(j)=ev(i)-cyc_data(j);
%         else
%             targets4_mod(j)=ev(i)-cyc_data(j)-ev(i-1);
%         end
%     end
%     counter=ev(i);
% end
% 
% figure, plot(targets4_mod); hold on; plot(eqdata(:,14),'ro'); hold on;

targets4=eqdata(:,14)';
idata4=eqdata(:,1:13)';
image=idata4(1,:);

%%
%random data division
hiddenLayerSize4 = [200];
net4 = fitnet(hiddenLayerSize4);

net4.divideParam.trainRatio = 70/100;
net4.divideParam.valRatio = 15/100;
net4.divideParam.testRatio = 15/100;
net4.trainFcn = 'trainlm';

[net4,tr4] = train(net4,idata4,targets4);

outputs4 = net4(idata4);
errors4 = gsubtract(outputs4,targets4);
performance4 = perform(net4,targets4,outputs4);

trOut4 = outputs4(tr4.trainInd);
tsOut4 = outputs4(tr4.testInd);
vout4=outputs4(tr4.valInd);
vtarg4=targets4(tr4.valInd);
trTarg4 = targets4(tr4.trainInd);
tsTarg4 = targets4(tr4.testInd);

figure, ploterrhist(errors4);

f4_r=figure;
plotregression(trTarg4, trOut4, 'Train',vtarg4, vout4, 'Validation' , tsTarg4, tsOut4, 'Testing',targets4,outputs4,'All');

error_ratio4=errors4./targets4;
all4=length(error_ratio4);
positive4=length(find(error_ratio4<=0.1));
success4=(positive4/all4)*100;
true_error4=norm(targets4-outputs4,2)/norm(targets4,2);

f4_fit=figure;
plot(image(:),targets4); hold on;
plot(image(tr4.trainInd), trOut4, 'o','markersize',6); hold on;
plot(image(tr4.testInd), tsOut4, 'x','markersize',6); hold on;
legend('Images to failure','Training', 'Testing')

disp(strcat('Image to EQ Prediction - Percentage of Success (<10% Error): ',num2str(success4),' %'));
disp(strcat('Image to EQ Prediction - True Error): ',num2str(true_error4),' %'));
%%
%index data division

%%
%random data division
hiddenLayerSize42 = [200];
net42 = fitnet(hiddenLayerSize42);

idata42=idata4;
% idata_42=

net42.divideFcn='divideind';
net42.divideParam.trainInd = 1:813;
net42.divideParam.valInd = 814:987;
net42.divideParam.testInd =  988:1162;
net42.trainFcn = 'trainlm';

[net42,tr42] = train(net42,idata42,targets4);

outputs42 = net42(idata42);
errors42 = gsubtract(outputs42,targets4);
performance42 = perform(net42,targets4,outputs42);

trOut42 = outputs42(tr42.trainInd);
tsOut42 = outputs42(tr42.testInd);
vout42=outputs4(tr42.valInd);
vtarg42=targets4(tr42.valInd);
trTarg42 = targets4(tr42.trainInd);
tsTarg42 = targets4(tr42.testInd);

figure, ploterrhist(errors42);

f42_r=figure;
plotregression(trTarg42, trOut42, 'Train',vtarg42, vout42, 'Validation' , tsTarg42, tsOut42, 'Testing',targets4,outputs42,'All');

error_ratio42=errors42./targets4;
all42=length(error_ratio42);
positive42=length(find(error_ratio42<=0.1));
success42=(positive42/all42)*100;
true_error42=norm(targets4-outputs42,2)/norm(targets4,2);

f42_fit=figure;
plot(image(:),targets4); hold on;
plot(image(tr42.trainInd), trOut42, 'o','markersize',6); hold on;
plot(image(tr42.testInd), tsOut42, 'x','markersize',6); hold on;
legend('Images to failure','Training', 'Testing')

disp(strcat('Image to EQ Prediction - Percentage of Success (<10% Error): ',num2str(success42),' %'));
disp(strcat('Image to EQ Prediction - True Error): ',num2str(true_error42),' %'));

