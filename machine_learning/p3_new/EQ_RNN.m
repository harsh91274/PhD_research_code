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
%%
targets5=eqdata(:,14)';
idata5=eqdata(:,2:13)';
image=eqdata(:,1)';

net5 = layrecnet(1:2,20);

[~,len]=size(idata5);
rat=round(len*0.7);

X = con2seq(idata5);
T = con2seq(targets5);

Xtrain=con2seq(idata5(:,1:rat));
Xtest=con2seq(idata5(:,rat+1:end));

Ttrain=con2seq(targets5(:,1:rat));
Ttest=con2seq(targets5(:,rat+1:end));

[Xs,Xi,Ai,Ts] = preparets(net5,Xtrain,Ttrain);
net5 = train(net5,Xs,Ts,Xi,Ai);
view(net5);
Y = net5(Xs,Xi,Ai);
Xs2=cell2mat(Xs);
Ts2=cell2mat(Ts);
Y2=cell2mat(Y);
% perf = perform(net5,Y,Ts);

rnn_predictions=sim(net5,Xtest);
% perf = perform(net5,,Ts);
plot(rnn_predictions,'x');
plot(Xs2,Ts2);

figure;
plot(image,targets5,'r'); hold on;
plot(Xs2, Y2,'o'); hold on;
plot (image(rat+1:end), rnn_predictions, 'x');  hold on
legend ('Images to Failure', 'Training Set', 'Testing Set'); 

outputs5=cat(1,Y,rnn_predictions);
true_error_rnn1=norm(targets5-outputs5,2)/norm(targets5,2);

