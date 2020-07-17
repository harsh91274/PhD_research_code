close all
clear all
clc
inputdata=importdata('test1_in.txt');   %import inputs
outputdata=importdata('test1_out.txt'); %import outputs

idata=inputdata.data;   %data only
odata=outputdata.data;  %data only

odata_cp=odata(:,1);
odata_s1=odata(:,2);
odata_strain=odata(:,3);

%predict cp
idata_scale=max(max(idata));    %inputs scaling factor
odata_cp_scale=max(odata_cp);

idata2=idata./idata_scale;                  %scaled inputs
odata2_cp=odata_cp./odata_cp_scale;      %scaled outputs

[len,~]=size(idata2);   
[train_ind, val_ind, test_ind]=dividerand(len,0.7,0,0.3);   %indices for training and testing sets

idata_train=idata2(train_ind)';  %input training set
idata_test=idata2(test_ind)';    %input testing set
odata_train=odata2_cp(train_ind);  %output training set
odata_test=odata_cp(test_ind);    %output testing set

md_cp=fitensemble(idata_train, odata_train,'LSBOOST', 100,'Tree');
cp_pre=predict(md_cp, idata_test);
cp_predicted=cp_pre.*odata_cp_scale;
cp_error=abs(odata_test-cp_predicted);

f1=figure(1);
plot(odata_test,'o','markersize',4); hold on;
plot(cp_predicted,'*','markersize',4); hold on;
plot(cp_error,'.-'); hold on;
legend('Confining Pressure Test Data', 'Confining Pressure Predicted','Error','Location','best');

%predict stress to failure
odata_s1_scale=max(odata_s1);
odata2_s1=odata_s1./odata_s1_scale;      %scaled outputs

odata_train2=odata2_s1(train_ind);  %output training set
odata_test2=odata_s1(test_ind);    %output testing set

md_s1=fitensemble(idata_train, odata_train2,'LSBOOST', 100,'Tree');

s1_pre=predict(md_s1, idata_test);
s1_predicted=s1_pre.*odata_s1_scale;
s1_error=abs(odata_test2-s1_predicted);

f2=figure(2);
plot(odata_test2,'o','markersize',4); hold on;
plot(s1_predicted,'*','markersize',4); hold on;
plot(s1_error,'.-'); hold on;
legend('Axial Stress to Failure Test Data', 'Axial Stress to Failure Predicted','Error','Location','best');

%predict strain to failure
odata_strain_scale=max(odata_strain);
odata2_strain=odata_strain./odata_strain_scale;      %scaled outputs

odata_train3=odata2_strain(train_ind);  %output training set
odata_test3=odata_strain(test_ind);    %output testing set

md_strain=fitensemble(idata_train, odata_train3,'LSBOOST', 100,'Tree');

strain_pre=predict(md_strain, idata_test);
strain_predicted=strain_pre.*odata_strain_scale;
strain_error=abs(odata_test3-strain_predicted);

f3=figure(3);
plot(odata_test3,'o','markersize',4); hold on;
plot(strain_predicted,'*','markersize',4); hold on;
plot(strain_error,'.-'); hold on;
legend('Axial Strain to Failure Test Data', 'Axial Strain to Failure Predicted','Error','Location','best');

