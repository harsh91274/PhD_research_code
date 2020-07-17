%master2
%1 = index
%2 = axial strain
%3 = time
%4 = stress
%5 = target time to failure
%6 = target stress to failure
%7 = TTF index (0 = no prediction, 1 = testing, 2 = training, 3 = validation)
%8 = STF index (0 = no prediction, 1 = testing, 2 = training, 3 = validation)
%9 = TTF prediction 
%10 = STF prediction

ss_data=xlsread('p3_data_binned4.xlsx','initial_S1');
i_data=xlsread('p3_data_binned4.xlsx','wi');

ttf0_all=[i_data(:,4); time_tf0];
ttf5_all=[i_data(:,7); time_tf5];
ttf10_all=[i_data(:,10); time_tf10];
ttf15_all=[i_data(:,13); time_tf15];
ttf20_all=[i_data(:,16); time_tf20];
ttf25_all=[i_data(:,19); time_tf25];
ttf30_all=[i_data(:,22); time_tf30];
ttf35_all=[i_data(:,25); time_tf35];
ttf40_all=[i_data(:,28); time_tf40];
ttf45_all=[i_data(:,31); time_tf45];
ttf50_all=[i_data(:,34); time_tf50];

stf0_all=[i_data(:,5); stress_tf0];
stf5_all=[i_data(:,8); stress_tf5];
stf10_all=[i_data(:,11); stress_tf10];
stf15_all=[i_data(:,14); stress_tf15];
stf20_all=[i_data(:,17); stress_tf20];
stf25_all=[i_data(:,20); stress_tf25];
stf30_all=[i_data(:,23); stress_tf30];
stf35_all=[i_data(:,26); stress_tf35];
stf40_all=[i_data(:,29); stress_tf40];
stf45_all=[i_data(:,32); stress_tf45];
stf50_all=[i_data(:,35); stress_tf50];
%%
index_pr=zeros(10,1);
index_btest=zeros(101,1); index_btest(:)=3;

flag1=zeros(819,1);
flag1(tr1.trainInd)=1;
% flag1(tr1.valInd)=2;
flag1(tr1.valInd)=3;
flag1(tr1.testInd)=3;
flag11=reshape(flag1,[91 9]);

ttf_flag0=[index_pr; flag11(:,1)];
ttf_flag5=[index_pr; flag11(:,2)];
ttf_flag10=[index_pr; flag11(:,3)];
ttf_flag15=[index_pr; flag11(:,4)];
ttf_flag20=[index_pr; flag11(:,5)];
ttf_flag25=index_btest;
ttf_flag30=[index_pr; flag11(:,6)];
ttf_flag35=[index_pr; flag11(:,7)];
ttf_flag40=[index_pr; flag11(:,8)];
ttf_flag45=[index_pr; flag11(:,9)];
ttf_flag50=index_btest;

flag2=zeros(819,1);
flag2(tr2.trainInd)=1;
% flag2(tr2.valInd)=2;
flag2(tr2.valInd)=3;
flag2(tr2.testInd)=3;
flag22=reshape(flag2,[91 9]);

stf_flag0=[index_pr; flag22(:,1)];
stf_flag5=[index_pr; flag22(:,2)];
stf_flag10=[index_pr; flag22(:,3)];
stf_flag15=[index_pr; flag22(:,4)];
stf_flag20=[index_pr; flag22(:,5)];
stf_flag25=index_btest;
stf_flag30=[index_pr; flag22(:,6)];
stf_flag35=[index_pr; flag22(:,7)];
stf_flag40=[index_pr; flag22(:,8)];
stf_flag45=[index_pr; flag22(:,9)];
stf_flag50=index_btest;

preds_ind=nan(10,1);
ttf_preds=reshape(outputs1,[91 9]);
ttf_pred0=[preds_ind; ttf_preds(:,1)];
ttf_pred5=[preds_ind; ttf_preds(:,2)];
ttf_pred10=[preds_ind; ttf_preds(:,3)];
ttf_pred15=[preds_ind; ttf_preds(:,4)];
ttf_pred20=[preds_ind; ttf_preds(:,5)];
ttf_pred25=[preds_ind; outputs1_25'];
ttf_pred30=[preds_ind; ttf_preds(:,6)];
ttf_pred35=[preds_ind; ttf_preds(:,7)];
ttf_pred40=[preds_ind; ttf_preds(:,8)];
ttf_pred45=[preds_ind; ttf_preds(:,9)];
ttf_pred50=[preds_ind; outputs1_50'];

stf_preds=reshape(outputs2,[91 9]);
stf_pred0=[preds_ind; stf_preds(:,1)];
stf_pred5=[preds_ind; stf_preds(:,2)];
stf_pred10=[preds_ind; stf_preds(:,3)];
stf_pred15=[preds_ind; stf_preds(:,4)];
stf_pred20=[preds_ind; stf_preds(:,5)];
stf_pred25=[preds_ind; outputs2_25'];
stf_pred30=[preds_ind; stf_preds(:,6)];
stf_pred35=[preds_ind; stf_preds(:,7)];
stf_pred40=[preds_ind; stf_preds(:,8)];
stf_pred45=[preds_ind; stf_preds(:,9)];
stf_pred50=[preds_ind; outputs2_50'];
%%
master0=[ss_data(:,1), ss_data(:,2), ss_data(:,3), ttf0_all, stf0_all, ttf_flag0, stf_flag0, ttf_pred0, stf_pred0];
master5=[ss_data(:,1), ss_data(:,2), ss_data(:,4), ttf5_all, stf5_all, ttf_flag5, stf_flag5, ttf_pred5, stf_pred5];
master10=[ss_data(:,1), ss_data(:,2), ss_data(:,5), ttf10_all, stf10_all, ttf_flag10, stf_flag10, ttf_pred10, stf_pred10];
master15=[ss_data(:,1), ss_data(:,2), ss_data(:,6), ttf15_all, stf15_all, ttf_flag15, stf_flag15, ttf_pred15, stf_pred15];
master20=[ss_data(:,1), ss_data(:,2), ss_data(:,7), ttf20_all, stf20_all, ttf_flag20, stf_flag20, ttf_pred20, stf_pred20];
master25=[ss_data(:,1), ss_data(:,2), ss_data(:,8), ttf25_all, stf25_all, ttf_flag25, stf_flag25, ttf_pred25, stf_pred25];
master30=[ss_data(:,1), ss_data(:,2), ss_data(:,9), ttf30_all, stf30_all, ttf_flag30, stf_flag30, ttf_pred30, stf_pred30];
master35=[ss_data(:,1), ss_data(:,2), ss_data(:,10), ttf35_all, stf35_all, ttf_flag35, stf_flag35, ttf_pred35, stf_pred35];
master40=[ss_data(:,1), ss_data(:,2), ss_data(:,11), ttf40_all, stf40_all, ttf_flag40, stf_flag40, ttf_pred40, stf_pred40];
master45=[ss_data(:,1), ss_data(:,2), ss_data(:,12), ttf45_all, stf45_all, ttf_flag45, stf_flag45, ttf_pred45, stf_pred45];
master50=[ss_data(:,1), ss_data(:,2), ss_data(:,13), ttf50_all, stf50_all, ttf_flag50, stf_flag50, ttf_pred50, stf_pred50];

masta1=[master0; master5; master10; master15; master20; master25; master30; master35; master40; master45; master50];

[lenn,~]=size(masta1);
masta_index=1:lenn;
index=masta_index';

masta2=[index, masta1];

f72=figure(72);
set(gcf,'OuterPosition',[159 78 1080 686]);

s1=subplot(3,1,1);
plot(masta2(:,1),masta2(:,4)./10^6,'k', 'linewidth', 1);
xlim([0 max(masta2(:,1))]);
xticks([]);
xlabel('Time (s)'); 
ylabel ('Axial Stress (MPa)'); 
box on

s2=subplot(3,1,2);
plot(masta2(:,1),masta2(:,5),'k', 'linewidth', 1);hold on;
xlim([0 max(masta2(:,1))]);
xticks([]);
for i=1:lenn
    if masta2(i,7)==1
        plot(masta2(i,1), masta2(i,9),'x','color','b','linewidth',0.9); hold on;
    elseif masta2(i,7)==2
        plot(masta2(i,1), masta2(i,9),'s','color','g','linewidth',1.2); hold on; 
    elseif masta2(i,7)==3
        plot(masta2(i,1), masta2(i,9),'o','color','r','linewidth',0.9); hold on;
    end
end

xlabel('Time (s)'); 
% ylim([-100 1500]);
ylabel ('Time to Failure (s)'); 
box on

s3=subplot(3,1,3);
plot(masta2(:,1),masta2(:,6),'k', 'linewidth', 1); hold on; 
% xticks([0 101 202 303 404 505 606 707 808 909 1010 1111]);
xlim([0 max(masta2(:,1))]);
xticks([]);
for i=1:lenn
    if masta2(i,8)==1
        plot(masta2(i,1), masta2(i,10),'x','color','b','linewidth',0.9); hold on;
    elseif masta2(i,8)==2
        plot(masta2(i,1), masta2(i,10),'s','color','g','linewidth',1.2); hold on; 
    elseif masta2(i,8)==3
        plot(masta2(i,1), masta2(i,10),'o','color','r','linewidth',0.9); hold on;
    end
end

xlabel('Time (s)'); 
% ylim([-200 100]);
ylabel ('Stress to Failure (MPa)'); 
box on

set(gcf,'color','w');
export_fig f72 f72_ANN_predictions -pdf -png -eps -tiff -dsvg -q101 -nocrop
