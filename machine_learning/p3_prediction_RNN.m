strain_counter=0;

while strain_counter~=1
close all

net1=layrecnet(1:2,10);

X = con2seq(idata1);
T = con2seq(targets1);

net1.divideParam.trainRatio = 70/100;
net1.divideParam.valRatio = 15/100;
net1.divideParam.testRatio = 15/100;
net1.trainFcn = 'trainlm';
net1.performfcn='msesparse';

[Xs,Xi,Ai,Ts] = preparets(net1,X,T);    %Ai and Xi are delays
net1 = train(net1,Xs,Ts,Xi,Ai);
view(net1);
Y = net1(Xs,Xi,Ai);
Xs2=cell2mat(Xs); %shifted inputs
Ts2=cell2mat(Ts); %shifted targets
Y2=cell2mat(Y);   
perf = perform(net5,Ts,Y);

Xtrain2=cell2mat(Xtrain);
Xtest2=cell2mat(Xtest);
T2=cell2mat(T);

rnn_predictions=sim(net5,Xtest);
rnn_predictions2=cell2mat(rnn_predictions);
% perf = perform(net5,,Ts);
plot(rnn_predictions2);
plot(Xs2,Ts2);

figure;
plot(image,targets5,'r'); hold on;
plot(Xs2, Y2,'o'); hold on;
plot (image(rat+1:end), rnn_predictions, 'x');  hold on
legend ('Strain to Failure', 'Training Set', 'Testing Set'); 

outputs5=cat(1,Y,rnn_predictions);
true_error_rnn1=norm(targets5-outputs5,2)/norm(targets5,2);

end
