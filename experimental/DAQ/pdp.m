function pdp

clc
clear all
close all

file_name='Voltage_data.txt';
data=importdata(file_name);

disp('PDP Test Start Time');
h1=input('Enter Hours: ');
m1=input('Enter minutes: ');
s1=input('Enter seconds: ');
start_time=h1*3600+m1*60+s1;

disp('PDP Test End Time');
h2=input('Enter Hours: ');
m2=input('Enter minutes: ');
s2=input('Enter seconds: ');
end_time=h2*3600+m2*60+s2;

% dt=input('Enter dt: ');
len=length(data);
time=linspace(start_time,end_time,len); %x axis of plot

%Voltage-Pressure Scaling

m1=input('Enter slope for Upstream Sensor (m1): ');
int1=input('Enter intercept for Upstream Sensor (int1): ');
m2=input('Enter slope for Downstream Sensor (m2): ');
int2=input('Enter intercept for Downstream Sensor (int2): ');

pdata=zeros(size(data));
pdata(:,1)=m1.*data(:,1)+int1;
pdata(:,2)=m2.*data(:,2)+int2;
p1_data=pdata(:,1);
p2_data=pdata(:,2);

%Plotting and limits   
figure
plot(time,pdata(:,1),'r','LineWidth',5);
hold on
plot(time,pdata(:,2),'b','LineWidth',5);
hold on
xlabel('Time (seconds)');ylabel('Pressure (psi)');
legend('Upstream Pressure','Downstream Pressure');

%Data Extraction
satisfied=0;
while satisfied==0    
    %Upstream Sensor
    disp('Select limits for Upstream sensor, lower limit first');
    buff1=waitforbuttonpress;
    [u_locs,~]=ginput(2);
    [t_locs1]=find(time>=u_locs(1));
    [t_locs2]=find(time<=u_locs(2));
    [t_locs]=intersect(t_locs1,t_locs2);
    s1_pressures=p1_data(t_locs);
    s1_times=(time(t_locs))';

    %Slope calculation
    coeff1=polyfit(s1_times,s1_pressures,1);
    slope1=coeff1(1);
    disp(strcat('Upstream Slope: ',num2str(slope1)));
    
    %Plot best fit lines
    plot_coeff1=polyval(coeff1,s1_times);
    bf1=plot(s1_times,plot_coeff1,'g--','LineWidth',3);
    
    %Downstream Sensor
    disp('Select limits for Downstream sensor, lower limit first');
    buff2=waitforbuttonpress;
    [d_locs,~]=ginput(2);
    [d_t_locs1]=find(time>=d_locs(1));
    [d_t_locs2]=find(time<=d_locs(2));
    [d_t_locs]=intersect(d_t_locs1,d_t_locs2);
    s2_pressures=p2_data(d_t_locs);
    s2_times=(time(d_t_locs))';

    %Slope calculation
    coeff2=polyfit(s2_times,s2_pressures,1);
    slope2=coeff2(1);
    disp(strcat('Downstream Slope: ',num2str(slope2)));
    
    %Plot best fit lines
    plot_coeff2=polyval(coeff2,s2_times);
    bf2=plot(s2_times,plot_coeff2,'g--','LineWidth',3);
    
    %Iterate ?
    satisfied=input('Satisfied with fit ? (Yes = 1, No = 0) : ');
    if satisfied==0
        delete(bf1);delete(bf2);
    end
end


%%
keyboard
end
