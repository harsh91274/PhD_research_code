ev=[43 829 1162];   %cycle for onset of EQ
cyc_data=eqdata(:,1);
targets4_mod=zeros(size(cyc_data));

counter=0;
counter2=0;
for i=1:length(ev)
    for j=counter+1:ev(i)
            targets4_mod(j)=ev(i)-cyc_data(j);
    end
    counter=ev(i);
end

figure, 
plot(eqdata(:,14),'o'); hold on;
plot(targets4_mod); hold on; 

