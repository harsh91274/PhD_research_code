output_bins2=[];
% bin_index=[1 11 21 31 41 51 61 71 81 91];   %bins start no
bin_index=1:91;
s1_buff1=dlmread('s1.dat');
for iiq=1:length(bin_index)
    close all
    output_bins1=[];
    n1=bin_index(iiq);
    n2=bin_index(iiq)+9;
    disp(strcat('CYCLE ',num2str (n1), 'to CYCLE ', num2str(n2)));
    output_bins1=[cp, n1, n2];    
    master2_bin=master2(n1:n2,:);
    bond_master2_bin=bond_master2(n1:n2,:);
    %%
    %S3 energy
    sample_area=newoutput(n1:n2,9);
    sample_vol=sample_area.*zave;
    sample_length=newoutput(n1:n2,4)-newoutput(n1:n2,3);
    area_s3=sample_length.*zave;
    
    y_top_initial=newoutput_i(1,6);
    y_bot_initial=newoutput_i(1,7);
    r_initial=y_top_initial-y_bot_initial;
    y_top=newoutput(n1:n2,6);
    y_bot=newoutput(n1:n2,7);
    r_change=y_top-y_bot;
    
    area_s1=zave*r_initial;
    
    r_strain=(r_change-r_initial)./r_initial;
    
    %S3 energy
    pe=master2_bin(:,8);
    % pe2=master2(:,8).*r_strain;
    pe3=trapz(r_strain,pe);
    
    %S1 energy
    sed=master2_bin(:,9);
    % sed2=master2(:,9).*master2(:,13);
    se=trapz(master2_bin(:,13),sed);
    %
    % ie=sed2+pe2;
    % ie2=ie'*sample_vol;
    
    pe4=master2_bin(:,8).*10^6.*area_s3;    %energy from S3
    pe5=trapz(r_strain,pe4);
    
    
    se4=master2_bin(:,9).*10^6.*area_s1;    %energy from S1
    se5=trapz(master2_bin(:,13),se4);
    
    ie0=se+pe3;
    ie=(pe5+se5);
    
    disp('INPUT ENERGY');
    disp(strcat('Input Strain Energy (from S1): ',num2str(se),' MJ/m2'));
    disp(strcat('Input Strain Energy (from S1): ',num2str(se5),' Joules'));
    disp(strcat('Input Potential Energy (from S3): ',num2str(pe3),' MJ/m2'));
    disp(strcat('Input Potential Energy (from S3): ',num2str(pe5),' Joules'));
    disp(strcat('Total Input Energy: ',num2str(ie0),' MJ/m2'));
    disp(strcat('Total Input Energy: ',num2str(ie),' Joules'));
    
    output_bins1=[output_bins1, se5, pe5, ie, ie0];
    %%
    %microcracking and fracture energy
    n_shear=sum(bond_master2_bin(:,4));
    e_shear=sum(bond_master2_bin(:,5));
    n_ts=sum(bond_master2_bin(:,6));
    e_ts=sum(bond_master2_bin(:,7));
    n_tensile=sum(bond_master2_bin(:,8));
    e_tensile=sum(bond_master2_bin(:,9));
    
    disp('BOND DATA');
    disp(strcat('Number of Shear MF: ',num2str(n_shear)));
    disp(strcat('Number of Tensile-Shear MF: ',num2str(n_ts)));
    disp(strcat('Number of Tensile MF: ',num2str(n_tensile)));
    n_sum=n_shear+n_ts+n_tensile; n_sf=(n_shear+n_ts)/n_sum;
    
    disp(strcat('Energy from Shear MF: ',num2str(e_shear),' J'));
    disp(strcat('Energy from Tensile-Shear MF: ',num2str(e_ts),' J'));
    disp(strcat('Energy from Tensile MF: ',num2str(e_tensile),' J'));
    e_sum=e_shear+e_ts+e_tensile; e_sf=(e_shear+e_ts)/e_sum;
    e_frac0=e_sum/ie;
    disp(strcat('Total Microcracking Energy: ',num2str(e_sum),' J'));
    disp(strcat('Microcracking Energy (J)/Input Energy (J): ',num2str(e_frac0)));
    
    fracture_energy1=master2_bin(:,5)./(sample_vol);
    G=sum(fracture_energy1); efrac1=G/(ie0*10^6);
    % G2=sum(master2(:,18)); efrac2=G2/(ie0*10^6);
    
    % disp(strcat('Total Input Energy: ',num2str(ie0),' MJ/m2'));
    % disp(strcat('Total Input Energy: ',num2str(ie),' Joules')ie0);
    
    disp(strcat('Fracture Energy from area calculations (G1): ',num2str(G),' J/m2'));
    disp(strcat('Input/Output energy fraction (G1/IE): ',num2str(efrac1)));
    % disp(strcat('Fracture Energy from microcrack calculations (G2): ',num2str(G2),' J/m2'));
    % disp(strcat('Input/Output energy fraction (G2/IE): ',num2str(efrac2)));
    
    output_bins1=[output_bins1, n_shear, n_ts, n_tensile, n_sum, n_sf, e_shear, e_ts, e_tensile, e_sum, e_sf, e_frac0, G];
    % output0=[output0, n_shear, n_ts, n_tensile, n_sum, n_sf, e_shear, e_ts,
    % e_tensile, e_sum, e_sf, e_frac0, G];
    
    %%
    %event statistics
    er_std=std(master2_bin(:,2));
    disp(strcat('Event Rate Standard Deviation: ',num2str(er_std)));
    er_var=var(master2_bin(:,2));
    disp(strcat('Event Rate Variance: ',num2str(er_var)));
    er_kur=kurtosis(master2_bin(:,2));
    disp(strcat('Event Rate Kurtosis: ',num2str(er_kur)));
    
    output_bins1=[output_bins1, er_std, er_var, er_kur];
    %event rate
%     fig23=figure;     %uncomment for plot
%     plot(master2_bin(:,13),master2_bin(:,2),'ro','LineWidth',2);     %uncomment for plot
%     hold on
%     xlabel('Axial Strain','fontweight','bold');ylabel('Number of microfractures','fontweight','bold');
%     axis tight

    er_fit=polyfit(master2_bin(:,13),master2_bin(:,2), 1);
    er_fit_plot=polyval(er_fit,master2_bin(:,13));
%     er_plot=plot(master2_bin(:,13),er_fit_plot,'k--','linewidth',3);   %uncomment for plot
    er_val=er_fit(1);
%     [x_range_er,~]=ginput(2);
%     
%     [locs1_er]=find(master2(:,13)>=x_range_er(1));
%     [locs2_er]=find(master2(:,13)<=x_range_er(2));
%     [locs_er]=intersect(locs1_er, locs2_er);
%     
%     xfit_er=master2(locs_er,13);
%     yfit_er=master2(locs_er,2);
%     er_fit=polyfit(xfit_er,yfit_er,1);
%     er_fit_plot=polyval(er_fit,xfit_er);
%     er_plot=plot(xfit_er,er_fit_plot,'k--','linewidth',3);
%     hold on
%     
    
%     legend('Number of Microfractures','best fit line', 'location','best');
    disp(strcat('Event Rate: ',num2str(er_val),' events/ 1% strain'));
    
    output_bins1=[output_bins1, er_val];
    %%
    %event energy statistics
    ee_std=std(master2_bin(:,5));
    disp(strcat('Event Eenergy Standard Deviation: ',num2str(ee_std)));
    ee_var=var(master2_bin(:,5));
    disp(strcat('Event Energy Variance: ',num2str(ee_var)));
    ee_kur=kurtosis(master2_bin(:,5));
    disp(strcat('Event Energy Kurtosis: ',num2str(ee_kur)));

    output_bins1=[output_bins1, ee_std, ee_var, ee_kur];
    
%     fig24=figure; %uncomment for plot
%         plot(master2_bin(:,13),master2_bin(:,5),'ro','LineWidth',2);
%     hold on
%     xlabel('Axial Strain','fontweight','bold');ylabel('Microfracture energy','fontweight','bold');
%     axis tight
    
    ee_fit=polyfit(master2_bin(:,13),master2_bin(:,5), 1);
    ee_fit_plot=polyval(ee_fit,master2_bin(:,13));
%     ee_plot=plot(master2_bin(:,13),ee_fit_plot,'k--','linewidth',3);     %uncomment for plot
    ee_val=ee_fit(1);
    
    legend('Microcrack Energy','best fit line', 'location','best');
    disp(strcat('Event Rate: ',num2str(ee_val),' events/ 1% strain'));
    
    output_bins1=[output_bins1, er_val];
    %%
    %outfile - microcrack statistics
    mfwrite=master2_bin(:,[13,2,5]);
    fname=strcat('mf_bin_statistics_',num2str(n1),'_',num2str(n2),'.dat');
    dlmwrite(fname,mfwrite,'delimiter','\t');
    %%
    %outputs from stress-strain
    max_S1=max(master2_bin(:,9));
    disp(strcat('Max S1: ',num2str(max_S1)));  
    axial_strain=master2_bin(end,13);
    disp(strcat('Max Axial Strain: ',num2str(axial_strain)));
    output_bins1=[output_bins1, max_S1, axial_strain];
    %%
    %D-value unclumped
    r_vector=(linspace(0,max_X,100))';   %change number of bins for speed
    cr=zeros(length(r_vector),1);
    n_count=zeros(length(r_vector),1);
    
    [a1,~]=find(master1(:,15)>=n1);
    [a2,~]=find(master1(:,15)<=n2);
    
    a3=intersect(a1,a2);
    master1_bin=master1(a3,:);
    [p1_new, ~]=size(master1_bin);
    if p1_new==0 || p1_new==1
        d_value=0;
    elseif p1_new==2
        d_value=1;
    else
        
        if length(a3)==n_sum
            disp('EVENTS MATCH FOR D-VALUE');
        else
            disp('CHECK master1_bin ccalculation FOR D-VALUE');
        end
        
        for i=1:length(cr)
            for j=1:p1_new
                for k=1:p1_new
                    r_len1=sqrt((master1_bin(j,6)-master1_bin(k,6))^2+(master1_bin(j,7)-master1_bin(k,7))^2);
                    r_len2=sqrt((master1_bin(j,9)-master1_bin(k,9))^2+(master1_bin(j,10)-master1_bin(k,10))^2);
                    r_len3=sqrt((master1_bin(j,6)-master1_bin(k,6))^2+(master1_bin(j,10)-master1_bin(k,10))^2);
                    r_len4=sqrt((master1_bin(j,9)-master1_bin(k,9))^2+(master1_bin(j,7)-master1_bin(k,7))^2);
                    if r_len1 <= r_vector(i) && r_len2 <= r_vector(i) && r_len3 <= r_vector(i) && r_len4 <= r_vector(i) && j~=k
                        n_count(i)=n_count(i)+1;
                    end
                end
            end
            cr(i)=(2*n_count(i))/(p1_new*(p1_new-1));
        end
        
        log_cr=log10(cr);
        log_r=log10(r_vector);
        
        %write output file for D-value
        dwrite=[r_vector,log_cr];
        fd=strcat('d-value_',num2str(n1),'_',num2str(n2),'.dat');
        dlmwrite(fd,dwrite,'delimiter','\t');
        
%         fig10=figure;
%         semilogx(r_vector,log_cr,'LineWidth',3);
%         hold on
%         axis tight
        
        r_x_range=[r_vector(3), r_vector(52)];
        
        [r_locs1]=find(r_vector(:)>=r_x_range(1));
        [r_locs2]=find(r_vector(:)<=r_x_range(2));
        [r_locs]=intersect(r_locs1,r_locs2);
        r_bfit=log_r(r_locs);
%         r_bfit_plot=r_vector(r_locs);
        cr_bfit=log_cr(r_locs);
        
        r_P=polyfit(r_bfit,cr_bfit,1);
        r_plotP=polyval(r_P,r_bfit);
        %     log_r_plot=log10(r_plotP);
%         r_bfp=semilogx(r_bfit_plot,r_plotP,'r--','LineWidth',5);
        %     r_bfp=semilogx(r_bfit,r_plotP,'r--');
%         hold on
%         legend('Spatial Damage Correlation relationship','Best fit line','Location','SouthEast');
%         hold on
        %     legend('Magnitude-Frequency relationship','Best fit line');
        %     hold on
%         xlabel('radius'); ylabel('log(C(r))');
%         title('D-value plot');
        d_value=r_P(1);
    end
    disp(strcat('D-value: ',num2str(d_value)));
    output_bins1=[output_bins1, d_value];      
    %%
    %moment statistics
    [m1,~]=find(clump_events(:,3)>=n1);
    [m2,~]=find(clump_events(:,3)<=n2);
    m3=intersect(m1,m2);
    clump_events_bin=clump_events(m3,:);
    
    [m1_new, ~]=size(master1_bin);
    
    avg_moment=mean(clump_events_bin(:,9));
    if isempty(avg_moment)==1
        avg_moment=0;
    end
    max_moment=max(clump_events_bin(:,9));
    if isempty(max_moment)==1
        max_moment=0;
    end
    moment_differential=max(clump_events_bin(:,9))-min(clump_events_bin(:,9));
    if isempty(moment_differential)==1
        moment_differential=0;
    end   
    disp(strcat('Mean Moment: ',num2str(avg_moment)));
    disp(strcat('Max Moment: ',num2str(max_moment)));
    disp(strcat('Max Moment Difference: ',num2str(moment_differential)));
     
    output_bins1=[output_bins1, avg_moment, max_moment, moment_differential];
    
    cm_std=std(clump_events_bin(:,9));
    disp(strcat('Moment Standard Deviation: ',num2str(cm_std)));
    cm_var=var(clump_events_bin(:,9));
    disp(strcat('Moment Variance: ',num2str(cm_var)));
    cm_kur=kurtosis(clump_events_bin(:,9));
    disp(strcat('Moment Kurtosis: ',num2str(cm_kur)));
    
    output_bins1=[output_bins1, cm_std, cm_var, cm_kur];
    %%
    %moment outfile
    mmwrite=clump_events(:,[1,4,9]);
    fname2=strcat('moment_statistics_',num2str(n1),'_',num2str(n2),'.dat');
    dlmwrite(fname2,mmwrite,'delimiter','\t');
        %%
    %b-value
    if m1_new==0
        b_value_c=0;
    else
        AE_vector_c=clump_events_bin(:,9);
        max_AE_c=max(AE_vector_c);
        min_AE_c=min(AE_vector_c);
        bins_c=100;
        stepsize_c=(max_AE_c-min_AE_c)/bins_c; %100 is the number of bins
        N_vector_c=zeros(bins_c,1);
        A_vector_c=(linspace(min_AE_c,max_AE_c,bins_c))';
        threshold_c=min_AE_c;
        
        for i=1:bins_c
            [loc_AE_c]=find(AE_vector_c(:)>=threshold_c);
            N_vector_c(i)=log10(length(loc_AE_c));
            threshold_c=threshold_c+stepsize_c;
        end
        
        %write output file for b-value
        bwrite=[A_vector_c,N_vector_c];
        fnameb=strcat('b_value_',num2str(n1),'_',num2str(n2),'.dat');
        dlmwrite(fnameb,bwrite,'delimiter','\t');
        
%         fig14=figure;
%         plot(A_vector_c,N_vector_c,'o-')
%         xlabel('AE Magnitude'); ylabel('N(>Amax)');
%         title('b-value plot - CLUMPED');
%         hold on

        x_range_c=[A_vector_c(1), A_vector_c(end)];%Range for fit
        
        [AE_locs1_c]=find(A_vector_c(:)>=x_range_c(1));
        [AE_locs2_c]=find(A_vector_c(:)<=x_range_c(2));
        [AE_locs_c]=intersect(AE_locs1_c,AE_locs2_c);
        bfit_N_c=N_vector_c(AE_locs_c);
        bfit_A_c=A_vector_c(AE_locs_c);
        P_c=polyfit(bfit_A_c,bfit_N_c,1);
        plotP_c=polyval(P_c,bfit_A_c);
%         bfp_c=plot(bfit_A_c,plotP_c,'r--','LineWidth',3);
%         hold on
%         legend('Magnitude-Frequency relationship','Best fit line','Location','SouthWest');
%         hold on
        
        b_value_c=-P_c(1);

%         set(gcf,'Color','w');
%         saveas(gcf,'b_value_clumped','epsc');
%         export_fig fig14 f14_bvalue -q101 -painters -nocrop -pdf -jpeg -eps
    end
    disp(strcat('b value: ',num2str(b_value_c)));
    output_bins1=[output_bins1, b_value_c];           
    %%
    % D-value
    
    [p1_c,~]=size(clump_events_bin);
    r_vector_c=(linspace(0,max_X,100))';
    cr_c=zeros(length(r_vector_c),1);
    n_count_c=zeros(length(r_vector_c),1);
    
    if p1_c==0 || p1_c==1
        d_value_c=0;
    elseif p1_c==2
        d_value_c=1;
    else
        for i=1:length(cr_c)
            for j=1:p1_c
                for k=1:p1_c
                    r_len1_c=sqrt((clump_events_bin(j,5)-clump_events_bin(k,5))^2+(clump_events_bin(j,6)-clump_events_bin(k,6))^2);
                    if r_len1_c <= r_vector_c(i) && j~=k
                        n_count_c(i)=n_count_c(i)+1;
                    end
                end
            end
            cr_c(i)=(2*n_count_c(i))/(p1_c*(p1_c-1));
        end
        
        log_cr_c=log10(cr_c);
        log_r_c=log10(r_vector_c);

        %write output file
        dwrite=[r_vector_c,log_cr_c];
        fd2=strcat('d-value_clumped_',num2str(n1),'_',num2str(n2),'.dat');
        dlmwrite(fd2,dwrite,'delimiter','\t');
        
        % use log_cr and r_vector for plotting
        %use log_cr and log_r for regression
        
%         fig15=figure;
%         semilogx(r_vector_c,log_cr_c,'LineWidth',3);
%         hold on
%         axis tight
                
        r_x_range_c=[r_vector_c(3), r_vector_c(52)];
        
        [r_locs1_c]=find(r_vector_c(:)>=r_x_range_c(1));
        [r_locs2_c]=find(r_vector_c(:)<=r_x_range_c(2));
        [r_locs_c]=intersect(r_locs1_c,r_locs2_c);
        r_bfit_c=log_r_c(r_locs_c);
        r_bfit_plot_c=r_vector_c(r_locs_c);
        cr_bfit_c=log_cr_c(r_locs_c);
        
        r_P_c=polyfit(r_bfit_c,cr_bfit_c,1);
        r_plotP_c=polyval(r_P_c,r_bfit_c);
%         r_bfp_c=semilogx(r_bfit_plot_c,r_plotP_c,'r--','LineWidth',5);
%         hold on
%         legend('Spatial Damage Correlation relationship','Best fit line','Location','SouthEast');
%         hold on
%         xlabel('radius (m)'); ylabel('log(C(r))');
%         title('D-value plot - CLUMPED');
        
        d_value_c=r_P_c(1);
    end
    disp(strcat('D-value Coupled: ',num2str(d_value_c)));
    output_bins1=[output_bins1, d_value_c];
    output_bins1(isnan(output_bins1)==1)=0;

    foutname=strcat('output_',num2str(n1),'_',num2str(n2),'.dat');
    dlmwrite(foutname,output_bins1,'delimiter','\t');
    
    s1_bin=s1_buff1(n2+1);
    
    output_bins1=[output_bins1, s1_bin];
    
    output_bins2=[output_bins2; output_bins1];
end
dlmwrite('output_bins_master.dat',output_bins2,'delimiter','\t');
    
    
    
  