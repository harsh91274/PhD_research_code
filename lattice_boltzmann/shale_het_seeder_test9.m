function shale_het_seeder 
clc
clear all
close all

path(path,genpath('C:\matlab_adds'));
%modify lengths in make_het_seed too
% smec_t=1; 
% smec_l=50; 
% 
% ch29; 
% ill_l=1000; 
% 
% kao_t=75; 
% kao_l=1500; 
% 
% chlo_t=20;            %chlorite thickness
% chlo_l=500;              %chlorite platelet length

smec_t=1; 
smec_l=5; 

ill_t=5; 
ill_l=30; 

kao_t=20; 
kao_l=50; 

chlo_t=7;            %chlorite thickness
chlo_l=25;              %chlorite platelet length

vol_smec=smec_l*smec_l*smec_t;
vol_ill=ill_l*ill_l*ill_t;
vol_kao=kao_l*kao_l*kao_t;
vol_chlo=chlo_l*chlo_l*chlo_t;
%%
sigma=input('Enter Intrabed Pore throat Diameter (2) : ');
eps=sigma/2;

disp('There are 10 platelets in the seed');
disp('1=Smectite; 2=Illite; 3=Chlorite; 4=Kaolinite'); 
disp ('PLATELET LENGTHS IN MAKESEEDER'); 
disp('MODIFY INITIAL SEED ARRAY IN CODE');
% row1=input('Enter First Seeder Row (3 platelets): ');

%permutations of platelets in input seed
layer_order0=[1 1 2 1 3 1];  %approximate to mineralogy fractions in core
layer_order00=[1 1 2 1 3 1];  %approximate to mineralogy fractions in core
layer_order000=[1 1 2 1 3 1];  %approximate to mineralogy fractions in core

order_sat=0;
while order_sat==0
X1=randperm(numel(layer_order0));
layer_order1=reshape(layer_order0(X1),size(layer_order0));   
X2=randperm(numel(layer_order00));
layer_order2=reshape(layer_order00(X2),size(layer_order00));
X3=randperm(numel(layer_order000));
layer_order3=reshape(layer_order000(X3),size(layer_order000));
seed_order=[layer_order1; layer_order2; layer_order3]
order_sat=input('Satisfied (1=Yes, 2=No): ');
end
dlmwrite('platelet_order.dat',seed_order,'\t');

no1=length(find(seed_order==1));
no2=length(find(seed_order==2));
no3=length(find(seed_order==3));
no4=length(find(seed_order==4));
%%
%make initial seed
[t, sh_init_seed, no_babies]=make_het_seed2(sigma,eps, seed_order);   %initial seed constructed
[nx_init,ny_init,nz_init]=size(sh_init_seed);           %size of seed

no1=no1+no_babies;
vol_t1=no1*vol_smec;
vol_t2=no2*vol_ill;
vol_t3=no3*vol_chlo;
vol_t4=no4*vol_kao;
vol_grains=vol_t1+vol_t2+vol_t3+vol_t4;

volf1=vol_t1/vol_grains; 
disp(strcat('Smectite Volume Fraction: ',num2str(volf1)));
volf2=vol_t2/vol_grains; 
disp(strcat('Illite Volume Fraction: ',num2str(volf2)));
volf3=vol_t3/vol_grains; 
disp(strcat('Chlorite Volume Fraction: ',num2str(volf3)));
volf4=vol_t4/vol_grains; 
disp(strcat('Kaolinite Volume Fraction: ',num2str(volf4)));
sh_init_phi=porosity_calc(sh_init_seed);    %porosity calculation    
disp(strcat('Porosity: ',num2str(sh_init_phi)));
keyboard
imaging(sh_init_seed);
printseeder1(sh_init_seed,sh_init_phi); %write Z file
fname1=strcat('sh_het_Z',num2str(round(sh_init_phi*100)),'.dat');
disp(fname1);
disp(strcat('nx = ',num2str(nx_init),', ny = ',num2str(ny_init),', nz = ',num2str(nz_init)));   %display seed size
rotateseed(sh_init_seed,sh_init_phi);
disp('Rotated Seeds Created !');

or_option=input('Reorient Seed (1=Yes): ');
if or_option==1
    eps_mat=zeros(nx_init,ny_init,eps); 
    orient=input('Enter orientation (degrees): ');
    seed_orient(sh_init_seed, orient, eps_mat, t, eps, sh_init_phi);
end

%%
%dilation 
again=1;

while again==1
    disp('3D Dilation');
    mag=input('Enter Dilation Magnitude: ');
    eps1=eps+mag;
    sigma1=sigma+mag;
    [t, sh_pore_3d, ~]=make_het_seed2(sigma1,eps1,seed_order);
    sh_phi_pore_3d=porosity_calc(sh_pore_3d);
    disp(strcat('Porosity: ',num2str(sh_phi_pore_3d)));
    [nx_3d,ny_3d,nz_3d]=size(sh_pore_3d);
    fname1=strcat('sh_seed_Z',num2str(round(sh_phi_pore_3d*100)),'.dat');
    disp(fname1);
    disp(strcat('nx = ',num2str(nx_3d),', ny = ',num2str(ny_3d),', nz = ',num2str(nz_3d)));
    printseeder1(sh_pore_3d, sh_phi_pore_3d);
    rotateseed(sh_pore_3d, sh_phi_pore_3d);
%     imaging(sh_pore_3d);
    disp('Rotated Seeds Created !');
    
    or_option1=input('Reorient Seed (1=Yes): ');
    if or_option1==1
        eps1_mat=zeros(nx_3d,ny_3d,eps1);
        orient=input('Enter orientation (degrees): ');
        seed_orient(sh_pore_3d, orient, eps1_mat, t, eps1,sh_phi_pore_3d);
    end

	again=input('Again(1)? :'); 

end
end
%%
function [t1, seed, no_babies]= make_het_seed2(sigma,eps, seed_order)
smec_t=1; 
smec_l=5; 

ill_t=5; 
ill_l=30; 

kao_t=5; 
kao_l=50; 

chlo_t=7;            %chlorite thickness
chlo_l=25;              %chlorite platelet length

% smec_t=1; 
% smec_l=50; 
% 
% ill_t=29; 
% ill_l=1000; 
% 
% kao_t=75; 
% kao_l=1500; 
% 
% chlo_t=20;            %chlorite thickness
% chlo_l=500;              %chlorite platelet length
%%

seed_order1=seed_order(1,:);
seed_order2=seed_order(2,:);
seed_order3=seed_order(3,:);
%%

%determine size of each bed
thick_ind=max(seed_order(1,:));

if thick_ind==1
    seed_bed_thickness=smec_t;
elseif thick_ind==2
    seed_bed_thickness=ill_t;
elseif thick_ind==3
    seed_bed_thickness=chlo_t;
else
    seed_bed_thickness=kao_t;
end

%determine length of each bed
seed_bed_length=0;
[~,num_plates]=size(seed_order);

for i=1:num_plates
    if seed_order1(i)==1
        seed_bed_length=seed_bed_length+smec_l;
    elseif seed_order1(i)==2
        seed_bed_length=seed_bed_length+ill_l;
    elseif seed_order1(i)==3
        seed_bed_length=seed_bed_length+chlo_l;
    else
        seed_bed_length=seed_bed_length+kao_l;
    end
end

seed_bed_length=seed_bed_length+((num_plates-1)*sigma);
t1=seed_bed_thickness;

bed1=zeros(seed_bed_length,seed_bed_length,seed_bed_thickness);
bed2=zeros(seed_bed_length,seed_bed_length,seed_bed_thickness);
bed3=zeros(seed_bed_length,seed_bed_length,seed_bed_thickness);

ipore1=zeros(seed_bed_length,seed_bed_length,eps, 1);
ipore2=zeros(seed_bed_length,seed_bed_length,eps, 1);
%%
%find proxy
smallest_grain=min(seed_order(1,:));
seq=findseq(seed_order1);
ind_seq= seq(:,1)==smallest_grain;
ind_seq=find(max(seq(:,4)));
seql=seq(ind_seq,4);
seq_start=seq(ind_seq,2);
% disp(strcat('Smallest Grain is: ',num2str(smallest_grain)));
% disp(strcat('Longest Sequence of filler grain: ',num2str(seql),' platelets'));
% if isempty(seql)~=1
%     [babyseed, baby_l, no_baby_beds1]=makefiller(smallest_grain, seql, seed_bed_thickness, sigma, eps);
%     no_gbabies1=seql*no_baby_beds1;
% else
[babyseed, baby_l, no_baby_beds1]=makefiller(smallest_grain, 1, seed_bed_thickness,sigma, eps);
% no_gbabies1=1*no_baby_beds1;
[babx1, baby1, babyz1]=size(babyseed);

repx=2*ceil(seed_bed_length/babx1);
repy=2*ceil(seed_bed_length/baby1);
bed1_proxy=repmat(babyseed,[repx repy 1]);
% imaging(bed1_proxy); 

repx2=ceil(seed_bed_length/babx1);
repy2=ceil(seed_bed_length/baby1);

no_babies_proxy1=no_baby_beds1*repx2*repy2;
baby_correction1=babycalc(seed_order1, sigma)*no_baby_beds1;
no_babies1=no_babies_proxy1-baby_correction1;

fc=0;
baby_ind=0;
for i=1:num_plates
    if seed_order1(i)==1
        fc=fc+smec_l;
    elseif seed_order1(i)==2
        fc=fc+ill_l;
    elseif seed_order1(i)==3
        fc=fc+chlo_l;
    elseif seed_order1(i)==4
        fc=fc+kao_l;       
    end
    if i~=num_plates
        fc=fc+sigma;
    end
    if i==seq_start
        baby_ind=fc+1;
    end
end
 
ind_adjust=rem(baby_ind,(baby_l+sigma));
bed1(:,:,1:babyz1)=bed1_proxy(ind_adjust+1:seed_bed_length+ind_adjust, ind_adjust+1:seed_bed_length+ind_adjust, 1:babyz1);
%%
fc1=0;  %fill counter
for i=1:num_plates
    if seed_order1(i)==1
        bed1(fc1+1: fc1+smec_l, fc1+1: fc1+smec_l, 1:smec_t)=1;
        fc1=fc1+smec_l;
    elseif seed_order1(i)==2
        bed1(fc1+1: fc1+ill_l, fc1+1: fc1+ill_l, 1:ill_t)=1;
        fc1=fc1+ill_l;
    elseif seed_order1(i)==3
        bed1(fc1+1: fc1+chlo_l, fc1+1: fc1+chlo_l, 1:chlo_t)=1;
        fc1=fc1+chlo_l;
    elseif seed_order1(i)==4
        bed1(fc1+1: fc1+kao_l, fc1+1: fc1+kao_l, 1:kao_t)=1;
        fc1=fc1+kao_l;       
    end
    if i~=num_plates
        fc1=fc1+sigma;
    end
end
% imaging(bed1)
%%
%find proxy
smallest_grain=min(seed_order(2,:));
seq=findseq(seed_order2);
% ind_seq= seq(:,1)==smallest_grain;
ind_seq=find(max(seq(:,4)));
seql=seq(ind_seq,4);
seq_start=seq(ind_seq,2);
% disp(strcat('Smallest Grain is: ',num2str(smallest_grain)));
% disp(strcat('Longest Sequence of filler grain: ',num2str(seql),' platelets'));
% if isempty(seql)~=1
%     [babyseed, baby_l, no_baby_beds2]=makefiller(smallest_grain, seql, seed_bed_thickness, sigma, eps);
%     no_gbabies2=seql*no_baby_beds2;
% else
    [babyseed, baby_l, no_baby_beds2]=makefiller(smallest_grain, 1, seed_bed_thickness, sigma, eps);
%     no_gbabies2=1*no_baby_beds2;
% end
[babx2, baby2, babyz2]=size(babyseed);

repx=2*ceil(seed_bed_length/babx2);
repy=2*ceil(seed_bed_length/baby2);
bed2_proxy=repmat(babyseed,[repx repy 1]);

repx3=ceil(seed_bed_length/babx2);
repy3=ceil(seed_bed_length/baby2);

no_babies_proxy2=no_baby_beds2*repx3*repy3;
baby_correction2=babycalc(seed_order2, sigma)*no_baby_beds2;
no_babies2=no_babies_proxy2-baby_correction2;

fc=0;
baby_ind=0;
for i=1:num_plates
    if seed_order2(i)==1
        fc=fc+smec_l;
    elseif seed_order2(i)==2
        fc=fc+ill_l;
    elseif seed_order2(i)==3
        fc=fc+chlo_l;
    elseif seed_order2(i)==4
        fc=fc+kao_l;       
    end
    if i~=num_plates
        fc=fc+sigma;
    end
    if i==seq_start
        baby_ind=fc+1;
    end
end

ind_adjust=rem(baby_ind,(baby_l+sigma));
bed2(:,:,1:babyz2)=bed2_proxy(ind_adjust+1:seed_bed_length+ind_adjust, ind_adjust+1:seed_bed_length+ind_adjust, 1:babyz2);

% (bed2);
%%
fc2=0;  %fill counter
for i=1:num_plates
    if seed_order2(i)==1
        bed2(fc2+1: fc2+smec_l, fc2+1: fc2+smec_l, 1:smec_t)=1;
        fc2=fc2+smec_l;
    elseif seed_order2(i)==2
        bed2(fc2+1: fc2+ill_l, fc2+1: fc2+ill_l, 1:ill_t)=1;
        fc2=fc2+ill_l;
    elseif seed_order2(i)==3
        bed2(fc2+1: fc2+chlo_l, fc2+1: fc2+chlo_l, 1:chlo_t)=1;
        fc2=fc2+chlo_l;
    elseif seed_order2(i)==4
        bed2(fc2+1: fc2+kao_l, fc2+1: fc2+kao_l, 1:kao_t)=1;
        fc2=fc2+kao_l;       
    end
    if i~=num_plates
        fc2=fc2+sigma;
    end
end
% imaging(bed2)
%%
%find proxy
smallest_grain=min(seed_order(3,:));
seq=findseq(seed_order3);
% ind_seq= seq(:,1)==smallest_grain;
ind_seq=find(max(seq(:,4)));
seql=seq(ind_seq,4);
seq_start=seq(ind_seq,2);
% disp(strcat('Smallest Grain is: ',num2str(smallest_grain)));
% disp(strcat('Longest Sequence of filler grain: ',num2str(seql),' platelets'));
% if isempty(seql)~=1
%     [babyseed, baby_l, no_baby_beds3]=makefiller(smallest_grain, seql, seed_bed_thickness, sigma, eps);
%     no_gbabies3=seql*no_baby_beds3;
% else
    [babyseed, baby_l, no_baby_beds3]=makefiller(smallest_grain, 1, seed_bed_thickness, sigma, eps);
%     no_gbabies3=1*no_baby_beds3;
% end
[babx3, baby3, babyz3]=size(babyseed);

repx=2*ceil(seed_bed_length/babx3);
repy=2*ceil(seed_bed_length/baby3);

bed3_proxy=repmat(babyseed,[repx repy 1]);

repx4=ceil(seed_bed_length/babx3);
repy4=ceil(seed_bed_length/baby3);

no_babies_proxy3=no_baby_beds3*repx4*repy4;
baby_correction3=babycalc(seed_order3, sigma)*no_baby_beds3;
no_babies3=no_babies_proxy3-baby_correction3;

fc=0;
baby_ind=0;
for i=1:num_plates
    if seed_order3(i)==1
        fc=fc+smec_l;
    elseif seed_order3(i)==2
        fc=fc+ill_l;
    elseif seed_order3(i)==3
        fc=fc+chlo_l;
    elseif seed_order3(i)==4
        fc=fc+kao_l;       
    end
    if i~=num_plates
        fc=fc+sigma;
    end
    if i==seq_start
        baby_ind=fc+1;
    end
end
 
ind_adjust=rem(baby_ind,(baby_l+sigma));
bed3(:,:,1:babyz3)=bed3_proxy(ind_adjust+1:seed_bed_length+ind_adjust, ind_adjust+1:seed_bed_length+ind_adjust, 1:babyz3);
% imaging(bed3);
%%
fc3=0;  %fill counter
for i=1:num_plates
    if seed_order3(i)==1
        bed3(fc3+1: fc3+smec_l, fc3+1: fc3+smec_l, 1:smec_t)=1;
        fc3=fc3+smec_l;
    elseif seed_order3(i)==2
        bed3(fc3+1: fc3+ill_l, fc3+1: fc3+ill_l, 1:ill_t)=1;
        fc3=fc3+ill_l;
    elseif seed_order3(i)==3
        bed3(fc3+1: fc3+chlo_l, fc3+1: fc3+chlo_l, 1:chlo_t)=1;
        fc3=fc3+chlo_l;
    elseif seed_order3(i)==4
        bed3(fc3+1: fc3+kao_l, fc3+1: fc3+kao_l, 1:kao_t)=1;
        fc3=fc3+kao_l;       
    end
    if i~=num_plates
        fc3=fc3+sigma;
    end
end

% imaging(bed3)
seed=cat(3, bed1, ipore1, bed2, ipore2, bed3);

% babies=[no_babies1 no_babies2 no_babies3];
no_babies=no_babies1+no_babies2+no_babies3;
end
%%
function [babyseed, baby_l, no_baby_beds]=makefiller(gtype,len,t,sigma, eps)
% smec_t=1; 
% smec_l=50; 
% 
% ill_t=29; 
% ill_l=1000; 
% 
% kao_t=75; 
% kao_l=1500; 
% 
% chlo_t=20;            %chlorite thickness
% chlo_l=500; 

smec_t=1; 
smec_l=5; 

ill_t=5; 
ill_l=30; 

kao_t=5; 
kao_l=50; 

chlo_t=7;            %chlorite thickness
chlo_l=25;              %chlorite platelet length

%%
if gtype==1
    baby_l=smec_l; baby_t=smec_t;
elseif gtype==2
    baby_l=ill_l; baby_t=ill_t;
elseif gtype==3
    baby_l=chlo_l; baby_t=chlo_t;
else 
    baby_l=kao_l; baby_t=kao_t;
end
cubex=baby_l + sigma;
cubey=baby_l + sigma;
cubez=t; 

bed0=zeros(cubex, cubey, baby_t+eps);
babyseed=zeros(cubex, cubey, cubez);
no_baby_beds=floor(t/(baby_t+eps));

bed0(1:baby_l,1:baby_l,1:baby_t)=1;
bed1=bed0;
if no_baby_beds>1
    bed_proxy=bed0;
    for i=2:no_baby_beds;
        bed2=circshift(bed_proxy,[sigma sigma 0]);
        bed1=cat(3, bed1,bed2);
        bed_proxy=bed2;
    end
end
% bed1=repmat(bed0, [1 1 no_baby_beds]);
% 

[~,~,nz]=size(bed1);
babyseed(:,:,1:nz)=bed1;

end

function baby_correction=babycalc(bed_order, sigma)   %works only for smectite
% smec_t=1; 
% smec_l=50; 
% 
% ill_t=29; 
% ill_l=1000; 
% 
% kao_t=75; 
% kao_l=1500; make
% 
% chlo_t=20;            %chlorite thickness
% chlo_l=500; 

% smec_t=1; 
smec_l=5; 

% ill_t=5; 
ill_l=30; 

% kao_t=5; 
kao_l=50; 

% chlo_t=7;            %chlorite thickness
chlo_l=25;              %chlorite platelet length

no_ill=length(find(bed_order==2));
no_chlo=length(find(bed_order==3));
no_kao=length(find(bed_order==4));
% baby_correction=round(no_ill*((ill_l^2 *ill_t)/((smec_l+sigma)^2)*(smec_t+eps)) + no_chlo*((chlo_l^2*chlo_t)/((smec_l+sigma)^2)*(smec_t+eps)) + no_kao*((kao_l^2*kao_t)/((smec_l+sigma)^2)*(smec_t+eps)));  
baby_corr1=no_ill*((ill_l/(smec_l+sigma))*(ill_l/(smec_l+sigma)));
baby_corr2=no_chlo*((chlo_l/(smec_l+sigma))*(chlo_l/(smec_l+sigma)));
baby_corr3=no_kao*((kao_l/(smec_l+sigma))*(kao_l/(smec_l+sigma)));
baby_correction=round(baby_corr1+baby_corr2+baby_corr3);
end

%%
function imaging(seed)
figure
% spy(seed);
p = patch(isosurface(seed,0.5));
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
xlabel('x'); ylabel('y'); zlabel('z');
daspect([1 1 1]);
view(3); 
axis tight;
camlight; lighting gouraud
end

function phi=porosity_calc(seed)
[nx,ny,nz]=size(seed);
phi=1-(sum(sum(sum(seed)))/(nx*ny*nz)); 
end

function printseeder1(seed, phi)
fname=strcat('sh_het_Z',num2str(round(phi*100)),'.dat');
fid = fopen(fname, 'w');
fprintf(fid, '%i\n', seed);
fclose(fid);
end

function rotateseed(seed, phi)
[nx,ny,nz]=size(seed);
    pseed1=zeros(nx,nz,ny);   
    for i=1:nx
        for j=1:ny
            for k=1:nz
                pseed1(i,k,j)=seed(i,j,k);
            end
        end
    end    
   pseed2=zeros(nz,ny,nx);    
    for i=1:nx
        for j=1:ny
            for k=1:nz
                pseed2(k,j,i)=seed(i,j,k);
            end
        end
    end
 
fname2=strcat('sh_het_Y',num2str(round(phi*100)),'.dat');
disp(fname2);


[nx1,ny1,nz1]=size(pseed1);
disp(strcat('nx = ',num2str(nx1),', ny = ',num2str(ny1),', nz = ',num2str(nz1)));

fid=fopen(fname2,'w');
fprintf(fid, '%i\n', pseed1);
fclose(fid);

imaging(pseed1);

[nx2,ny2,nz2]=size(pseed2);
fname3=strcat('sh_het_X',num2str(round(phi*100)),'.dat');
disp(fname3);

disp(strcat('nx = ',num2str(nx2),', ny = ',num2str(ny2),', nz = ',num2str(nz2)));

fid=fopen(fname3,'w');
fprintf(fid, '%i\n', pseed2);

fclose(fid);
imaging(pseed2);
end
%%
function seed_orient(blob_old,orient, eps_mat, t, eps, phi)
[slice1, slice2]=cutseed2(blob_old, t, eps);
blob=cat(3, slice1, eps_mat, blob_old, eps_mat, slice2, eps_mat);
% imaging(blob);

theta=orient*(pi/180);

blob_center = (size(blob) + 1) / 2;
T1 = [1 0 0 0 ; 0 1 0 0; 0 0 1 0; -blob_center 1];
T2 = [cos(theta) 0 -sin(theta) 0; 0 1 0 0; sin(theta) 0  cos(theta)   0; 0 0 0 1];
T3 = [1 0 0 0; 0 1 0 0; 0 0 1 0; blob_center 1];
T = T1 * T2 * T3;
tform = maketform('affine', T);
R = makeresampler('nearest', 'circular');
TDIMS_A = [1 2 3];
TDIMS_B = [1 2 3];
TSIZE_B = size(blob_old);
TMAP_B = [];
F = 0;
blob2 = tformarray(blob, tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);

imaging(blob2);
check=unique(blob2);
if length(check)==2
    disp('Transformed values are 0 and 1');
else 
    disp('Check transformed values for numbers besides 0 and 1');
end

por_old=1-((sum(sum(sum(blob_old))))/(size(blob_old,1)*size(blob_old,2)*size(blob_old,3)));
por1=1-((sum(sum(sum(blob))))/(size(blob,1)*size(blob,2)*size(blob,3)));
por2=1-((sum(sum(sum(blob2))))/(size(blob2,1)*size(blob2,2)*size(blob2,3)));

disp(strcat('old seed porosity: ',num2str(por_old)));
disp(strcat('transformation porosity: ',num2str(por1)));
disp(strcat('final seed porosity: ',num2str(por2)));

orientseeder(blob2,orient, phi, por2);

end

function [cut_seed1, cut_seed2]= cutseed2(seed,t,eps)
[nx,ny,nz]=size(seed);
    pseed1=zeros(nx,ny,nz);
    pseed2=zeros(nx,ny,nz); 
    for i=1:nx
        for j=1:ny
            for k=1:nz
                pseed1(i,j,k)=seed(nx-i+1,j,k);
                pseed2(i,j,k)=seed(i,ny-j+1,k);
            end
        end
    end    
cut_seed1=pseed1(:,:,t+eps+1:2*t+eps);
cut_seed2=pseed2(:,:,t+eps+1:2*t+eps);
% cut_makeseeder1(cut_seed);
% cut_phi=porosity_calc(cut_seed);
% [cut_nx,cut_ny,cut_nz]=size(cut_seed);
% disp(strcat('Porosity: ',num2str(cut_phi)));
% disp('cut_seed_Z.dat');
% disp(strcat('nx = ',num2str(cut_nx),', ny = ',num2str(cut_ny),', nz = ',num2str(cut_nz)));
%  imaging(cut_seed);
end


function orientseeder(seed,orient, phi, phi2)
[nx,ny,nz]=size(seed);
fnamo=(strcat('sh_seed_Z',num2str(round(100*phi)),'-',num2str(round(100*phi2)),'_orient',num2str(orient),'.dat'));
disp(fnamo);
disp(strcat('nx = ',num2str(nx),', ny = ',num2str(ny),', nz = ',num2str(nz)));
fid = fopen(fnamo, 'w');
fprintf(fid, '%i\n', seed);
fclose(fid);
orient_rotateseed(seed,orient,phi,phi2);
end

function orient_rotateseed(seed,orient, phi, phi2)
[nx,ny,nz]=size(seed);
    pseed1=zeros(nx,nz,ny);   
    for i=1:nx
        for j=1:ny
            for k=1:nz
                pseed1(i,k,j)=seed(i,j,k);
            end
        end
    end    
   pseed2=zeros(nz,ny,nx);    
    for i=1:nx
        for j=1:ny
            for k=1:nz
                pseed2(k,j,i)=seed(i,j,k);
            end
        end
    end

% fname2=(strcat('sh_seed_Y_orient',num2str(orient),'.dat'));    
% disp(fname2);
fnamo2=(strcat('sh_seed_Y',num2str(round(100*phi)),'-',num2str(round(100*phi2)),'_orient',num2str(orient),'.dat'));
disp(fnamo2);
[nx1,ny1,nz1]=size(pseed1);
disp(strcat('nx = ',num2str(nx1),', ny = ',num2str(ny1),', nz = ',num2str(nz1)));

% if choice==2
%     fid=fopen('sh_seed_frac_Y.dat','w');
%     disp('sh_seed_frac_Y');
%     disp(strcat('nx = ',num2str(nx1),', ny = ',num2str(ny1),', nz = ',num2str(nz1)));
% elseif choice==1
%     fid=fopen('pl_seed_3d_Y.dat','w');
%     disp('sh_seed_3d_Y');
%     disp(strcat('nx = ',num2str(nx1),', ny = ',num2str(ny1),', nz = ',num2str(nz1)));
% end

fid=fopen(fnamo2,'w');
fprintf(fid, '%i\n', pseed1);
fclose(fid);

% figure
% imaging(pseed1);

[nx2,ny2,nz2]=size(pseed2);
fnamo3=(strcat('sh_seed_X',num2str(round(100*phi)),'-',num2str(round(100*phi2)),'_orient',num2str(orient),'.dat'));
disp(fnamo3);
% fname3=(strcat('sh_seed_X_orient',num2str(orient),'.dat'));
% disp(fname3);
disp(strcat('nx = ',num2str(nx2),', ny = ',num2str(ny2),', nz = ',num2str(nz2)));

% if choice==2
%     fid=fopen('sh_seed_frac_X.dat','w');
%     disp('sh_seed_frac_X');
%     disp(strcat('nx = ',num2str(nx2),', ny = ',num2str(ny2),', nz = ',num2str(nz2)));
% elseif choice==1
%     fid=fopen('sh_seed_3d_X.dat','w');
%     disp('sh_seed_3d_X');
%     disp(strcat('nx = ',num2str(nx2),', ny = ',num2str(ny2),', nz = ',num2str(nz2)));
% end

fid=fopen(fnamo3,'w');
fprintf(fid, '%i\n', pseed2);
fclose(fid);

% imaging(pseed2);
end
