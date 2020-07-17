function shale_het_seeder 
clc
clear all
close all

%modify lengths in make_het_seed too
% smec_t=1; 
% smec_l=50; 
% 
% ill_t=29; 
% ill_l=1000; 
% 
% kao_t=75; 
% kao_l=1500; 
% 
% chlo_t=50;            %chlorite thickness
% chlo_l=1200;              %chlorite platelet length

% smec_t=2; 
% smec_l=100; 
% 
% ill_t=58; 
% ill_l=2000; 
% 
% kao_t=150; 
% kao_l=3000; 
% 
% chlo_t=50;            %chlorite thickness
% chlo_l=1200;              %chlorite platelet length

smec_t=2; 
smec_l=10; 

ill_t=6; 
ill_l=200; 

kao_t=15; 
kao_l=300; 

chlo_t=10;            %chlorite thickness
chlo_l=250;  
%%
sigma=input('Enter Intrabed Pore throat Diameter (2) : ');
eps=sigma/2;

disp('There are 10 platelets in the seed');
disp('1=Smectite; 2=Illite; 3=Chlorite; 4=Kaolinite'); 
disp ('PLATELET LENGTHS IN MAKESEEDER'); 
disp('MODIFY INITIAL SEED ARRAY IN CODE');
% row1=input('Enter First Seeder Row (3 platelets): ');

%permutations of platelets in input seed
layer_order0=[1 1 2 1 1];  %approximate to mineralogy fractions in core

order_sat=0;
while order_sat==0
X1=randperm(numel(layer_order0));
layer_order1=reshape(layer_order0(X1),size(layer_order0));   
X2=randperm(numel(layer_order1));
layer_order2=reshape(layer_order1(X2),size(layer_order1));
X3=randperm(numel(layer_order2));
layer_order3=reshape(layer_order2(X3),size(layer_order2));
seed_order=[layer_order1; layer_order2; layer_order3]
order_sat=input('Satisfied (1=Yes, 2=No): ');
end
dlmwrite('platelet_order.dat',seed_order,'\t');
%%
%make initial seed
[t, sh_init_seed]=make_het_seed2(sigma,eps, seed_order);   %initial seed constructed
[nx_init,ny_init,nz_init]=size(sh_init_seed);           %size of seed
imaging(sh_init_seed);
sh_init_phi=porosity_calc(sh_init_seed);    %porosity calculation    
printseeder1(sh_init_seed,sh_init_phi); %write Z file
fname1=strcat('sh_het_Z',num2str(round(sh_init_phi*100)),'.dat');     
disp(strcat('Porosity: ',num2str(sh_init_phi)));
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
    [t, sh_pore_3d]=make_het_seed2(sigma1,eps1,seed_order);
    sh_phi_pore_3d=porosity_calc(sh_pore_3d);
    disp(strcat('Porosity: ',num2str(sh_phi_pore_3d)));
    [nx_3d,ny_3d,nz_3d]=size(sh_pore_3d);
    fname1=strcat('sh_seed_Z',num2str(round(sh_phi_pore_3d*100)),'.dat');
    disp(fname1);
    disp(strcat('nx = ',num2str(nx_3d),', ny = ',num2str(ny_3d),', nz = ',num2str(nz_3d)));
    printseeder1(sh_pore_3d, sh_phi_pore_3d);
    rotateseed(sh_pore_3d, sh_phi_pore_3d);
    imaging(sh_pore_3d);
    disp('Rotated Seeds Created !');

	again=input('Again(1)? :'); 

end
end
%%
function [t1, seed]= make_het_seed2(sigma,eps, seed_order)
% smec_t=2; 
% smec_l=100; 
% 
% ill_t=58; 
% ill_l=2000; 
% 
% kao_t=150; 
% kao_l=3000; 
% 
% chlo_t=50;            %chlorite thickness
% chlo_l=1200;           %chlorite platelet length

% smec_t=1; 
% smec_l=50; 
% 
% ill_t=29; 
% ill_l=1000; 
% 
% kao_t=75; 
% kao_l=1500; 
% 
% chlo_t=50;            %chlorite thickness
% chlo_l=1200;              %chlorite platelet length


smec_t=2; 
smec_l=10; 

ill_t=6; 
ill_l=200; 

kao_t=15; 
kao_l=300; 

chlo_t=10;            %chlorite thickness
chlo_l=250;  
%%

seed_order1=seed_order(1,:);
seed_order2=seed_order(2,:);
seed_order3=seed_order(3,:);

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
% 
% phi_bed=porosity_calc(bed1);
% disp(strcat('Bed Porosity: ',num2str(phi_bed)));

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

seed=cat(3, bed1, ipore1, bed2, ipore2, bed3);
end

function imaging(seed)
figure
% spy(seed);
p = patch(isosurface(seed,0.5));
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
daspect([1 1 1]);
view(3); axis tight;
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

[nx2,ny2,nz2]=size(pseed2);
fname3=strcat('sh_het_X',num2str(round(phi*100)),'.dat');
disp(fname3);

disp(strcat('nx = ',num2str(nx2),', ny = ',num2str(ny2),', nz = ',num2str(nz2)));

fid=fopen(fname3,'w');
fprintf(fid, '%i\n', pseed2);

fclose(fid);
end

function seed_orient(blob_old,orient, eps_mat, t, eps, phi)
slice=cutseed(blob_old, t, eps);
blob=cat(3, slice, eps_mat, blob_old, eps_mat, slice, eps_mat);
imaging(blob);

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

function [cut_seed]= cutseed(seed,t,eps)
[nx,ny,nz]=size(seed);
    pseed0=zeros(ny,nx,nz);   
    for i=1:nx
        for j=1:ny
            for k=1:nz
                pseed0(j,i,k)=seed(i,j,k);
            end
        end
    end    
cut_seed=pseed0(:,:,t+eps+1:2*t+eps);
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
imaging(pseed1);

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

imaging(pseed2);
end
