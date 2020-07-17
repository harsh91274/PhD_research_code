function shale_bed_seeder

clc
clear all
close all

disp('Enter even values !')
sigma=input('Enter thickness of plates (4) : ');
eps=sigma/4;
t=sigma;

factor=input('Enter factor for length (10) : ');
w=factor*sigma;

[sh_init_seed]=makeseed1(sigma,eps,t,w);
sh_init_phi=porosity_calc(sh_init_seed);
disp(strcat('Porosity: ',num2str(sh_init_phi)));
[nx_init,ny_init,nz_init]=size(sh_init_seed);
makeseeder(sh_init_seed);
disp('sh_seed_initial');
disp(strcat('nx = ',num2str(nx_init),', ny = ',num2str(ny_init),', nz = ',num2str(nz_init)));
% imaging(sh_init_seed);
rotateseed(sh_init_seed);
disp('Rotated Seeds Created !');

again=1;

while again==1
    disp('Select option:');
    disp('1. 3D Dilation');
    disp('2. Microfractures parallel to bedding');
    disp('3. Microfractures perpendicular to bedding');
    disp('4. Fracture parallel to bedding');
    disp('5. Fracture perpendicular to bedding');
    choice=input('Enter choice: ');
    
    if choice==1
        mag=input('Enter Dilation Magnitude: ');
        eps1=eps+mag;
        sigma1=sigma+mag;
        [sh_pore_3d]=makeseed1(sigma1,eps1,t,w);
        sh_phi_pore_3d=porosity_calc(sh_pore_3d);
        disp(strcat('Porosity: ',num2str(sh_phi_pore_3d)));
        [nx_3d,ny_3d,nz_3d]=size(sh_pore_3d);
        disp('sh_seed_Z');
        disp(strcat('nx = ',num2str(nx_3d),', ny = ',num2str(ny_3d),', nz = ',num2str(nz_3d)));
        makeseeder1(sh_pore_3d);
        rotateseed(sh_pore_3d);
%         imaging(sh_pore_3d);
        disp('Rotated Seeds Created !');
     
    elseif choice==2
        mag=input('Enter magnitude of microfracture growth along bedding plane: ');
        eps1=eps+mag;
        [sh_pore_mf_parallel]=makeseed1(sigma,eps1,t,w);
        sh_phi_mf_parallel=porosity_calc(sh_pore_mf_parallel);
        disp(strcat('Porosity: ',num2str(sh_phi_mf_parallel)));
        [nx_mf_parallel,ny_mf_parallel,nz_mf_parallel]=size(sh_pore_mf_parallel);
        disp('sh_seed_Z');
        disp(strcat('nx = ',num2str(nx_mf_parallel),', ny = ',num2str(ny_mf_parallel),', nz = ',num2str(nz_mf_parallel)));
        makeseeder1(sh_pore_mf_parallel);
        rotateseed(sh_pore_mf_parallel);
%         imaging(sh_pore_mf_parallel);
        disp('Rotated Seeds Created !');
        
    elseif choice==3
        mag=input('Enter magnitude of microfracture growth perpendicular to bedding plane: ');
        sigma1=sigma+mag;
        [sh_pore_mf_perp]=makeseed1(sigma1,eps,t,w);
        sh_phi_mf_perp=porosity_calc(sh_pore_mf_perp);
        disp(strcat('Porosity: ',num2str(sh_phi_mf_perp)));
        [nx_mf_perp,ny_mf_perp,nz_mf_perp]=size(sh_pore_mf_perp);
        disp('sh_seed_Z');
        disp(strcat('nx = ',num2str(nx_mf_perp),', ny = ',num2str(ny_mf_perp),', nz = ',num2str(nz_mf_perp)));
        makeseeder1(sh_pore_mf_perp);
        rotateseed(sh_pore_mf_perp);
        %         imaging(sh_pore_mf_perp);
        disp('Rotated Seeds Created !');
        cut_choice=input('Cut seed (1=Yes)? ');
        if cut_choice==1
            cut_mf_seed=cutseed(sh_pore_mf_perp,t,eps);
            cut_rotateseed(cut_mf_seed);
            disp('Rotated Cut Seeds Created !');
        end
        
    elseif choice==4
        mag=input('Enter magnitude of fracture parallel to bedding: ');
        [sh_pore_fr_parallel]=fracture_parallel(eps, t,sh_init_seed,mag); 
        sh_phi_pore_fr_parallel=porosity_calc(sh_pore_fr_parallel);
        disp(strcat('Porosity: ',num2str(sh_phi_pore_fr_parallel)));
        [nx_fr_para,ny_fr_para,nz_fr_para]=size(sh_pore_fr_parallel);
        disp('sh_seed_Z');
        disp(strcat('nx = ',num2str(nx_fr_para),', ny = ',num2str(ny_fr_para),', nz = ',num2str(nz_fr_para)));
        makeseeder1(sh_pore_fr_parallel);
        rotateseed(sh_pore_fr_parallel);
%         imaging(sh_pore_fr_parallel);
        disp('Rotated Seeds Created !');
        
    elseif choice==5
        mag=input('Enter magnitude of fracture perpendicular to bedding: ');
        [sh_pore_fr_perp]=fracture_perp(sigma, eps, t, w, sh_init_seed,mag); 
        sh_phi_pore_fr_perp=porosity_calc(sh_pore_fr_perp);
        disp(strcat('Porosity: ',num2str(sh_phi_pore_fr_perp)));
        [nx_fr_perp,ny_fr_perp,nz_fr_perp]=size(sh_pore_fr_perp);
        disp('sh_seed_Z');
        disp(strcat('nx = ',num2str(nx_fr_perp),', ny = ',num2str(ny_fr_perp),', nz = ',num2str(nz_fr_perp)));
        makeseeder1(sh_pore_fr_perp);
        rotateseed(sh_pore_fr_perp);
%         imaging(sh_pore_fr_perp);
        disp('Rotated Seeds Created !');
        if cut_choice==1
            cut_f_seed=cutseed(sh_pore_fr_perp,t,eps);
            cut_rotateseed(cut_f_seed);
            disp('Rotated Cut Seeds Created !');
        end
    end    
    again=input('Again(1)? :');
end


end

function [cut_seed]= cutseed(seed,t,eps)
cut_seed=seed(:,:,t+eps+1:2*t+eps);
cut_makeseeder1(cut_seed);
cut_phi=porosity_calc(cut_seed);
[cut_nx,cut_ny,cut_nz]=size(cut_seed);
disp(strcat('Porosity: ',num2str(cut_phi)));
disp('cut_seed_Z.dat');
disp(strcat('nx = ',num2str(cut_nx),', ny = ',num2str(cut_ny),', nz = ',num2str(cut_nz)));
imaging(cut_seed);
end


function [frac_pore]=fracture_perp(sigma, eps, t, w, init_seed,mag)

[nx,ny,nz]=size(init_seed);
nx_new=nx+mag;

frac_pore=zeros(nx_new,ny,nz);

cube1=zeros(nx_new, ny, t);
cube2=zeros(nx_new, ny, t);
cube3=zeros(nx_new, ny, t);

cube1(1:w,:,:)=init_seed(1:w,:,1:t);
cube1(w+sigma+mag+1:end,:,:)=init_seed(w+sigma+1:end,:,1:t);

cube2(1:1.5*w+0.5*sigma,:,:)=init_seed(1:1.5*w+0.5*sigma,:,t+eps+1:2*t+eps);
cube2(1.5*w+1.5*sigma+mag+1:end,:,:)=init_seed(1.5*w+1.5*sigma+1:end,:,t+eps+1:2*t+eps);

cube3(1:2*w+sigma,:,:)=init_seed(1:2*w+sigma,:,2*t+2*eps+1:3*t+2*eps);
cube3(2*w+2*sigma+mag+1:end,:,:)=init_seed(2*w+2*sigma+1:end,:,2*t+2*eps+1:3*t+2*eps);

frac_pore(:,:,1:t)=cube1;
frac_pore(:,:,t+eps+1:2*t+eps)=cube2;
frac_pore(:,:,2*t+2*eps+1:end)=cube3;
end

function [frac_pore]=fracture_parallel(eps,t,init_seed,mag)
[nx,ny,nz]=size(init_seed);
nz_new=nz+mag;

frac_pore=zeros(nx,ny,nz_new);

frac_pore(:,:,1:2*t+2*eps)=init_seed(:,:,1:2*t+2*eps);
frac_pore(:,:,2*eps+2*t+mag+1:2*eps+3*t+mag)=init_seed(:,:,2*eps+2*t+1:2*eps+3*t);
end

function [seed]= makeseed1(sigma,eps, t,w)

% cubex=3*sigma+3*w;
% cubey=3*sigma+3*w;
% cubez=3*eps+3*t;
% % cubex=2*sigma+3*w;
% % cubey=2*sigma+3*w;
% % cubez=2*eps+3*t;

cubex=3*w+2*sigma;
cubey=3*w+2*sigma;
cubez=3*t+2*eps;

seed=zeros(cubex,cubey,cubez);

seed1=zeros(cubex,cubey,t);
seed2=zeros(cubex,cubey,t);

%seed1
for i=1:w
    for j=1:w
        seed1(i,j,:)=1;
    end
    for j=w+sigma+1:2*w+sigma
        seed1(i,j,:)=1;
    end
    for j=2*w+2*sigma+1:3*w+2*sigma
        seed1(i,j,:)=1;
    end
end

for i=w+sigma+1:2*w+sigma
    for j=1:w
        seed1(i,j,:)=1;
    end
    for j=w+sigma+1:2*w+sigma
        seed1(i,j,:)=1;
    end
    for j=2*w+2*sigma+1:3*w+2*sigma
        seed1(i,j,:)=1;
    end
end
for i=2*w+2*sigma+1:3*w+2*sigma
    for j=1:w
        seed1(i,j,:)=1;
    end
    for j=w+sigma+1:2*w+sigma
        seed1(i,j,:)=1;
    end
    for j=2*w+2*sigma+1:3*w+2*sigma
        seed1(i,j,:)=1;
    end
end

% for i=1:w
%     for j=1:w
%         seed1(i,j,:)=1;
%     end
%     for j=sigma+w+1:sigma+2*w
%         seed1(i,j,:)=1;
%     end
%     for j=2*sigma+2*w+1:2*sigma+3*w
%         seed1(i,j,:)=1;
%     end
% end
% 
% for i=sigma+w+1:sigma+2*w
%     for j=1:w
%         seed1(i,j,:)=1;
%     end
%     for j=sigma+w+1:sigma+2*w
%         seed1(i,j,:)=1;
%     end
%     for j=2*sigma+2*w+1:2*sigma+3*w
%         seed1(i,j,:)=1;
%     end
% end
% 
% for i=2*sigma+2*w+1:2*sigma+3*w
%     for j=1:w
%         seed1(i,j,:)=1;
%     end
%     for j=sigma+w+1:sigma+2*w
%         seed1(i,j,:)=1;
%     end
%     for j=2*sigma+2*w+1:2*sigma+3*w
%         seed1(i,j,:)=1;
%     end
% end
% imaging(seed1)

%seed 2

for i=1:(w-sigma)/2
    for j=1:(w-sigma)/2
        seed2(i,j,:)=1;
    end
    for j=(w+sigma)/2+1:1.5*w+0.5*sigma
        seed2(i,j,:)=1;
    end
    for j=1.5*w+1.5*sigma+1:2.5*w+1.5*sigma
        seed2(i,j,:)=1;
    end
    for j=2.5*w+2.5*sigma+1:3*w+2*sigma
        seed2(i,j,:)=1;
    end
end
for i=(w+sigma)/2+1:1.5*w+0.5*sigma
    for j=1:(w-sigma)/2
        seed2(i,j,:)=1;
    end
    for j=(w+sigma)/2+1:1.5*w+0.5*sigma
        seed2(i,j,:)=1;
    end
    for j=1.5*w+1.5*sigma+1:2.5*w+1.5*sigma
        seed2(i,j,:)=1;
    end
    for j=2.5*w+2.5*sigma+1:3*w+2*sigma
        seed2(i,j,:)=1;
    end
end
for i=1.5*w+1.5*sigma+1:2.5*w+1.5*sigma
    for j=1:(w-sigma)/2
        seed2(i,j,:)=1;
    end
    for j=(w+sigma)/2+1:1.5*w+0.5*sigma
        seed2(i,j,:)=1;
    end
    for j=1.5*w+1.5*sigma+1:2.5*w+1.5*sigma
        seed2(i,j,:)=1;
    end
    for j=2.5*w+2.5*sigma+1:3*w+2*sigma
        seed2(i,j,:)=1;
    end
end
for i=2.5*w+2.5*sigma+1:3*w+2*sigma
    for j=1:(w-sigma)/2
        seed2(i,j,:)=1;
    end
    for j=(w+sigma)/2+1:1.5*w+0.5*sigma
        seed2(i,j,:)=1;
    end
    for j=1.5*w+1.5*sigma+1:2.5*w+1.5*sigma
        seed2(i,j,:)=1;
    end
    for j=2.5*w+2.5*sigma+1:3*w+2*sigma
        seed2(i,j,:)=1;
    end
end
% for i=1:w/2 
%     for j=1:w/2
%         seed2(i,j,:)=1;
%     end
%     for j=w/2+sigma+1:1.5*w+sigma
%         seed2(i,j,:)=1;
%     end
%     for j=1.5*w+2*sigma+1:2.5*w+2*sigma
%         seed2(i,j,:)=1;
%     end
%     for j=2.5*w+3*sigma+1:3*w+3*sigma
%         seed2(i,j,:)=1;
%     end
% end
% 
% for i=w/2+sigma+1:1.5*w+sigma
%     for j=1:w/2
%         seed2(i,j,:)=1;
%     end
%     for j=w/2+sigma+1:1.5*w+sigma
%         seed2(i,j,:)=1;
%     end
%     for j=1.5*w+2*sigma+1:2.5*w+2*sigma
%         seed2(i,j,:)=1;
%     end
%     for j=2.5*w+3*sigma+1:3*w+3*sigma
%         seed2(i,j,:)=1;
%     end
% end
% 
% for i=1.5*w+2*sigma+1:2.5*w+2*sigma
%     for j=1:w/2
%         seed2(i,j,:)=1;
%     end
%     for j=w/2+sigma+1:1.5*w+sigma
%         seed2(i,j,:)=1;
%     end
%     for j=1.5*w+2*sigma+1:2.5*w+2*sigma
%         seed2(i,j,:)=1;
%     end
%     for j=2.5*w+3*sigma+1:3*w+3*sigma
%         seed2(i,j,:)=1;
%     end
% end
% 
% for i=2.5*w+3*sigma+1:3*w+3*sigma
%     for j=1:w/2
%         seed2(i,j,:)=1;
%     end
%     for j=w/2+sigma+1:1.5*w+sigma
%         seed2(i,j,:)=1;
%     end
%     for j=1.5*w+2*sigma+1:2.5*w+2*sigma
%         seed2(i,j,:)=1;
%     end
%     for j=2.5*w+3*sigma+1:3*w+3*sigma
%         seed2(i,j,:)=1;
%     end
% end
%     
% imaging(seed2)

%seed

% seed(:,:,eps/2+1:eps/2+t)=seed1;
% seed(:,:,1.5*eps+t+1:1.5*eps+2*t)=seed2;
% seed(:,:,2.5*eps+2*t+1:2.5*eps+3*t)=seed1;  

% seed(:,:,1:t)=seed1;
% seed(:,:,eps+t+1:eps+2*t)=seed2;
% seed(:,:,2*eps+2*t+1:2*eps+3*t)=seed1;  

seed(:,:,1:t)=seed1;
seed(:,:,t+eps+1:2*t+eps)=seed2;
seed(:,:,2*eps+2*t+1:3*t+2*eps)=seed1;

% imaging(seed)

end

function imaging(seed)
figure
p=patch(isosurface(seed,0.5),'FaceColor','blue');
isonormals(seed,p)
view(3)
axis vis3d image
xlabel('X');ylabel('Y');zlabel('Z');
camlight; lighting phong
end

function phi=porosity_calc(seed)
[nx,ny,nz]=size(seed);
phi=1-(sum(sum(sum(seed)))/(nx*ny*nz)); 
end

function makeseeder(seed)
[nx,ny,nz]=size(seed);
fid = fopen('sh_seed_initial.dat', 'w');
for i=1:nx
    for j=1:ny
        for k=1:nz
            fprintf(fid, '%i\n', seed(i,j,k));
        end
    end
end
fclose(fid);
end

function makeseeder1(seed)
[nx,ny,nz]=size(seed);
fid = fopen('sh_seed_Z.dat', 'w');
for i=1:nx
    for j=1:ny
        for k=1:nz
            fprintf(fid, '%i\n', seed(i,j,k));
        end
    end
end
fclose(fid);
end

function cut_makeseeder1(seed)
[nx,ny,nz]=size(seed);
fid = fopen('cut_seed_Z.dat', 'w');
for i=1:nx
    for j=1:ny
        for k=1:nz
            fprintf(fid, '%i\n', seed(i,j,k));
        end
    end
end
fclose(fid);
end

function rotateseed(seed)
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
    
disp('sh_seed_Y');
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

fid=fopen('sh_seed_Y.dat','w');
for i=1:nx1
    for j=1:ny1
        for k=1:nz1
            fprintf(fid, '%i\n', pseed1(i,j,k));
        end
    end
end

fclose(fid);

[nx2,ny2,nz2]=size(pseed2);
disp('sh_seed_X');
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
fid=fopen('sh_seed_X.dat','w');
for i=1:nx2
    for j=1:ny2
        for k=1:nz2
            fprintf(fid, '%i\n', pseed2(i,j,k));
        end
    end
end

fclose(fid);
end

function cut_rotateseed(seed)
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
    
disp('cut_seed_Y');
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

fid=fopen('cut_seed_Y.dat','w');
for i=1:nx1
    for j=1:ny1
        for k=1:nz1
            fprintf(fid, '%i\n', pseed1(i,j,k));
        end
    end
end

fclose(fid);

[nx2,ny2,nz2]=size(pseed2);
disp('cut_seed_X');
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
fid=fopen('cut_seed_X.dat','w');
for i=1:nx2
    for j=1:ny2
        for k=1:nz2
            fprintf(fid, '%i\n', pseed2(i,j,k));
        end
    end
end

fclose(fid);
end
