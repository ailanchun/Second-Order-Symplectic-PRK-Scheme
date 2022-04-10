clc;
clear all;
close all;
dx=0.0025;
dz=0.0025;
x=0:dx:1;
z=-0.2:dx:0.8;
z1=-0.2:dz:0;
z2=-0.2:dz:0.2;
z3=-0.2:dz:0.4;
nx=length(x);
nz=length(z);
z1nz=length(z1);
z2nz=length(z2);
z3nz=length(z3);

epp(1:nx,1:nz)=1;
ep(1:nx,1:nz)=1;
mu(1:nx,1:nz)=1;
chiyu(1:nx,1:nz)=0;
sig(1:nx,1:nz)=0;

epp(1:nx,z1nz:nz)=6;
ep(1:nx,z1nz:nz)=6;
mu(1:nx,z1nz:nz)=1;
chiyu(1:nx,z1nz:nz)=0;
sig(1:nx,z1nz:nz)=0.001;

figure; 
subplot(2,2,1);
imagesc(x,z,ep'); axis image; colorbar
xlabel('x (m)'); ylabel('z (m)');
title('Original \ep matrix');
subplot(2,2,2)
imagesc(x,z,epp'); axis image; colorbar
xlabel('x (m)'); ylabel('z (m)');
title('Interpolated \epp matrix');
subplot(2,2,3)
imagesc(x,z,chiyu'); axis image; colorbar
xlabel('x (m)'); ylabel('z (m)');
title('Interpolated \chiyu');

disp(['Using dx = ',num2str(dx),' m, dz = ',num2str(dz),' m']);
dt =5e-12;
disp(['Using dt = ',num2str(dt/1e-9),' ns']);
disp(' ');

t=0:dt:12e-9;                                   
srcpulse = ricker(dt,1000e6,12e-9);

srcx=0.5; 
srcz=0*ones(size(srcx));
recx =srcx+0.02;
recz =srcz;
srcloc = [srcx srcz];
recloc = [recx recz];

% set some output and plotting parameters
outstep = 1;
plotopt = [1 250 0.002];
disp('Press any key to begin simulation...');
disp(' ');
pause;
close all

% run the simulation
tic;
xprop=x;
zprop=z;

% determine true permittivity and permeability matrices from supplied relative ones
ep0 = 8.8541878176e-12;% dielectric permittivity of free space
mu0 = 1.2566370614e-6;% magnetic permeability of free space
ep=ep*ep0;
epp=epp*ep0;
mu = mu*mu0;% true permeability matrix

nx = length(xprop);
nz = length(zprop) ;         % maximum number of field nodes in the z-direction

dx = (xprop(2)-xprop(1));        % electric and magnetic field spatial discretization in x (m)
dz = (zprop(2)-zprop(1));         % electric and magnetic field spatial discretization in z (m)

% x and z position vectors corresponding to Hx, Hz, and Ey field matrices
% (these field matrices are staggered in both space and time, and thus have different coordinates)
xHx = xprop(1):dx:xprop(end);
zHx = zprop(1):dz:zprop(end);
xHz = xprop(1):dx:xprop(end);
zHz = zprop(1):dz:zprop(end);
xEy = xHx;
zEy = zHz;

% determine source and receiver (i,j) indices in Ey field matrix,
% and true coordinates of sources and receivers in numerical model (after discretization)
nsrc = size(srcloc,1);                        % number of sources
nrec = size(recloc,1);                          % number of receivers
for s=1:nsrc                                  
    [temp,srci(s)] = min(abs(xEy-srcloc(s,1))); % source x index in Ey field matrix
    [temp,srcj(s)] = min(abs(zEy-srcloc(s,2))); % source z index in Ey field matrix
    srcx(s) = xEy(srci(s));                     % true source x location
    srcz(s) = zEy(srcj(s));                     % true source z location
end
for r=1:nrec;                                   
    [temp,reci(r)] = min(abs(xEy-recloc(r,1))); % receiver x index in Ey field matrix
    [temp,recj(r)] = min(abs(zEy-recloc(r,2))); % receiver z index in Ey field matrix
    recx(r) = xEy(reci(r));                    % true receiver x location
    recz(r) = zEy(recj(r));                  % true receiver z location
end

% determine time stepping parameters from supplied time vector
dt = t(2)-t(1) ;                                % temporal discretization
numit = length(t);                              % number of iterations

c = 2.996e8;
m = 1:nx;
n = 1:nz;

	 vp=zeros(nx,nz);
     s=zeros(nx,nz);
     t1=zeros(nx,nz);
     t2=zeros(nx,nz);
     t3=zeros(nx,nz);
     t21=zeros(nx,nz);
     t22=zeros(nx,nz);
     t23=zeros(nx,nz);
     t24=zeros(nx,nz);
     t25=zeros(nx,nz);
  
vp(m,n) = c ./ (sqrt(ep(m,n)./ep0)); 

  s(m,n)=vp(m,n).*dt/dx;  
  t1(m,n)=(s(m,n)-1).*(s(m,n)-2)./2;
  t2(m,n)=(-s(m,n)).*(s(m,n)-2);
  t3(m,n)=(s(m,n)-1).*s(m,n)./2 ;
  t21(m,n)=t1(m,n).*t1(m,n);
  t22(m,n)=2.*t1(m,n).*t2(m,n);
  t23(m,n)=2.*t1(m,n).*t3(m,n)+t2(m,n).*t2(m,n);
  t24(m,n)=2.*t2(m,n).*t3(m,n);
  t25(m,n)=t3(m,n).*t3(m,n) ;
	 
     
  cp=zeros(nx,nz);%CP(m)
  cq=zeros(nx,nz);%CQ(m)
  cp(m,n)=1;
  cq(m,n)= dt./ (mu(m,n).*dx) ;
  c1=zeros(nx,nz);
  c2=zeros(nx,nz);
  c3=zeros(nx,nz);
  c4=zeros(nx,nz);
  c5=zeros(nx,nz);
  c6=zeros(nx,nz);
  c1(m,n)=sig(m,n)*0.5+(sig(m,n).*chiyu(m,n)+ep(m,n))./dt+epp(m,n).*chiyu(m,n)./(dt*dt);
  c2(m,n)=sig(m,n)*0.5-(sig(m,n).*chiyu(m,n)+ep(m,n))./dt-2*epp(m,n).*chiyu(m,n)./(dt*dt);
 
  c3(m,n)=epp(m,n).*chiyu(m,n)./(dt*dt);
  c4(m,n)=chiyu(m,n)/(dt*dt)+(1/dt);
  c5(m,n)=-2*chiyu(m,n)/(dt*dt)-(1/dt);
  c6(m,n)=chiyu(m,n)/(dt*dt);

disp('Beginning FDTD simulation...')


% initialize gather matrix where data will be stored
gather = zeros(fix((numit-1)/outstep)+1,nrec,nsrc);
% loop over number of sources
for s=1:nsrc;

    % zero all field matrices
    Ey = zeros(nx,nz);          % Ey component of electric field
       Ey1 = zeros(nx,nz);          % Ey component of electric field
       Ey2 = zeros(nx,nz);          % Ey component of electric field
      Dy = zeros(nx,nz); 
     
      Dy_old=zeros(nx,nz);
      Dy_old_old=zeros(nx,nz);
      Hx = zeros(nx,nz);            % Hx component of magnetic field
     Hy = zeros(nx,nz); 
    % Hz component of magnetic field

    % time stepping loop
    for it=1:numit; 
     i = 1:nx-1;  j = 1:nz-1;                % indices for all components in Hx matrix to update
                   
       Hx(i,j) =cp(i,j).*Hx(i,j) - cq(i,j).*(Ey(i,j+1)-Ey(i,j));       
       Hy(i,j) = cp(i,j).*Hy(i,j) + cq(i,j).*(Ey(i+1,j)-Ey(i,j));
           
        
        i = 2:nx;  j = 2:nz; % indices for all components in Ey matrix to update update to be applied to the whole Ey grid
         Dy_old_old(1:end,1:end)=Dy_old(1:end,1:end); 
         Dy_old(1:end,1:end)=Dy(1:end,1:end); 
       Dy(i,j) = Dy_old(i,j)+(dt/dx)*((Hy(i,j)-Hy(i-1,j))-(Hx(i,j)-Hx(i,j-1)));
  
       
 Ey2(1:end,1:end)=Ey1(1:end,1:end) ;
        Ey1(1:end,1:end)=Ey(1:end,1:end)  ;           
                 
         
%   Ey(i,j)=CU(i,j)*Ey(i,j)+CV(i,j)*Dy(i,j)-CW(i,j)*Dy_old(i,j)
Ey(i,j)=(1./c1(i,j)).*(-c2(i,j).*Ey1(i,j)-c3(i,j).*Ey2(i,j)+c4(i,j).*Dy(i,j)+c5(i,j).*Dy_old(i,j)+c6(i,j).*Dy_old_old(i,j));
     
        
        % indices for all components in Ey matrix to update
        % add source pulse to Ey at source location
        % (emulates infinitesimal Ey directed line source with current = srcpulse)
       
        x = srci(s);
        y = srcj(s);
        Ey(x,y) = Ey(x,y) + srcpulse(it); 

        
    
        Ey(1,j)=2*(t1(1,j).*Ey1(1,j) + t2(1,j).*Ey1(2,j) + t3(1,j).*Ey1(3,j))- ...
     (t21(1,j).*Ey2(1,j) + t22(1,j).*Ey2(2,j) + t23(1,j).*Ey2(3,j) + t24(1,j).* ...
     Ey2(4,j) + t25(1,j).*Ey2(5,j));
       
%right
     
        Ey(nx,j)=2.*(t1(nx,j).*Ey1(nx,j) + t2(nx,j).*Ey1(nx-1,j) + t3(nx,j).* ...
     Ey1(nx-2,j)) - (t21(nx,j).*Ey2(nx,j) + t22(nx,j).*Ey2(nx-1,j) + t23(nx,j).* ...
     Ey2(nx-2,j) + t24(nx,j).*Ey2(nx-3,j) + t25(nx,j).*Ey2(nx-4,j));

%!bottom
        
        Ey(i,1)=2.*(t1(i,1).*Ey1(i,1) + t2(i,1).*Ey1(i,2) + t3(i,1).*Ey1(i,3))- ...
     (t21(i,1).*Ey2(i,1) + t22(i,1).*Ey2(i,2) + t23(i,1).*Ey2(i,3) + t24(i,1).* ...
     Ey2(i,4) + t25(i,1).*Ey2(i,5));
   
%!top
    
       Ey(i,nz)=2.*(t1(i,nz).*Ey1(i,nz) + t2(i,nz).*Ey1(i,nz-1) + t3(i,nz).* ...
     Ey1(i,nz-2))-(t21(i,nz).*Ey2(i,nz) + t22(i,nz).*Ey2(i,nz-1) + t23(i,nz).* ...
     Ey2(i,nz-2) + t24(i,nz).*Ey2(i,nz-3) + t25(i,nz).*Ey2(i,nz-4));
  
%!corner

       Ey(1,1)=2*(t1(1,1)*Ey1(1,1)+t2(1,1)*Ey1(2,2)+t3(1,1)*Ey1(3,3))- ...
     (t21(1,1)*Ey2(1,1)+t22(1,1)*Ey2(2,2)+t23(1,1)*Ey2(3,3)+t24(1,1) ...
     *Ey2(4,4)+t25(1,1)*Ey2(5,5));


       Ey(1,nz)=2*(t1(1,nz)*Ey1(1,nz)+t2(1,nz)*Ey1(2,nz-1)+t3(1,nz)* ...
     Ey1(3,nz-2))-(t21(1,nz)*Ey2(1,nz)+t22(1,nz)*Ey2(2,nz-1)+t23(1,nz) ...
     *Ey2(3,nz-2)+t24(1,nz)*Ey2(4,nz-3)+t25(1,nz)*Ey2(5,nz-4));

       Ey(nx,1)=2*(t1(nx,1)*Ey1(nx,1)+t2(nx,1)*Ey1(nx-1,2)+t3(nx,1)* ...
     Ey1(nx-2,3))-(t21(nx,1)*Ey2(nx,1)+t22(nx,1)*Ey2(nx-1,2)+t23(nx,1) ...
     *Ey2(nx-2,3)+t24(nx,1)*Ey2(nx-3,4)+t25(nx,1)*Ey2(nx-4,5));

       Ey(nx,nz)=2*(t1(nx,nz)*Ey1(nx,nz)+t2(nx,nz)*Ey1(nx-1,nz-1)+ ...
     t3(nx,nz)*Ey1(nx-2,nz-2))-(t21(nx,nz)*Ey2(nx,nz)+t22(nx,nz)* ...
     Ey2(nx-1,nz-1)+t23(nx,nz)*Ey2(nx-2,nz-2)+t24(nx,nz)*Ey2(nx-3,nz-3) ...
     +t25(nx,nz)*Ey2(nx-4,nz-4));	 
	 
        % record the results in gather matrix if necessary
        if mod(it-1,outstep)==0
            tout((it-1)/outstep+1) = t(it);
            for r=1:nrec
                gather((it-1)/outstep+1,r,s) = Ey(reci(r),recj(r));
            end
        end
  
    end
end

xlabel('Position (m)');
ylabel('Time (ns)');
set(get(gca,'XLabel'),'FontSize',15);
set(get(gca,'YLabel'),'FontSize',15);
set(gca,'linewidth',1.5);
set(gca,'FontSize',15);
fid=fopen('outfile.txt','w');
fprintf(fid,'%g\n',tout*1e9,co_data)
fclose(fid)

