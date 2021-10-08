%%
T1 = 400;
T2 = 500;
T3 = T2;%500;%1500;
dt = 1;
eps = 70;
eps2 = 1.2;
eps3 = 0.6;
c1 = 1;
L = 50;
c2 = L;
delta = 0.0625;
et = eps*dt;
c1t = dt*c1;
c2t = dt*c2;


im = imread('Star_inpaint_18.png');
im = rgb2gray(im);
im = im2double(im);
%s = imread('Star_inpaint_11_0.png');
%im = rgb2gray(im);
%s = im2double(im);
f = im;%2*im-1;
u = f;
nx = size(f,2);
ny = size(f,1);
N = nx*ny;
u = f;

%imshow(u, 'InitialMagnification', 1000);

lmbda = L*ones(ny,nx);
for i = 1:nx
    for j=1:ny
        if f(j,i)<0.9 && f(j,i)>0.1%-0.9
            lmbda(j,i) = 0;
        end
    end
end

M=zeros(ny,nx);
for k=0:(nx-1)
    for l=0:(ny-1)
        M(l+1,k+1) = -2*(cos(2*pi*l/ny)-1)-2*(cos(2*pi*k/nx)-1); %
    end
end

%% shock filter pocetni
% % [ux,uy] = gradient(u);
% % nably = sqrt(ux.^2+uy.^2);
% % lapl = del2(u);
% % u = -sign(lapl).*nably;
% % 
% % figure;
% % subplot(2,1,1);
% % imshow((f+1)/2);
% % subplot(2,1,2);
% % for t = 0:dt:T1
% %     imshow((u+1)/2, 'InitialMagnification', 1000);
% %     title(["t = " t]);
% %     pot1 = u.^3-u;
% %     drawnow;
% %     ftu = (dt*(1/eps*fft2(u)+fft2(lmbda.*(f-u))-c1*M.*fft2(u)+c2*fft2(u))+fft2(u))./(1+eps*dt*M.^2-c1t*M+c2t);
% %     u = real(ifft2(ftu));
% % end

%% shock filter potential

for t = 0:dt:T1
    %imshow(u, 'InitialMagnification', 1000);
    %title(["t = " t]);
    %drawnow;
    [ux,uy] = gradient(u);
    nably = sqrt(ux.^2+uy.^2);
    lapl = del2(u);
    sfilt = -lapl./(abs(lapl)+delta).*nably;
    %sfilt = 4*workim.^3 - 6*workim.^2 + 2*workim; %DOUBLE WELL u^2*(u-1)^2
    ftu = (dt*(-1/eps*M.*fft2(sfilt)+fft2(lmbda.*(f-u))+c1*M.*fft2(u)+c2*fft2(u))+fft2(u))./(1+eps*dt*M.^2+c1t*M+c2t);
    u = real(ifft2(ftu));
end
imshow(u, 'InitialMagnification', 1000);
title('prvi korak');
drawnow;
for t = T1:dt:T2
    %imshow(u, 'InitialMagnification', 1000);
    %title(["t = " t]);
    %drawnow;
    [ux,uy] = gradient(u);
    nably = sqrt(ux.^2+uy.^2);
    lapl = del2(u);
    sfilt = -lapl./(abs(lapl)+delta).*nably;
    %sfilt = 4*workim.^3 - 6*workim.^2 + 2*workim; %DOUBLE WELL u^2*(u-1)^2
    ftu = (dt*(1/eps2*M.*fft2(sfilt)+fft2(lmbda.*(f-u))-c1*M.*fft2(u)+c2*fft2(u))+fft2(u))./(1+eps2*dt*M.^2-c1t*M+c2t);
    u = real(ifft2(ftu));
end
imshow(u, 'InitialMagnification', 1000);
title('drugi korak');
drawnow;
%L2 = sum((s-u).^2,'all');
%disp(L2)
for i = 1:nx
    for j=1:ny
        if u(j,i)<0.5
            u(j,i) = 0;
        end
        
        if u(j,i) >= 0.5
            u(j,i) = 1;
        end
    end
end
imshow(u, 'InitialMagnification', 1000);
title('treći korak');
drawnow;

f = u;
lmbda = L*ones(ny,nx);
for i = 1:nx
    for j=1:ny
        if f(j,i)<0.9 && f(j,i)>0.1
            lmbda(j,i) = 0;
            %u(j,i) = 0;
        end
    end
end


for t = T2:dt:T3
    %imshow(u, 'InitialMagnification', 1000);
    %title(["t = " t]);
    %drawnow;
    [ux,uy] = gradient(u);
    nably = sqrt(ux.^2+uy.^2);
    lapl = del2(u);
    sfilt = -lapl./(abs(lapl)+delta).*nably;
    %sfilt = 4*workim.^3 - 6*workim.^2 + 2*workim; %DOUBLE WELL u^2*(u-1)^2
    ftu = (dt*(1/eps3*M.*fft2(sfilt)+fft2(lmbda.*(f-u))-c1*M.*fft2(u)+c2*fft2(u))+fft2(u))./(1+eps3*dt*M.^2-c1t*M+c2t);
    u = real(ifft2(ftu));
end
%imshow(u, 'InitialMagnification', 1000);
%title('treći korak');
%drawnow;
%f = u;
%lambda2 = L*ones(ny,nx);
%for i = 1:nx
%    for j=1:ny
%        if f(j,i)<0.9 && f(j,i)>-0.9
%            lmbda2(j,i) = 0;
%            %u(j,i) = 0;
%        end
%    end
%end
