%%
T1 = 200;
dt = 1;
eps = 0.01;
c1 = 20;
L = 100000000;
c2 = L;
et = eps*dt;
c1t = dt*c1;
c2t = dt*c2;

im = imread('Star_inpaint_1.png');
im = rgb2gray(im);
im = im2double(im);
f = 2*im-1;
u = f;
nx = size(f,2);
ny = size(f,1);
N = nx*ny;
u = f;

lmbda = L*ones(ny,nx);
for i = 1:nx
    for j=1:ny
        if f(j,i)<0.9 && f(j,i)>-0.9 
            lmbda(j,i) = 0;
        end
    end
end

M=zeros(ny,nx);
for k=0:(nx-1)
    for l=0:(ny-1)
        M(l+1,k+1) = 2*ny^2*(cos(2*pi*l/ny)-1)+2*nx^2*(cos(2*pi*k/nx)-1); %
    end
end

%% shock filter pocetni
[ux,uy] = gradient(u);
nably = sqrt(ux.^2+uy.^2);
lapl = del2(u);
u = -sign(lapl).*nably;

figure;
subplot(2,1,1);
imshow((f+1)/2);
subplot(2,1,2);
for t = 0:dt:T1
    imshow((u+1)/2, 'InitialMagnification', 1000);
    title(["t = " t]);
    pot1 = u.^3-u;
    drawnow;
    ftu = (dt*(1/eps*fft2(u)+fft2(lmbda.*(f-u))-c1*M.*fft2(u)+c2*fft2(u))+fft2(u))./(1+eps*dt*M.^2-c1t*M+c2t);
    u = real(ifft2(ftu));
end

%% shock filter potential
for t = 0:dt:T1
    imshow((u+1)/2, 'InitialMagnification', 1000);
    title(["t = " t]);
    drawnow;
    [ux,uy] = gradient(u);
    nably = sqrt(ux.^2+uy.^2);
    lapl = del2(u);
    sfilt = -sign(lapl).*nably;
    ftu = (dt*(1/eps*M.*fft2(sfilt)+fft2(lmbda.*(f-u))-c1*M.*fft2(u)+c2*fft2(u))+fft2(u))./(1+eps*dt*M.^2-c1t*M+c2t);
    u = real(ifft2(ftu));
end
