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

%%
im = imread('Star_inpaint_18.png');
im = rgb2gray(im);
im = im2double(im);
%s = imread('Star_inpaint_11_0.png');
%im = rgb2gray(im);
%s = im2double(im);
f = im;%2*im-1;
u = f;
%{
Bespotrebno koristenje 2 poziva funkcije, al ionako nije toliko presudno
kao neke for petlje nize. Takodjer maknut uduplani u=f.
%}
[ny, nx] = size(f);
N = nx*ny;

%imshow(u, 'InitialMagnification', 1000);

%%
%{
MATLAB ima multiline komentare, kao sto vidis sve sto je uglavljeno izmedju
%{ i %} je multiline komentar.
Dakle ovo ispod bi trebalo napraviti lambda matricu preko maski, naime
f<0.1 i f>0.9 su maske koje su 1 kad je uvjet ispunjen, a 0 drugdje, sto
upravo i zelimo.
%}
lmbda = L*(f<0.1 + f>0.9);

M=zeros(ny,nx);
for k=0:(nx-1)
    for l=0:(ny-1)
        M(l+1,k+1) = -2*(cos(2*pi*l/ny)-1)-2*(cos(2*pi*k/nx)-1); %
    end
end
M_2 = M.^2;

%% shock filter pocetni
%{
% Multiline komentar
[nably, nablydirection] = imgradient(u, 'central');
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
    fft2_u = fft2(u);
    ftu = (dt*(1/eps*fft2(pot1)+fft2(lmbda.*(f-u))-c1*M.*fft2_u+c2*fft2_u)+fft2_u)./(1+eps*dt*M.^2-c1t*M+c2t);
    u = real(ifft2(ftu));
end
%}
%% shock filter potential

for t = 0:dt:T1
    %imshow(u, 'InitialMagnification', 1000);
    %title(["t = " t]);
    %drawnow;
    %{
    Ovdje se moze koristiti druga vrsta gradijenta, koja je brza. Memorija
    me ne zanima, brzina mi je bitna. Zato sam ovaj dio promjenio. Naime,
    ovo odmah racuna absolutnu vrijednost gradijenta. Ovo 'central' samo
    znaci da koristim onu centralnu verziju za gradijent, mozes koristit i
    onu koja je derivacija prema naprijed ili natrag, kako zelis.
    Takodjer racunam fft2(u) samo jednom, ne tri puta kao sto si ti radila.
    Za kraj izbacujem M.^2 van petlje i zamjenjujem ga sa M_2, tako da je i
    to brze.
    %}
    [nably, nablydirection] = imgradient(u, 'central');
    lapl = del2(u);
    sfilt = -lapl./(abs(lapl)+delta).*nably;
    %sfilt = 4*workim.^3 - 6*workim.^2 + 2*workim; %DOUBLE WELL u^2*(u-1)^2
    fft2_u = fft2(u);
    ftu = (dt*(-1/eps*M.*fft2(sfilt)+fft2(lmbda.*(f-u))+c1*M.*fft2_u+c2*fft2_u)+fft2_u)./(1+eps*dt*M_2+c1t*M+c2t);
    u = real(ifft2(ftu));
end
imshow(u, 'InitialMagnification', 1000);
title('prvi korak');
drawnow;

%%
for t = T1:dt:T2
    %imshow(u, 'InitialMagnification', 1000);
    %title(["t = " t]);
    %drawnow;
    [nably, nablydirection] = imgradient(u, 'central');
    lapl = del2(u);
    sfilt = -lapl./(abs(lapl)+delta).*nably;
    %sfilt = 4*workim.^3 - 6*workim.^2 + 2*workim; %DOUBLE WELL u^2*(u-1)^2
    fft2_u = fft2(u);
    ftu = (dt*(-1/eps*M.*fft2(sfilt)+fft2(lmbda.*(f-u))+c1*M.*fft2_u+c2*fft2_u)+fft2_u)./(1+eps*dt*M_2+c1t*M+c2t);
    u = real(ifft2(ftu));
end
imshow(u, 'InitialMagnification', 1000);
title('drugi korak');
drawnow;

%%
%L2 = sum((s-u).^2,'all');
%disp(L2)
%{
Opet, zamjena koristenjem maski umjesto for petlji.
%}
u(u<0.5) = 0;
u(u>=0.5) = 1;
imshow(u, 'InitialMagnification', 1000);
title('treći korak');
drawnow;

f = u;
lmbda = L*(f<0.1 + f>0.9);

for t = T2:dt:T3
    %imshow(u, 'InitialMagnification', 1000);
    %title(["t = " t]);
    %drawnow;
    [nably, nablydirection] = imgradient(u, 'central');
    lapl = del2(u);
    sfilt = -lapl./(abs(lapl)+delta).*nably;
    %sfilt = 4*workim.^3 - 6*workim.^2 + 2*workim; %DOUBLE WELL u^2*(u-1)^2
    fft2_u = fft2(u);
    ftu = (dt*(-1/eps*M.*fft2(sfilt)+fft2(lmbda.*(f-u))+c1*M.*fft2_u+c2*fft2_u)+fft2_u)./(1+eps*dt*M_2+c1t*M+c2t);
    u = real(ifft2(ftu));
end
%imshow(u, 'InitialMagnification', 1000);
%title('treći korak');
%drawnow;
%{
% Ovo ne znam za sta sluzi pa cu ostaviti...
f = u;
lambda2 = L*(f<0.1 + f>0.9);
%}
