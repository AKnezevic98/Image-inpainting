filespath = [pwd '\'];
imagename = 'Star_inpaint.png';

im = imread([filespath imagename]);
im = rgb2gray(im);
im = im2double(im);
im = 2*im-1;
[n,m] = size(im);
N = n*m;

%%
epsilon = 0.1;
C1 = 10;
C2 = 1;
dt = 0.001;
l0 = 1;

workim = im(:);

lambda = l0*ones([N,1]);
for i=1:N
    if workim(i) ~= 1 || workim(i)~=-1
        lambda(i) = 0;
    end
end

d0 = (1 + 20*epsilon*dt + 4*C1*dt + C2*dt)*ones([N,1]);
dy1 = (-8*epsilon*dt - C1*dt)*ones([N-1,1]);
dy2 = 1*epsilon*dt*ones([N-2,1]);
dx1 = (-8*epsilon*dt - C1*dt)*ones([N-n,1]);
dx2 = 1*epsilon*dt*ones([N-2*n,1]);
dx1y1 = 2*epsilon*dt*ones([N-1-n, 1]);
dy1x1 = 2*epsilon*dt*ones([N-n+1, 1]);


M = diag(d0);
M = M + diag(dy1, 1) + diag(dy1, -1);
M = M + diag(dy2, 2) + diag(dy2, -2);
M = M + diag(dx1, n) + diag(dx1, -n);
M = M + diag(dx2, 2*n) + diag(dx2, -2*n);
M = M + diag(dx1y1, n+1) + diag(dx1y1, -n-1);
M = M + diag(dy1x1, n-1) + diag(dy1x1, -n+1);

d0 = -4*dt/epsilon*ones([N,1]);
d1 = 1*dt/epsilon*ones([N-1,1]);
d2 = 1*dt/epsilon*ones([N-n,1]);
Lapl = diag(d0) + diag(d1, 1) + diag(d1, -1) + diag(d2, n) + diag(d2, -n);

%%

figure;
subplot(2,1,1);
imshow((im+1)/2);
subplot(2,1,2);
for i=1:10
    Wim = workim.^3 - workim; 
    rhs = Lapl*Wim + workim - C1*epsilon*Lapl*workim + C2*workim + lambda.*(im(:) - workim);
    workim = M\rhs;
    workim = reshape(workim, [n,m]);
    imshow((workim+1)/2);
    drawnow;
    workim = workim(:);
end
