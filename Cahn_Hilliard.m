filespath = [pwd '\'];
imagename = 'Star_inpaint.png';

im = imread([filespath imagename]);
im = rgb2gray(im);
im = im2double(im);
im = 2*im-1;
[n,m] = size(im);
N = n*m;

%%
epsilon = 0.8;
l0 = 0;
C1 = 300;
C2 = 3*l0;
dt = 1;

workim = im(:);

lambda = l0*ones([N,1]);
for i=1:N
    if workim(i) < 0.8 && workim(i) > -0.8
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

%{
for i=1:N
    for j=1:N
        if i==j+1 || j==i+1 || i==j+(n-1) || j==i+(n-1) || i==j+n+1 || j==i+n+1
            if mod(i,n)==0 || mod(j,n)==0
                M(i,j) = 0;
            end
        elseif i==j+2 || j==i+2
            if mod(i,n)==0 || mod(j,n)==0 || mod(i,n)==n-1 || mod(j,n)==n-1
                M(i,j) = 0;
            end
        end
    end
end
%}

d0 = -4*dt*ones([N,1]);
d1 = 1*dt*ones([N-1,1]);
d2 = 1*dt*ones([N-n,1]);
Lapl = diag(d0) + diag(d1, 1) + diag(d1, -1) + diag(d2, n) + diag(d2, -n);

%%

figure;
subplot(2,1,1);
imshow((im+1)/2);
subplot(2,1,2);
for i=1:2400
    Wim = workim.^3 - workim; 
    rhs = Lapl*(Wim/epsilon - C1*workim) + (1+C2*dt)*workim + dt*lambda.*(im(:) - workim);
    workim = M\rhs;
    workim = reshape(workim, [n,m]);
    imshow((workim+1)/2);
    title(["t = " num2str(i*dt)]);
    drawnow;
    workim = workim(:);
end
