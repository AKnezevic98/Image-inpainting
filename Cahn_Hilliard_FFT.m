filespath = [pwd '\'];
imagename = 'Star_inpaint_5.png';

im = imread([filespath imagename]);
im = rgb2gray(im);
im = im2double(im);
[n,m] = size(im);
N = n*m;

%%
T = 1000;
swtch_t = 100;
epsilon = 1;
swtch_eps = 0.001;
l0 = 10000;
C1 = 3333;
C2 = 3*l0;
dt = 1;

workim = im(:);
lambda = l0*ones([N,1]);
for i=1:N
    if workim(i) < 0.8 && workim(i) > 0.2
        lambda(i) = 0;
    end
end
lambda = reshape(lambda, [n,m]);
workim = reshape(workim, [n,m]);

K2 = zeros([n,m]);
for i=1:n
    for j=1:m
        K2(i,j) = 4*pi^2 *(((i)/n)^2 + ((j)/m)^2);
    end
end
K4 = K2.^2;

%%
figure;
subplot(2,1,1);
imshow(im);
subplot(2,1,2);
for i=0:dt:T
    if i>swtch_t
        epsilon = swtch_eps;
    end
    workim_tilda = fft2(workim);
    F = 4*workim.^3 - 6*workim.^2 + 2*workim;
    F_tilda = fft2(F);
    lambda_u = lambda.*(im-workim);
    lambda_u_tilda = fft2(lambda_u);
    rhs = workim_tilda - (dt/epsilon)*K2.*F_tilda + dt*(lambda_u_tilda + C1*K2.*workim_tilda + C2*workim_tilda);
    lhs = 1*ones([n,m]) + dt*(epsilon*K4 + C1*K2 + C2*ones([n,m]));
    workim_tilda = rhs./lhs;
    workim = real(ifft2(workim_tilda));
    imshow(workim);
    title(["t = " num2str(i*dt)]);
    drawnow;
end
