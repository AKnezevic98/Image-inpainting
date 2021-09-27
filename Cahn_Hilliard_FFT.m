%% IMAGE IMPORTING
filespath = [pwd '\'];
imagename = 'Star_inpaint.png';
im = imread([filespath imagename]);
im = rgb2gray(im);
im = im2double(im);
[n,m] = size(im);
N = n*m;

%% PARAMETERS
T = 4000;
swtch_t = T/5;
epsilon = 1;
delta = 0.0625;
swtch_eps = 0.01;
l0 = 500;
C1 = 100;
C2 = 1*l0;
dt = 1;

%% REDEFINING IN CASE OF LOG POTENTIAL
%{
im = im(:);
for i=1:N
    if im(i)==1
        im(i)=0.999;
    elseif im(i)==0
        im(i)=0.001;
    end
end
im = reshape(im, [n,m]);
%}

%% FIDELITY TERM MASK
workim = im(:);
lambda = l0*ones([N,1]);
for i=1:N
    if workim(i) < 0.8 && workim(i) > 0.2
        lambda(i) = 0;
    end
end
lambda = reshape(lambda, [n,m]);
workim = reshape(workim, [n,m]);

%% K^2 AND K^4 MATRICES
K2 = zeros([n,m]);
for i=1:n
    for j=1:m
        K2(i,j) = 2*(2 - cos(2*pi*(i-1)/n) - cos(2*pi*(j-1)/m));
    end
end
K4 = K2.^2;

%% MASKS FOR SHOCK FILTER

%AVERAGED DERIVATIVE
gradx = [0,0,0; -1,0,1; 0,0,0]/2;
grady = [0,1,0; 0,0,0; 0,-1,0]/2;

%FORWARD DERIVATIVE
gradpx = [0,0,0; 0,-1,1; 0,0,0];
gradpy = [0,1,0; 0,-1,0; 0,0,0];

%BACKWARD DERIVATIVE
gradmx = [0,0,0; -1,1,0; 0,0,0];
gradmy = [0,0,0; 0,1,0; 0,-1,0];

%LAPLACE
laplace = [0,1,0; 1,-4,1; 0,1,0];

%% SHOCK FILTER ON INITIAL IMAGE
%{
gradxim = imfilter(im, gradx);
gradyim = imfilter(im, grady);
agim = sqrt(gradxim.^2 + gradyim.^2);

laplaceim = imfilter(im, laplace);
slim = sign(laplaceim);
F = -agim.*slim;

F_tilda = fft2(F);
%}

%% ALGORITHM
swtch_t = floor(swtch_t);

%FIGURE DEFINING
figure;
subplot(1,3,1);
imshow(im);
subplot(1,3,3);

for i=0:dt:T
    
    %SWITCHING AFTER CERTAIN t
    if i==swtch_t
        epsilon = swtch_eps;
        subplot(1,3,2);
        imshow(workim);
        title(["t = " num2str(i*dt)]);
        drawnow;
        subplot(1,3,3);
    end
    
    %TRANSFORM
    workim_tilda = fft2(workim);
    
    %POTENTIALS
    %F = 4*workim.^3 - 6*workim.^2 + 2*workim; %DOUBLE WELL u^2*(u-1)^2
    %F = -0.5*ones([n,m]).*heaviside(-workim) + (-workim/2 + ...
    %ones([n,m])/4).*heaviside(workim).*heaviside(ones([n,m])-workim) + ...
    %0.5*ones([n,m]).*heaviside(workim-ones([n,m])); %DOUBLE OBSTACLE 1/4*u*(1-u)
    %F = log((workim)./(ones([n,m])-workim))/4 - 2*workim + ones([n,m]); %LOG 1/4*(u*log(u/(1-u)) + log(1-u) + log(2) - 4*u^2 + 4*u -1)
    
    
    %SHOCK FILTER POTENTIAL
    gradxworkim = imfilter(workim, gradx);
    gradyworkim = imfilter(workim, grady);
    agworkim = sqrt(gradxworkim.^2 + gradyworkim.^2);
   
    laplaceworkim = imfilter(workim, laplace);
    funclaplace = laplaceworkim./(abs(laplaceworkim) + delta);
    F = -agworkim.*funclaplace;
    
    
    F_tilda = fft2(F);
    
    %FIDELITY TERM
    lambda_u = lambda.*(im-workim);
    lambda_u_tilda = fft2(lambda_u);
    
    %RIGHT HAND SIDE OF EQUATION
    rhs = workim_tilda - (dt/epsilon)*K2.*F_tilda + dt*(lambda_u_tilda + C1*K2.*workim_tilda + C2*workim_tilda);
    lhs = 1*ones([n,m]) + dt*(epsilon*K4 + C1*K2 + C2*ones([n,m]));
    
    %EVOLUTION
    workim_tilda = rhs./lhs;
    workim = real(ifft2(workim_tilda));
    
    %PLOT
    imshow(workim);
    title(["t = " num2str(i*dt)]);
    drawnow;
end
