function snake(im,x,y,niter,greymatter_flag)


%%Load Input Image
%im = imread('inp6.png');

%%Resize the Image if it is very big
if(size(im,1) >= 1000 || size(im,2) >= 1000)
    im = imresize(im,0.5);
end

%%Convert to grayscale
if size(im,3) ~= 1
    im = rgb2gray(im);
end

%%Create a copy of the original Image
origIm = im;

%%Adjust contrast in image to get better edges
if(greymatter_flag == 1)
    im = imadjust(im,[0.6 0.8],[]);
    %im = imadjust(im,[0.4 0.5],[]);
end

%sx sy] = size(origIm);
%origIm = origIm(1:sx/2,1:sy/2);
%im = im(1:sx/2,1:sy/2);

%%Convert Image to double from uint8 for computations
im = double(im);

%%Display Image and take user input for initial mask
imshow(origIm,[]);
%[x,y] = ginput;

%%Resample points on the snake contour
[x,y] = resamplePoints(x,y);

%%Set parameters for the snake energy
alpha = 0.001; %%Controls the continuity of the snake
beta = 0.00001; %%Controls the curvature of the snake
gamma = 0.002; %%Time step for snake modification
iterations = niter;%300; %%Number of iterations

%%Create the pentadiagonal matrix P
N = length(x);
a = +1*(2*alpha+6*beta);
b = -1*(alpha+4*beta);
c = +1*beta;
P = diag(repmat(a,1,N));
P = P + diag(repmat(b,1,N-1), 1) + diag(   b, -N+1);
P = P + diag(repmat(b,1,N-1),-1) + diag(   b,  N-1);
P = P + diag(repmat(c,1,N-2), 2) + diag([c,c],-N+2);
P = P + diag(repmat(c,1,N-2),-2) + diag([c,c], N-2);

P = P + gamma*eye(N);
P = inv(P);


%%Set the weights for Line Energy, Edge Energy and the Termination Energy
%%Image forces
Wline = 0.04;
Wedge = 4.0;
Wterm = 0.001;
%%Sigma value for the average filter
avgSigma = 15;
%%Get the image forces for contrast adjusted image
Eimage = imageForces(im,Wline,Wedge,Wterm,avgSigma);

%%First derivatives (Gradient) for the Image Energy
[gx, gy] = imgradientxy(Eimage);

%%Normalize the first derivatives
gx = gx./norm(gx);
gy = gy./norm(gy);

%%Get the GVF forces
[gx, gy] = getGVF(gx,gy);

%%Sigma value for the average filter
avgSigma = 20;
%%Get the image forces for the original image
Eimage2 = imageForces(origIm,Wline,Wedge,Wterm,avgSigma);

%%First derivatives (Gradient) for the Image Energy
[gx2, gy2] = imgradientxy(Eimage2);

%%Normalize the first derivatives
gx2 = gx2./norm(gx2);
gy2 = gy2./norm(gy2);

%%Get the GVF forces
[gx2, gy2] = getGVF(gx2,gy2);

%%New Image forces in the x (gx) and y (gy) direction
gx = gx + 0.5*gx2;
gy = gy + 0.5*gy2;
%gy = ones(size(gy));
%subplot(1,2,1);imshow(gx,[]);
%subplot(1,2,2);
imshow(origIm,[]);hold on;h = plot([x;x(1)],[y;y(1)]);
pause;

for ii = 1:iterations

    % Calculate external force by interpolatig on the x,y points
    Fext1(:,1)=interp2(gx,x,y);
    Fext1(:,2)=interp2(gy,x,y);
    % Interp2, can give nan's if contour close to border
    Fext1(isnan(Fext1))=0;
    
    %%Find the normal to the snake contour
    %%a is used to set the neighborhood distance for finding the derivative
    a=1;    

    % Derivatives of contour
    n=length(x);
    f=(1:n)+a; f(f>n)=f(f>n)-n;
    b=(1:n);
    
    dx=x(f)-x(b);
    dy=y(f)-y(b);    
    % Normals of contourpoints
    l=sqrt(dx.^2+dy.^2);
    nx = -dy./l;
    ny =  dx./l;
    
    nx = 0.0022*nx;
    ny = 0.0022*ny;    
    
    %%To give direction to the baloon forces
%    bs = ones(nPoints,1);
%     for i=1:100
%         if(im(round(x(i)),round(y(i))) == 0)
%             bs(i) = bs(i)*-1;
%         end
%     end
        
    
    %%Iteratively changing the x and y positions of the snake contour
    x = P*(gamma*x - Fext1(:,1) + nx);
    y = P*(gamma*y - Fext1(:,2) + ny);  
    
    %%Border conditions for the snake
    for i=1:length(x)
        if x(i)<=0
            x(i) = 1;
        elseif x(i) > size(im,1)
            x(i) = size(im,1);
        end
        
        if y(i)<=0
            y(i) = 1;
        elseif y(i) > size(im,2)
            y(i) = size(im,2);
        end
    end    
    
    %%Resample equidistant points on the snake contour
    [x,y] = resamplePoints(x,y);
    
    %%Remove loops from the snake
     [xx,yy] = removeloops(x,y,im);
%     
%     %%Resample equidistant points on the snake contour
     if(length(xx) > 500)
        [x,y] = resamplePoints(xx,yy);
     end
    
    %%Show the snake movement
    if(ishandle(h)), delete(h), end
    h=plot([x;x(1)],[y;y(1)],'r'); drawnow
    
end
%%Show the final position of snake
if(ishandle(h)), delete(h), end
h=plot([x;x(1)],[y;y(1)],'r'); drawnow
end