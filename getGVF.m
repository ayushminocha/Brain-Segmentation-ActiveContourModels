function [u,v] = getGVF(Fx,Fy)
%%getGVF(Fx,Fy) return the GVF forces in the X and Y direction
%%Input Arguments:
%%Fx: Forces in the X direction
%%Fy: Forces in the Y direction
%%Output:
%%u: GVF in x direction
%%v: GVF in y direction

%Squared magnitude of the forces in X and Y direction
sMag = Fx.^2+ Fy.^2;

%%Initialize u and v with initial forces in the X and Y directions
u=Fx;  v=Fy;
  
for i=1:50
% Calculate Laplacian

%%Different ways to get the second derivative
%%1. Using imgradient
%   [Ux,Uy] = imgradientxy(u);
%   [Uxx,Uyx] = imgradientxy(Ux);
%   [Uyy,Uxy] = imgradientxy(Uy);
%  
%   [Vx,Vy] = imgradientxy(v);
%   [Vxx,Vyx] = imgradientxy(Vx);
%   [Vxy,Vyy] = imgradientxy(Vy);
% 
%   Uxx(find(isnan(Uxx) == 1)) = 0;
%   Uyy(find(isnan(Uyy) == 1)) = 0;  
%   Vxx(find(isnan(Vxx) == 1)) = 0;
%   Vyy(find(isnan(Vyy) == 1)) = 0;
%%----------------

%%2. using derivative of gaussian kernel
%   DGaussxx = 1/(2*pi*sigma^4) * (xx.^2/sigma^2 - 1) .* exp(-(xx.^2 + yy.^2)/(2*sigma^2));
%   DGaussyy = 1/(2*pi*sigma^4) * (yy.^2/sigma^2 - 1) .* exp(-(xx.^2 + yy.^2)/(2*sigma^2));
%   
%   Uxx = imfilter(u,DGaussxx,'conv','symmetric');
%   Uyy = imfilter(u,DGaussyy,'conv','symmetric');
%   Vxx = imfilter(v,DGaussxx,'conv','symmetric');
%   Vyy = imfilter(v,DGaussyy,'conv','symmetric');
%   
%   Uxx(find(isnan(Uxx) == 1)) = 0;
%   Uyy(find(isnan(Uyy) == 1)) = 0;  
%   Vxx(find(isnan(Vxx) == 1)) = 0;
%   Vyy(find(isnan(Vyy) == 1)) = 0;
%%----------------

%%3. Using del2(U) which returns a discrete approximation of Laplace's 
%%differential operator applied to U using the default spacing, h = 1, 
%%between all points.
%%----------------

% Update the vector field
%   u = u + 0.8*(Uxx+Uyy) - sMag.*(u-Fx);
%   v = v + 0.8*(Vxx+Vyy) - sMag.*(v-Fy);
  u = u + 0.8*del2(u) - sMag.*(u-Fx);
  v = v + 0.8*del2(v) - sMag.*(v-Fy);
end
end