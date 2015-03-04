function E = imageForces(I,Wline,Wedge,Wterm,avgSigma)
%%imageForces(I,Wline,Wedge,Wterm,avgSigma) returns the Image Energy
%%Input Arguments:
%%I: Input image
%%Wline: Weight for the Line energy 
%%Wedge: Weight for the Edge energy 
%%Wterm: Weight for the Termination energy
%%avgSigma: Sigma value for average filter

    %%Create a kernel for averaging
    H = fspecial('average', avgSigma);
    
    %%Find the first derivative of the Image after averaging it
    [X_firstDer,Y_firstDer] = imgradientxy(imfilter(I,H));   
    %%Magnitude for the first derivative
    delI = (X_firstDer.^2 + Y_firstDer.^2);
    
    %%One of the first methods to get the second derivative
%     [X_secondDer,YX_secondDer] = imgradientxy(X_firstDer);
%     [Y_secondDer,XY_secondDer] = imgradientxy(Y_firstDer);
%     delSquareI = X_secondDer + XY_secondDer + YX_secondDer + Y_secondDer;
     

    %%Create a gaussian kernel for smoothing
    H = fspecial('gaussian', 10, 10);
    
    %%Find the first and second derivatives of the image
    [Cx,Cy] = imgradientxy(imfilter(I,H));        
    [Cxx,Cyx] = imgradientxy(Cx);
    [Cxy,Cyy] = imgradientxy(Cy); 
    
    %%Termination Energy
    Eterm = (Cyy.*(Cx.^2) - 2*Cxy.*Cx.*Cy + Cxx.*(Cy.^2)) ./ (Cx.^2 + Cy.^2).^1.5;
    %%Remove NaN's
    Eterm(find(isnan(Eterm) == 1)) = 0;
    
    %%Edge Energy
    Eedge = -1*(delI);
    %%Remove NaN's
    Eedge(find(isnan(Eedge) == 1)) = 0;
    
    %%Line Energy
    Eline = double(I);
    %%Remove NaN's
    Eline(find(isnan(Eline) == 1)) = 0;

    %%Total Image Energy
    E = Wedge*Eedge + Wline*Eline + double(Wterm*uint8(Eterm));
end