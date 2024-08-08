clear,clc,close all;

R = 0.25;
featureS = 10e-3;

r = linspace(-0.25,0.25,50);
[X,Y] = meshgrid(r,r);
Circ = (X.^2+Y.^2)<=R^2;
Ring1= (X.^2+Y.^2)<=(R/2)^2;


S = zeros(2000,2000);



for m = 1:11
    for n = 1:11
        
        if m==6 && n==6
           S = drawdisk(S,Circ-Ring1,m*100+400,n*100+400,R,featureS);
        else
            S = drawdisk(S,Circ,m*100+400,n*100+400,R,featureS);
        end
    end
end
S = uint8(S)*255;

S2 = 255-S;
S3 = 255-S2;
imshow(S2);



function K = drawdisk(K,Circ,i,j,R,featureS)
    K(i-R/featureS:i+R/featureS-1,j-R/featureS:j+R/featureS-1) = Circ;
end







