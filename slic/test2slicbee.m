%5,8为人像
%4为闪电图9为汽车
im=imread('4.jpg');
[C,Cn,Cfit,l,l2, Am] = slicbee(im, 300, 3000, 10, 1, 'median',0,2);

%%
[mask,imwithmask] = drawClusterMask(l,Cfit, im,0.65);
%%
figure(1)
im0=drawClustercenter(C,Cn,im);
im1=drawregionboundaries(l, im0, [255 255 255]);
imshow(im1);

% figure(2)
% lc = spdbscan(l, C, Am, 5);
% imshow(drawregionboundaries(lc, im0, [255 255 255]))

figure(3)
im0=drawClustercenter(C,Cn,imwithmask);
im1=drawregionboundaries(l, im0, [255 255 255]);
imshow(im1);

figure(4)
im0=drawClustercenter(C,Cn,mask);
im1=drawregionboundaries(l, im0, [255 255 255]);
imshow(im1);