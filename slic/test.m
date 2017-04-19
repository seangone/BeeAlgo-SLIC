im=imread('9.jpg');
[l,l2, Am, C] = slic(im, 3000, 10, 1, 'median',0,2);

figure(1)
imshow(drawregionboundaries(l, im, [255 255 255]))
figure(2)
imshow(drawregionboundaries(l2, im, [255 255 255]))
figure(3)
lc = spdbscan(l, C, Am, 5);
imshow(drawregionboundaries(lc, im, [255 255 255]))