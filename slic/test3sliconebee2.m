%%
imrgb=imread('13.jpg');
% Convert image to L*a*b* colourspace.  This gives us a colourspace that is
% nominally perceptually uniform. This allows us to use the euclidean
% distance between colour coordinates to measure differences between
% colours.  Note the image becomes double after conversion.  We may want to
% go to signed shorts to save memory.
imlab = rgb2lab(imrgb);
[rows, cols, ~] = size(imlab);

dengrid = Cn/(rows*cols/10000);
mask=ones(rows,cols);
[C1,Cfit1,l1,dengrid] = sliconebee(imlab, mask, 1, dengrid, 15, 'median',0,3);
l1=cleanupregionsbyadjecentpx(l1);
%%
% 
% figure(2)
% imrgb0=drawClustercenter(C1,imrgb);
% imrgb1=drawregionboundaries(l1, imrgb0, [255 255 255]);
% imshow(imrgb1);