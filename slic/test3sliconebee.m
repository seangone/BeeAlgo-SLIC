clc
clear all;
%%
imrgb=imread('13.jpg');
% Convert image to L*a*b* colourspace.  This gives us a colourspace that is
% nominally perceptually uniform. This allows us to use the euclidean
% distance between colour coordinates to measure differences between
% colours.  Note the image becomes double after conversion.  We may want to
% go to signed shorts to save memory.
imlab = rgb2lab(imrgb);
[rows, cols, ~] = size(imlab);
%%
dengrid = 1;
mask=ones(rows,cols);
[C1,Cfit1,l,dengrid] = sliconebee(imlab, mask, 1, dengrid, 15, 'median',0,3);
C=C1;
Cfit=Cfit1;
for i=1:3
    [~,Cn]=size(C);
    Cfitsorted=sortrows(Cfit,-2);
    for kk=1:round(Cn*0.2+0.5)
        lk=Cfitsorted(kk,1);
        mask=(l==lk);
        [C2,Cfit2,l2,~] = sliconebee(imlab, mask, Cfit(lk,5)+1,Cfit(lk,4)*4, 15, 'median',0,2);
        if Cfit2==-1
            Cfit(lk,2)=0;
            continue;
        end
        [~,Cn] = size(C);
        [~,Cn2] =size(C2);
        if Cn2>2
            %C的合并
            C = horzcat(C,C2);
            C(:,lk)=C2(:,Cn2);
            C=C(:,1:Cn+Cn2-1);
            %Cfit的合并
            Cfit = vertcat(Cfit,Cfit2);
            Cfit(lk,:)=Cfit2(Cn2,:);
            Cfit=Cfit(1:Cn+Cn2-1, :);
            Cfit(1:Cn+Cn2-1, 1) = 1:Cn+Cn2-1;
            %l的合并
            l2(mask)=l2(mask)+Cn;
            l2(l2==Cn+Cn2)=lk;
            l(mask)=l2(mask);
        end
    end
end

l=cleanupregionsbyadjecentpx(l);

%%
[colorfitness,imwithfitness] = drawClusterfitness(l,Cfit, imrgb, 0.55);
%%
figure(1)
imrgb0=drawClustercenter(C,imrgb);
imrgb1=drawregionboundaries(l, imrgb0, [255 255 255]);
imshow(imrgb1);


% figure(2)
% lc = spdbscan(l, C, Am, 5);
% imshow(drawregionboundaries(lc, im0, [255 255 255]))

figure(3)
imrgb0=drawClustercenter(C,imwithfitness);
imrgb1=drawregionboundaries(l, imrgb0, [255 255 255]);
imshow(imrgb1);

figure(4)
imrgb0=drawClustercenter(C,colorfitness);
imrgb1=drawregionboundaries(l, imrgb0, [255 255 255]);
imshow(imrgb1);