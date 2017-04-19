% SLIC Simple Linear Iterative Clustering SuperPixels
%
% Implementation of Achanta, Shaji, Smith, Lucchi, Fua and Susstrunk's
% SLIC Superpixels??
%
% Usage:   [l, Am, Sp, d] = slic(im, k, m, seRadius, colopt, mw)
%
% Arguments:  im - Image to be segmented.  欲分割的图像
%              k - Number of desired superpixels. Note  that this is nominal    
%                  the actual number of superpixels generated will generally      
%                  be a bit larger, espiecially if parameter m is small. 
%k   想要的超像素个数   生成的实际上的超像素个数会稍微多，在参数太小的情况
%              m - Weighting factor between colour and spatial   
%                  differences. Values from about 5 to 40 are useful.  Use a
%                  large value to enforce superpixels with more regular and
%                  smoother shapes. Try a value of 10 to start with.  
%m   为在颜色和空间上的差异的权重参数。范围是[5,40]，使用一个较大的值去保证超像素有更规则柔滑的形状
%       seRadius - Regions morphologically smaller than this are merged with
%                  adjacent regions. Try a value of 1 or 1.5.  Use 0 to
%                  disable. 
%seRadius  形态学上比这个参数更小的区域被与临近的区域合并
%         colopt - String 'mean' or 'median' indicating how the cluster
%                  colour centre should be computed. Defaults to 'mean'
%colopt   决定聚类中心是如何计算的，默认使用平均方法，可以选择平均或者中位数
%             mw - Optional median filtering window size.  Image compression
%                  can result in noticeable artifacts in the a*b* components
%                  of the image.  Median filtering can reduce this. mw can be
%                  a single value in which case the same median filtering is
%                  applied to each L* a* and b* components.  Alternatively it
%                  can be a 2-vector where mw(1) specifies the median
%                  filtering window to be applied to L* and mw(2) is the
%                  median filtering window to be applied to a* and b*.
%mw    可选择的过滤窗的大小。
% Returns:     l - Labeled image of superpixels. Labels range from 1 to k.  
%返回值：    l是标签好的超像素，标签取值[1,k]
%             Am - Adjacency matrix of segments.  Am(i, j) indicates whether  
%                  segments labeled i and j are connected/adjacent   
%   Am    各部分的相邻关系矩阵。   
%  Am(i,j)说明标签为i和j的部分是否相连
%             Sp - Superpixel attribute structure array with fields:  
%  Sp    超像素属性结构体数组
%                   L  - Mean L* value  
%                   a  - Mean a* value
%                   b  - Mean b* value
%                   r  - Mean row value
%                   c  - Mean column value
%                   stdL  - Standard deviation of L* 
%                   stda  - Standard deviation of a* 
%                   stdb  - Standard deviation of b* 
%                   N - Number of pixels
%                   edges - List of edge numbers that bound each
%                           superpixel. This field is allocated, but not set,
%                           by SLIC. Use SPEDGES for this.
%              d - Distance image giving the distance each pixel is from its
%                  associated superpixel centre.
%
% It is suggested that use of this function is followed by SPDBSCAN to perform a
% DBSCAN clustering of superpixels.  This results in a simple and fast
% segmentation of an image.  
%建议使用这个函数之后，使用SPDBSCAN函方法去显示DBSCAN超像素聚类，这会让图像分割更快更简单
% Minor variations from the original algorithm as defined in Achanta et al's
% paper:
%
% - SuperPixel centres are initialised on a hexagonal grid rather than a square
%   one. This results in a segmentation that will be nominally 6-connected
%   which hopefully facilitates any subsequent post-processing that seeks to
%   merge superpixels.
% - Initial cluster positions are not shifted to point of lowest gradient
%   within a 3x3 neighbourhood because this will be rendered irrelevant the
%   first time cluster centres are updated.
%
% Reference: R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua and
% S. Susstrunk. "SLIC Superpixels Compared to State-of-the-Art Superpixel
% Methods"  PAMI. Vol 34 No 11.  November 2012. pp 2274-2281.
%
% See also: SPDBSCAN, MCLEANUPREGIONS, REGIONADJACENCY, DRAWREGIONBOUNDARIES, RGB2LAB

% Copyright (c) 2013 Peter Kovesi
% Centre for Exploration Targeting
% School of Earth and Environment
% The University of Western Australia
% peter.kovesi at uwa edu au
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% Feb  2013
% July 2013 Super pixel attributes returned as a structure array

% Note that most of the computation time is not in the clustering, but rather
% in the region cleanup process.


function [C, Cfit, l, dengrid] = sliconebee(im, mask, segtimes, denmin, m, colopt, mw, nItr, eim, We)
	%设置参数的默认值
    if ~exist('colopt','var') || isempty(colopt), colopt = 'mean'; end
    if ~exist('mw','var')     || isempty(mw),         mw = 0;      end
    if ~exist('nItr','var')   || isempty(nItr),     nItr = 10;     end
    
    if exist('eim', 'var'), USEDIST = 0; else, USEDIST = 1; end
        
    MEANCENTRE = 1;
    MEDIANCENTRE = 2;
    
    if strcmp(colopt, 'mean')
        centre = MEANCENTRE;
    elseif strcmp(colopt, 'median')
        centre = MEDIANCENTRE;        
    else
        error('Invalid colour centre computation option');
    end
    
    [rows, cols, chan] = size(im);
    if chan ~= 3
        error('Image must be colour');
    end
%     m=m*(log(segtimes)+1);
%     m=m*3^(segtimes-1);

    
    %%
    l = -ones(rows, cols);   % Pixel labels.
    %%
    [trmin trmax tcmin tcmax]=findtarget(mask);
    imtarget   = im   (trmin:trmax, tcmin:tcmax, :);
    masktarget = mask (trmin:trmax, tcmin:tcmax, :);
    ltarget    = l    (trmin:trmax, tcmin:tcmax, :);
    save('a.mat','imtarget','masktarget');
    [rows, cols, chan] = size(imtarget);
% imtarget=im;
% masktarget=mask;
% ltarget=l;
    %%
    % Apply median filtering to colour components if mw has been supplied
    % and/or non-zero
    if mw
        if length(mw) == 1
            mw(2) = mw(1);  % Use same filtering for L and chrominance
        end
        for n = 1:3
            imtarget(:,:,n) = medfilt2(imtarget(:,:,n), [mw(1) mw(1)]);
        end
    end
    
    %%
    [C Cfit dengrid S vSpacing] = initclustercenter(imtarget, masktarget, denmin, segtimes);
    if Cfit==-1
        l=0;lclear=0;Am=0;Sp=0;
        C=0;
        return;
    end
    if C==0
        l=0;lclear=0;Am=0;Sp=0;
        C=0;
        return;
    end

    [~,Ck]=size(C);
    d = inf(rows, cols);     % Pixel distances from cluster centres.
    dc2 = inf(rows, cols);     % Pixel distances from cluster centres.
    %% 开始聚类算法，nItr是聚类算法迭代次数
    % Now perform the clustering.  10 iterations is suggested but I suspect n
    % could be as small as 2 or even 1
    S = round(S);  % We need S to be an integer from now on
	vSpacing=round(vSpacing);
    for n = 1:nItr
       for kk = 1:Ck  % for each cluster
           
           % Get subimage around cluster
 		   rmin = max(C(4,kk)-S, 1);          rmax = min(C(4,kk)+S, rows); 
           cmin = max(C(5,kk)-S, 1);          cmax = min(C(5,kk)+S, cols); 
           subim = imtarget(rmin:rmax, cmin:cmax, :);  
           assert(numel(subim) > 0)
           % Compute distances D between C(:,kk) and subimage
		   % 以下是自定义的计算距离的函数
		   % m为颜色值的权重
           if USEDIST
               [D DC2] = dist(C(:, kk), subim, rmin, cmin, S, m);
           else
               D = dist2(C(:, kk), subim, rmin, cmin, S, m, eim, We);
           end

           % If any pixel distance from the cluster centre is less than its
           % previous value, update its distance and label
		   %把像素与聚类中心的距离矩阵d更新为较小值
		   %把像素对应的聚类中心矩阵l更新为聚类标号
           subd =  d(rmin:rmax, cmin:cmax);
           subdc2= d(rmin:rmax, cmin:cmax);
           subl =  ltarget(rmin:rmax, cmin:cmax);
           updateMask = D < subd;
           subd(updateMask) = D(updateMask);
           subdc2(updateMask) = DC2(updateMask);
           subl(updateMask) = kk;
           
           Cfit(kk,3)=sum(subdc2(updateMask));
%            [gradxsubdc gradysubdc]=gradient(sqrt(subdc2));
%            Cfit(kk,1)=sum(gradxsubdc(updateMask).^2+gradysubdc(updateMask).^2);
           
           d(rmin:rmax, cmin:cmax) = subd;
           dc2(rmin:rmax, cmin:cmax) = subdc2;
           ltarget(rmin:rmax, cmin:cmax) = subl;
       end
       
       % Update cluster centres with mean values
	   %使用均值来更新聚类中心
       C(:) = 0;
       for r = 1:rows
           for c = 1:cols
              if ltarget(r,c)>0
                  tmp = [imtarget(r,c,1); imtarget(r,c,2); imtarget(r,c,3); r; c; 1];
                  C(:, ltarget(r,c)) = C(:, ltarget(r,c)) + tmp;
              end
           end
       end
       
       % Divide by number of pixels in each superpixel to get mean values
       for kk = 1:Ck 
           C(1:5,kk) = round(C(1:5,kk)/C(6,kk));
           Cfit(kk,3)= sqrt(Cfit(kk,3)/C(6,kk));
           Cfit(kk,2)= Cfit(kk,3)/(1+1*log10(Cfit(kk,4)));
%            Cfit(kk,2)= Cfit(kk,3)*1*log(C(6,kk))/(1+1*log10(Cfit(kk,4)));
       end
        
       % Note the residual error, E, is not calculated because we are using a
       % fixed number of iterations 
    end %n

%     ltarget = cleanupregionsbyadjecentpx(ltarget);
%%
% 统一处理目标区域的复原
%     [trmin trmax tcmin tcmax]
    im   (trmin:trmax, tcmin:tcmax, :) = imtarget;
    mask (trmin:trmax, tcmin:tcmax, :) = masktarget;
    ltarget(~masktarget)=-1;
    l    (trmin:trmax, tcmin:tcmax, :) = ltarget;

    for kk = 1:Ck 
       C(4,kk) = C(4,kk) + trmin; 
       C(5,kk) = C(5,kk) + tcmin;        
    end
%     im=imtarget;
%     mask=masktarget;
%     l=ltarget;
    %%
    
%-- dist -------------------------------------------
%
% Usage:  D = dist(C, im, r1, c1, S, m)
% 
% Arguments:   C - Cluster being considered
%             im - sub-image surrounding cluster centre
%         r1, c1 - row and column of top left corner of sub image within the
%                  overall image.
%              S - grid spacing
%              m - weighting factor between colour and spatial differences.
%
% Returns:     D - Distance image giving distance of every pixel in the
%                  subimage from the cluster centre
%
% Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% where:
% dc = sqrt(dl^2 + da^2 + db^2)  % Colour distance
% ds = sqrt(dx^2 + dy^2)         % Spatial distance
%
% m is a weighting factor representing the nominal maximum colour distance
% expected so that one can rank colour similarity relative to distance
% similarity.  try m in the range [1-40] for L*a*b* space
%
% ?? Might be worth trying the Geometric Mean instead ??
%  Distance = sqrt(dc * ds)
% but having a factor 'm' to play with is probably handy

% This code could be more efficient

function [D dc2 ds2] = dist(C, im, r1, c1, S, m)

    % Squared spatial distance
    %    ds is a fixed 'image' we should be able to exploit this
    %    and use a fixed meshgrid for much of the time somehow...
    if(exist('m','var'))
        [rows, cols, chan] = size(im);
        [x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
        x = x-C(5);  % x and y dist from cluster centre
        y = y-C(4);
        ds2 = x.^2 + y.^2;
        ds2 = ds2/S^2*m^2;
    else
        ds2 = 0;
    end
    
    % Squared colour difference
    for n = 1:3
        im(:,:,n) = (im(:,:,n)-C(n)).^2;
    end
    dc2 = sum(im,3);
%     dc = sqrt(dc2);
%     ds = sqrt(ds2);
    D = sqrt(dc2 + ds2);
    
    
    
%--- dist2 ------------------------------------------
%
% Usage:  D = dist2(C, im, r1, c1, S, m, eim)
% 
% Arguments:   C - Cluster being considered
%             im - sub-image surrounding cluster centre
%         r1, c1 - row and column of top left corner of sub image within the
%                  overall image.
%              S - grid spacing
%              m - weighting factor between colour and spatial differences.
%            eim - Edge strength sub-image corresponding to im
%
% Returns:     D - Distance image giving distance of every pixel in the
%                  subimage from the cluster centre
%
% Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% where:
% dc = sqrt(dl^2 + da^2 + db^2)  % Colour distance
% ds = sqrt(dx^2 + dy^2)         % Spatial distance
%
% m is a weighting factor representing the nominal maximum colour distance
% expected so that one can rank colour similarity relative to distance
% similarity.  try m in the range [1-40] for L*a*b* space
%

function D = dist2(C, im, r1, c1, S, m, eim, We)

    % Squared spatial distance
    %    ds is a fixed 'image' we should be able to exploit this
    %    and use a fixed meshgrid for much of the time somehow...
	%计算每个像素点到聚类中心的距离
    [rows, cols, chan] = size(im);
    [r,c] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
    dr = r-C(4);
    dc = c-C(5);
    ds2 = dr.^2 + dc.^2;
    
    % Squared colour difference
	%计算每个像素点与聚类中心的颜色值距离
    for n = 1:3
        im(:,:,n) = (im(:,:,n)-C(n)).^2;
    end
    dc2 = sum(im,3);
    
    % Combine colour and spatial distance measure
    D = sqrt(dc2 + ds2/S^2*m^2);
    
    % for every pixel in the subimage call improfile to the cluster centre
    % and use the largest value as the 'edge distance'
    rCentre = C(4)-r1;   % Cluster centre coords relative to this sub-image
    cCentre = C(5)-c1;
    de = zeros(rows,cols);
    for r = 1:rows
        for c = 1:cols
            v = improfile(eim,[c cCentre], [r rCentre]);
            de(r,c) = max(v);
        end
    end

    % Combine edge distance with weight, We with total Distance.
    D = D + We * de;
    
