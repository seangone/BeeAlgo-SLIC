BW = logical([1 1 1 0 0 0 0 0
              1 1 1 0 1 1 0 0
              1 1 1 0 1 1 0 0
              1 0 1 0 0 0 1 0
              1 1 1 0 0 0 1 0
              1 1 1 0 0 0 1 0
              1 1 1 0 0 1 1 0
              1 1 1 0 0 0 0 0]);
[L,c] = bwlabel(BW,4)
[x,y] = find(L == 2);
se=[0 1 0;1 1 1;0 1 0];
% BW1=imclose(BW,se)
% imdilate(L,se)
L=cleanupregionsbyadjecentpx(L)