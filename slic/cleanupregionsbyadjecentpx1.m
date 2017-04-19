function [seg] = cleanupregionsbyadjecentpx1(seg,THEROLD1,THEROLD2)

% clc 
% clear all
% seg=[
%      1     1     1     0     0     0     1     0
%      1     1     1     0     0     1     0     0
%      1     1     1     0     1     1     0     0
%      1     0     1     0     1     1     0     0
%      1     1     1     0     1     0     1     0
%      1     1     1     0     0     0     1     0
%      1     1     1     0     0     1     1     0
%      1     1     1     1     0     0     0     0]  ;
se=[1 1 1;1 1 1;1 1 1];
if ~exist('THEROLD1','var'),THEROLD1=10;end
if ~exist('THEROLD2','var'),THEROLD2=4;end
for l = 0:1
    [conl,num] = bwlabel(seg==l, 4);  
        for n = 1:num
            spur = conl==n;
            if sum(sum(spur))<=THEROLD1
%                 spurdil = imdilate(spur,se);
                spurdil = edge(spur,'prewitt',0.04);
                spurdil = seg(spurdil);       
                spurdil2 = reshape(spurdil',1,length(spurdil));
                spurdillabels = unique(spurdil2);
                spurdillabels = setdiff(spurdillabels,l);
                maxcount=0;
                replacelabel=l;
                for l2 = spurdillabels
                    count=length(find(spurdil2==l2));
                    if maxcount <= count;
                        replacelabel=l2;
                        maxcount=count;
                    end
                end
                seg(spur)=replacelabel;
            end
        end
end