function [seg] = cleanupregionsbyadjecentpx(seg,THEROLD1,THEROLD2)
% seg=[ 1 1 1 4 4 4 1 4
%       1 1 1 4 2 2 4 4
%       1 1 1 4 2 2 4 4
%       1 4 1 4 4 4 3 4
%       1 1 1 4 1 4 3 4
%       1 1 1 4 4 0 3 4
%       1 1 1 4 4 3 3 4
%       1 1 1 3 4 4 4 4]
% seg=[
%      1     1     1     0     0     0     1     0
%      1     1     1     0     1     1     0     0
%      1     1     1     0     1     1     0     0
%      1     0     1     0     0     0     1     0
%      1     1     1     0     1     0     1     0
%      1     1     1     0     0     0     1     0
%      1     1     1     0     0     1     1     0
%      1     1     1     1     0     0     0     0]  
se=[1 1 1;1 1 1;1 1 1];
if ~exist('THEROLD1','var'),THEROLD1=16;end
if ~exist('THEROLD2','var'),THEROLD2=16;end
labels = unique(seg(:))';
for l = labels
    [conl,num] = bwlabel(seg==l, 4);  
    for n = 1:num
        spur = conl==n;
        if sum(sum(spur))<=THEROLD1
            
            spurdil = imdilate(spur,se);
%             spurdil = edge(spur,'prewitt',0.04);
            spurdil = seg(spurdil);
            spurdil2 = reshape(spurdil',1,length(spurdil));
            spurdillabels = unique(spurdil2);
            spurdillabels = setdiff(spurdillabels,l);
            maxcount=0;
            replacelabel=l;
            for l2 = spurdillabels
%                 if length(find(seg==l2))>THEROLD2
                    count=length(find(spurdil2==l2));
                    if maxcount <= count;
                        replacelabel=l2;
                        maxcount=count;
                    end
%                 end
            end
            seg(spur)=replacelabel;
            
        end
    end
end