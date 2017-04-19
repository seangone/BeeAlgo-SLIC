function [rmin rmax cmin cmax]=findtarget(mask)
    [rows,cols]=size(mask);
%     masktarget=zeros(rows,cols);
    rmin=1;cmin=1;rmax=rows;cmax=cols;
    while(sum(mask(rmin,:))==0)
        rmin=rmin+1;
    end
    while(sum(mask(:,cmin))==0)
        cmin=cmin+1;
    end
    while(sum(mask(rmax,:))==0)
        rmax=rmax-1;
    end
    while(sum(mask(:,cmax))==0)
        cmax=cmax-1;
    end
%     masktarget(rmin:rmax,cmin:cmax)=1;