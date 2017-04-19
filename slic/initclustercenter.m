function [C Cfit dengrid S vSpacing] = initclustercenter(im, mask, denmin, segtimes)
    [rows, cols, chan] = size(im);
    gridpxcnt=rows*cols;
    maskpxcnt=sum(sum(mask));
    if maskpxcnt<1 
        C=0;
        Cfit=0;
        dengrid=0;
        S=0;
        vSpacing=0;
        return;
    end
    %%
    %den为每10000个像素拥有的聚类个数
    %k为聚类个数
    dengrid = denmin;
    den=0;
    trygridtimes=0;
    while(den<=denmin)
        if trygridtimes>5
            C=0;
            Cfit=-1;
            return;
        end
        trygridtimes=trygridtimes+1;

        kgrid   = gridpxcnt/10000 * dengrid;

        
        % Nominal spacing between grid elements assuming hexagonal grid
        %竖放六边形的占据宽度
        %k为超像素个数
        S = sqrt(gridpxcnt / (kgrid * sqrt(3)/2));
        if S<6 
            C=0;
            Cfit=0;
            vSpacing=0;
            return;
        end
        % Get nodes per row allowing a half column margin at one end that alternates
        % from row to row   nodeCols为每一行可以放下的六边形个数
        nodeCols = round(cols/S - 0.5);
        % Given an integer number of nodes per row recompute S
        % 调整S为整数
        S = cols/(nodeCols + 0.5); 
        % Get number of rows of nodes allowing 0.5 row margin top and bottom
        % nodeRows为每一列可以放下的六边形个数
        nodeRows = round(rows/(sqrt(3)/2*S));
        % vSpacing为六边形平均占据高度
        vSpacing = rows/nodeRows;
        % Recompute k
        % 调整k为整数
        kgrid = nodeRows * nodeCols;
        dengrid = kgrid/(gridpxcnt/10000);
        
        maskcluster=zeros(rows,cols);
        r = vSpacing/2;
        for ri = 1:nodeRows
            if mod(ri,2), c = S/2; else, c = S;  end
            for ci = 1:nodeCols
                cc = round(c); rr = round(r);
                if mask(rr,cc)==1
                    maskcluster(rr,cc)=1;
                end
                c = c+S;
            end
            r = r+vSpacing;
        end
        
        Ck = sum(sum(maskcluster));
        den = Ck/(maskpxcnt/10000);
        
        
        dengrid = dengrid*1.3;
    end %while

    %% 初始化聚类中心
    % Allocate memory and initialise clusters, labels and distances.
	% 分配内存，初始化聚类、标签、距离
	% C为聚类中心数据，1:3是平均Lab值，4:5是中心行列，6是像素序号
	%
	% l是每个像素的标签   d是每个像素和聚类中心的距离
    C = zeros(6,Ck);          % Cluster centre data  1:3 is mean Lab value,
                             % 4:5 is row, col of centre, 6 is count of pixels
    Cfit = zeros(Ck,5);     %Cfit为适应度数组，用于蜂群对聚类的分配
    
    % Initialise clusters on a hexagonal grid
	% r初始化为每列的首个六边形中心列坐标，初始化为半个六边形占据高度 
	% kk初始化为1，为像素的标签循环记号
    kk = 1;
    r = vSpacing/2;
    
    for ri = 1:nodeRows
        % Following code alternates the starting column for each row of grid
        % points to obtain a hexagonal pattern. Note S and vSpacing are kept
        % as doubles to prevent errors accumulating across the grid.
		% c初始化为每行的首个六边形中心横坐标
        if mod(ri,2), c = S/2; else, c = S;  end
        
        for ci = 1:nodeCols
            cc = round(c); rr = round(r);
             if mask(rr,cc)==1
                %聚类中心初始化为六边形的一个点的色彩数据
                %这里matlab会有问题，如果给的是lab，会出现负的色彩，就自动把C定为int16
                C(1:5, kk) = [squeeze(im(rr,cc,:)); rr; cc];
                Cfit(kk,1)=kk;
                Cfit(kk,4)=dengrid;
                Cfit(kk,5)=segtimes;
                kk = kk+1;
             end
             c = c+S;
        end
        
        r = r+vSpacing;
    end
