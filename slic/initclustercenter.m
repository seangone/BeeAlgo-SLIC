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
    %denΪÿ10000������ӵ�еľ������
    %kΪ�������
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
        %���������ε�ռ�ݿ��
        %kΪ�����ظ���
        S = sqrt(gridpxcnt / (kgrid * sqrt(3)/2));
        if S<6 
            C=0;
            Cfit=0;
            vSpacing=0;
            return;
        end
        % Get nodes per row allowing a half column margin at one end that alternates
        % from row to row   nodeColsΪÿһ�п��Է��µ������θ���
        nodeCols = round(cols/S - 0.5);
        % Given an integer number of nodes per row recompute S
        % ����SΪ����
        S = cols/(nodeCols + 0.5); 
        % Get number of rows of nodes allowing 0.5 row margin top and bottom
        % nodeRowsΪÿһ�п��Է��µ������θ���
        nodeRows = round(rows/(sqrt(3)/2*S));
        % vSpacingΪ������ƽ��ռ�ݸ߶�
        vSpacing = rows/nodeRows;
        % Recompute k
        % ����kΪ����
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

    %% ��ʼ����������
    % Allocate memory and initialise clusters, labels and distances.
	% �����ڴ棬��ʼ�����ࡢ��ǩ������
	% CΪ�����������ݣ�1:3��ƽ��Labֵ��4:5���������У�6���������
	%
	% l��ÿ�����صı�ǩ   d��ÿ�����غ;������ĵľ���
    C = zeros(6,Ck);          % Cluster centre data  1:3 is mean Lab value,
                             % 4:5 is row, col of centre, 6 is count of pixels
    Cfit = zeros(Ck,5);     %CfitΪ��Ӧ�����飬���ڷ�Ⱥ�Ծ���ķ���
    
    % Initialise clusters on a hexagonal grid
	% r��ʼ��Ϊÿ�е��׸����������������꣬��ʼ��Ϊ���������ռ�ݸ߶� 
	% kk��ʼ��Ϊ1��Ϊ���صı�ǩѭ���Ǻ�
    kk = 1;
    r = vSpacing/2;
    
    for ri = 1:nodeRows
        % Following code alternates the starting column for each row of grid
        % points to obtain a hexagonal pattern. Note S and vSpacing are kept
        % as doubles to prevent errors accumulating across the grid.
		% c��ʼ��Ϊÿ�е��׸����������ĺ�����
        if mod(ri,2), c = S/2; else, c = S;  end
        
        for ci = 1:nodeCols
            cc = round(c); rr = round(r);
             if mask(rr,cc)==1
                %�������ĳ�ʼ��Ϊ�����ε�һ�����ɫ������
                %����matlab�������⣬���������lab������ָ���ɫ�ʣ����Զ���C��Ϊint16
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
