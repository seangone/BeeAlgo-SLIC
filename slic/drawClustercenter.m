function im = drawClustercenter(C, im)
if exist('im', 'var') 
%     if ~exist('col', 'var'), col = [255 255 0]; end    
%     model=[ 0 1 0;
%             1 1 1;
%             0 1 0];
  
    [~,Cn]=size(C);
    for k=1:Cn
%       centercolor = 255-255*lab2rgb(C(1:3,k)');
        centercolor = [255 255 0];
        for c=1:3
            if C(4,k)~=0&&C(5,k)~=0
                im(C(4,k),C(5,k),c)=centercolor(c);
            end
        end
    end
end