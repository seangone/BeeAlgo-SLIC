function [colorfitness imwithfitness] = drawClusterfitness(l,Cfit, im, a)
    if ~exist('a','var')||a<0||a>1,a=0.3;end
    [rows,cols]=size(l);
    [Cn,~]=size(Cfit);
    if max(max(l))>Cn 
        error('Draw Mask error');
    end
    
    [Cfit]=normlizefitness(Cfit);
    
    Cfit=sortrows(Cfit,1);
    colorfitness=zeros(rows,cols,3);
    imwithfitness=zeros(rows,cols,3);
    mtimes=max(Cfit(:,5))-1;
    for r=1:rows
        for c=1:cols
            if l(r,c)
                fitnesscolor1=[0 0 255]*Cfit(l(r,c),2);fitnesscolor2=0;
                if mtimes>1
                    fitnesscolor2=[0 255 0]*(Cfit(l(r,c),5)-1)/mtimes;
                end
                colorfitness(r,c,:)=fitnesscolor1+fitnesscolor2;
                imwithfitness(r,c,:)=[255 255 255]-[0 255 255]*Cfit(l(r,c),2);
            end
        end
    end
    colorfitness=uint8(colorfitness);
    imwithfitness=uint8(imwithfitness);
    imwithfitness=imwithfitness*a+im*(1-a);

