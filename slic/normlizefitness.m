function [Cfit]=normlizefitness(Cfit)
    fitness=Cfit(:,2);
    fitness = setdiff(fitness,0);
    fitmin=min(fitness);
    fitmax=max(fitness);
    fitrange=fitmax-fitmin;
    Cfit=sortrows(Cfit,-2);
    if(fitrange>0)
        Cfit(:,2)=(Cfit(:,2)-fitmin)/fitrange;
    else
        Cfit(:,2)=0;
    end