% turn to the 4 quarter range angle
function alpha360_m=LcalAccAngle4Direction(curPtsonLine,GC_x,GC_y)

    alpha_m=atan2(GC_y(curPtsonLine), GC_x(curPtsonLine));
    alpha_m=-alpha_m;
    alpha360_m=alpha_m;
    temp=find(alpha_m<0);
    for k=1:length(temp)
        alpha360_m(temp(k))=pi+alpha_m(temp(k))+pi;
    end
    
end