function curIrre=LcalIrregularity(Centroid,curSP,sizeIM)
[y,x]=ind2sub(sizeIM,curSP);

distance=sqrt((Centroid(1)-x).^2+(Centroid(2)-y).^2);

curIrre=mean(abs(diff(distance)));

end