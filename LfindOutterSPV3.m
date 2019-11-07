%% V2 is more accurate
%  V3 consider the cases that the SP on other nuclei /not use thought,add
%  regularize

function AllSP=LfindOutterSPV3(Centroid,Area,im,bw,bw_allNuclei,GradienMapTpye,MeanIntensity, Regularization,curNucleiInd,shown)
%% decide the radi_maxal lines
% sample circular points around the center
RlineNO=50;
theta = linspace(0,2*pi,RlineNO);
theta(end)=[];

C_x=Centroid(1);
C_y=Centroid(2);

radi_max =ceil(2.6*sqrt(Area/pi));
% radi_max = 20;% should change with the size of the object

%%%  Threshold used to check if this SP to the centroid has a low mean
%%% value, if so disable this SP
TMI4RL=MeanIntensity*1.3;% threshold of the mean intensity for radial line

TcheckIfSPonOtherNuclei=3;% threshold to check if the SP is on the other nuclei region.
AllNucleiPixelIndx=find(bw_allNuclei==1);% for latter checking

%% supress the weak and noisy gradient by gaussian bluring
if strcmp(GradienMapTpye,'GaussianBlur')
    %     disp(' Compute edge map ...');
    f = double(im)/255;
    f0 = gaussianBlur(f,1);
    %     disp(' Comute the traditional external force ...');
    [px,py] = gradient(f0);
    %%% make the border gradient to be zeros
    px([1:2, end-1:end],:)=0;py([1:2, end-1:end],:)=0;
    px(:,[1:2, end-1:end])=0;py(:,[1:2, end-1:end])=0;
    
    %%% supress the weak gradient
    G_mag=sqrt(px.^2+ py.^2);
    % TGMag=mean(G_mag(:));   %max(G_mag(:));%figure(23);hist(G_mag,50);
    %       allMag=G_mag(:);allMag(find(allMag>150/255))=[];
    %       allMag(find(allMag<10/255))=[];
    %% !!!! important thrshold for  supressing the weak gradient
    TGMag=mean(G_mag(:));
    ValidGM=G_mag>=TGMag;
    px(~ValidGM)=0;   py(~ValidGM)=0;
    if shown
        % display the results
        figure(2);imshow(im,'InitialMagnification','fit');hold on;
        quiver(px,py,'y');
        hold off;
    end
    
    % changed gradient map
    GC_y=py; GC_x=px;
end

%% supress the weak gradient, need better method here
if strcmp(GradienMapTpye,'Threshold')
    [G_x,G_y]=gradient(double(im));
    G_mag=sqrt(G_x.^2+ G_y.^2);
    % TGMag=mean(G_mag(:));   %max(G_mag(:));%figure(23);hist(G_mag,50);
    allMag=G_mag(:);allMag(find(allMag>150))=[];
    allMag(find(allMag<10))=[];
    %% !!!! important thrshold for  supressing the weak gradient
    TGMag=mean(allMag)/3;
    
    %figure(23);hist(allMag,50);
    % TGMag=graythresh(allMag)*255;
    ValidGM=G_mag>=TGMag;
    
    G_x4shown=G_x;G_y4shown=G_y;
    G_x4shown(~ValidGM)=0;G_y4shown(~ValidGM)=0;
    
    if shown
        figure(1);imshow(im,'InitialMagnification','fit');hold on;
        quiver(G_x4shown,G_y4shown,'y');
        hold off;
    end
    
    
    % changed gradient map
    GC_y=G_y4shown; GC_x=G_x4shown;
end

% radi_min=10;% should be change with the obj boundary
%%
xi =C_x + radi_max * cos( theta );
yi =C_y + radi_max * sin( theta );
% plot circle points
if shown
    figure(1);imshow(im,'InitialMagnification','fit');hold on;
    quiver(GC_x,GC_y,'y'); %hold off;
    for i=1:length(theta)
        %     plot( C_x, C_y, xi(i),yi(i), 'y' );
        plot( [C_x,xi(i)], [C_y,yi(i)], 'y' );
    end
    hold off;
end
%% convert the continous line to discrete points
curSPt=[C_y,C_x];
for i=1:length(theta)
    curEPt=[yi(i),xi(i)];
    
    curPts=LgetLineSegmentbyTwoPts_light(curSPt,curEPt,size(bw));
    
    showbw=im;showbw(:)=0;
    curPtsInd=sub2ind(size(showbw),curPts(:,1),curPts(:,2));
    curPtsInd=unique(curPtsInd);
    PtsonLineSegment_full{i}=curPtsInd;
    
    showbw(curPtsInd)=1;
    %%% the pts outside the obj region is useful
    showbw=showbw&~bw;
    %%%% check if the line is being break into two, if so, keep the near
    %%%% centroid one
    cc=bwconncomp(showbw,8);
    stats=regionprops(cc,'Centroid');
    if cc.NumObjects>1
        %         show(showbw,1);
        curdistkk=[];
        for kk=1:cc.NumObjects
            curC=stats(kk).Centroid;
            curdistkk(kk)=norm(curC-[C_x,C_y]);
        end
        [minval,minind]=min(curdistkk);
        showbw = ismember(labelmatrix(cc), minind);
    end
    
    curValidPtsInd=find(showbw==1);
    %     show(~bw,1);
    
    %%% order them based on the distance to the centroid, the first element
    %%% is the most close one
    [SubInd_y,SubInd_x]=ind2sub(size(bw),curValidPtsInd);
    % cal the distance to the centroid
    distance=sqrt((C_y-SubInd_y).^2+(C_x-SubInd_x).^2);
    [sval,sortind]=sort(distance,'ascend');
    curValidPtsInd_s=curValidPtsInd(sortind);
    
    PtsonLineSegment{i}=curValidPtsInd_s;
    
    %check
    %     show(showbw,1);
    %     pause();
end
% hold off;

%%  for each radi_maxal line check out the cos(theta_n-alpha_m) for the mth pts on the line

for i=1:length(theta)
    curPtsonLine=PtsonLineSegment{i};
    if length(curPtsonLine)>1
        % turn to the 4 quarter range angle
        alpha360_m=LcalAccAngle4Direction(curPtsonLine,GC_x,GC_y);
        
        magnitude_m=sqrt(GC_x(curPtsonLine).^2+ GC_y(curPtsonLine).^2);
        if max(magnitude_m)~=0
            NormMag_m=magnitude_m/max(magnitude_m);
            ZeroMagInd= find(NormMag_m==0);
        else
            ZeroMagInd=1:length(magnitude_m);
        end
        
        if i==1
            angleDiff=alpha360_m-(theta(i));
            angleDiff(ZeroMagInd)=NaN;   % there is no angle different for the zero magnitude pts
        else
            angleDiff=alpha360_m-(2*pi-theta(i));% be careful of this
            angleDiff(ZeroMagInd)=NaN;   % there is no angle different for the zero magnitude pts
        end
        
        imBK=im;
        imBK( curPtsonLine)=0;
        if shown
            figure(1);imshow(imBK,'InitialMagnification','fit');hold on;
            quiver(GC_x,GC_y,'y');
            hold off;
        end
        
        %%% cal the parameter cos(angleDiff)+magnitude for all the pts have opposit gradient direction
        %%% get the opposit gradient direction
        if i==1
            OGInd=find(angleDiff>pi/2&angleDiff<pi*3/2);
        else
            %             OGInd=find(abs(angleDiff)>pi/2&abs(angleDiff)<pi*3/2);
            OGInd=find(abs(angleDiff)>pi/2);
        end
        
        if ~isempty(OGInd)
            %         AllPV=cos(abs(angleDiff))+NormMag_m;
            %         [maxVal,maxInd]=max(AllPV(OGInd));
            %%% consider the magnitude only
            %             [maxVal,maxInd]=max(NormMag_m(OGInd));
            %%% consider the combination of magnitude and direction
            AllPV=abs(cos(abs(angleDiff)))+NormMag_m;% all values for the pts on radial line
            [maxVal,maxInd]=max(AllPV(OGInd));
            %%% consider the combination of magnitude and direction and the
            %%% distance to the centroid
            %             curdist=linspace(1,0,length(angleDiff));
            %
            %             AllPV=abs(cos(abs(angleDiff)))+NormMag_m.*curdist';
            %             [maxVal,maxInd]=max(AllPV(OGInd));
            
            %%  further check if this SP to the centroid has a low mean
            %%% value, if so disable this SP
            curI4RL=im(curPtsonLine(1:OGInd(maxInd)-1));% get the intensity for the pontential supporting line
            curMI4RL=mean(curI4RL);
            if curMI4RL<TMI4RL
                maxInd=1;
                curSP=curPtsonLine(1);
            else
                curSP=curPtsonLine(OGInd(maxInd));% the support pts at this radial line
            end
            %% further check if the SP is on other nuclei, if so disable this SP and find another one*
%             if maxInd~=1
%                 [curSP_y,curSP_x]=ind2sub(size(bw),curSP);
%                 %%% exclude the cur nuclei
%                 AllNucleiPixelIndxExceptItself=setdiff(AllNucleiPixelIndx,curNucleiInd);
%                 [SubInd4Nuclei_y,SubInd4Nuclei_x]=ind2sub(size(bw),AllNucleiPixelIndxExceptItself);
%                 SP2NucleiDist=sqrt((curSP_y-SubInd4Nuclei_y).^2+(curSP_x-SubInd4Nuclei_x).^2);
%                 [minval,minind]=min(SP2NucleiDist);
%                 % if the SP is on other nuclei regions, disable it
%                 if minval<TcheckIfSPonOtherNuclei
%                     curSP=curPtsonLine(1);
%                 end
%             end
        else% if there is no support pt, set the support pt to the one nearest to the centroid
            %             [SubInd_y,SubInd_x]=ind2sub(size(bw),curPtsonLine);
            %             % cal the distance to the centroid
            %             distance=sqrt((C_y-SubInd_y).^2+(C_x-SubInd_x).^2);
            %             [minval,minind]=min(distance);
            curSP=curPtsonLine(1);
        end
    else
        if isempty(curPtsonLine)
            if i~=1
                curSP=AllSP(i-1);
            else
                temp=PtsonLineSegment_full{i};
                curSP=temp(end);
            end
        else
            curSP=curPtsonLine(1);
        end
    end
    
    %     imBK=im;
    if shown
        imBK(curSP)=160;
        figure(1);imshow(imBK,'InitialMagnification','fit');hold on;
        quiver(GC_x,GC_y,'y');
        hold off;
    end
    
    AllSP(i)=curSP;
end

%% do the regularization here
if Regularization==1
    if shown
        figure(11);imshow(im,'InitialMagnification','fit');hold on;
        for i=1:length(AllSP)
            [curSubSP_r,curSubSP_c]=ind2sub(size(bw),AllSP);
            curSubSP_r=[curSubSP_r curSubSP_r(1)];
            curSubSP_c=[curSubSP_c curSubSP_c(1)];
            
            plot(curSubSP_c,curSubSP_r,'y');
            
        end
        hold off;
    end
    
    %%% turn to the distance to centroid
    [curSubSP_r,curSubSP_c]=ind2sub(size(bw),AllSP);
    distance2C=sqrt( (C_x-curSubSP_c).^2+(C_y-curSubSP_r).^2);
    distance2C_Tri=[distance2C,distance2C,distance2C];
    % distance2C_s=smooth(distance2C,'moving');
    distance2C_s_Tri=smooth(distance2C_Tri,'rlowess');
    distance2C_s=distance2C_s_Tri(length(distance2C)+1:2*length(distance2C));
    
    % shown=1;
    if shown
        figure(12);plot(distance2C);hold on; plot(distance2C_s,'r');hold off;
    end
    %%%% map the smooth distance back to the radial line
    TdiffbetweenSmoothandOri=3;
    for i=1:length(theta)
        curdist=distance2C(i);
        curdist_s=distance2C_s(i);
        if abs(curdist-curdist_s)>TdiffbetweenSmoothandOri
            curPtsonLine=PtsonLineSegment{i};
            if ~isempty(curPtsonLine)
                [curSubSP_r,curSubSP_c]=ind2sub(size(bw),curPtsonLine);
                distance2Ct=sqrt( (C_x-curSubSP_c).^2+(C_y-curSubSP_r).^2);
                [minval,minind]=min(abs(curdist_s-distance2Ct));                
                AllSP_s(i)=curPtsonLine(minind);
            else
                AllSP_s(i)=AllSP(i);
%                 disp('find small problem here');
            end
        else
            AllSP_s(i)=AllSP(i);
        end       
    end
    AllSP=AllSP_s;
    
    if shown
        figure(11);imshow(im,'InitialMagnification','fit');hold on;
        for i=1:length(AllSP)
            [curSubSP_r,curSubSP_c]=ind2sub(size(bw),AllSP_s);
            curSubSP_r=[curSubSP_r curSubSP_r(1)];
            curSubSP_c=[curSubSP_c curSubSP_c(1)];            
            plot(curSubSP_c,curSubSP_r,'y');            
        end
        hold off;
    end
end
end