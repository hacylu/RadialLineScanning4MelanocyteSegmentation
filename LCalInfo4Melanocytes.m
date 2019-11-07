%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is for calculating the potential useful info for melanocyte detection
% Input:
%   -ROI_GC,ROI_bw the green channel image and binary mask of the nuclei respectively
%   -AllSP the point index list that specifiy the supporting region(SR) for all the nuclei.
%
% Output:
%   -AllareaofSR:  the area of the SR
%   -AllmeanIntensityofSR:  the mean intensity of the SR
%   -AllContrast:  the contrast of the SR with respect to the nuclei region,
%                   i.e., AllmeanIntensityofSR(i)/mean(ColorinNuclei);
%   -AllAreaRatio: the ratio of the area of SR and the area of the nuclei
%                   region.
%   -AllmuDiff: the mean differences (mu_1-mu_2) of the whole SR+nuclei regions by the
%               two GMM.



% (c) Edited by Cheng Lu,
% Deptment of Eletrical and Computer Engineering,
% University of Alberta, Canada.  Aug, 2011
% If you have any problem feel free to contact me.
% Please address questions or comments to: hacylu@yahoo.com

% Terms of use: You are free to copy,
% distribute, display, and use this work, under the following
% conditions. (1) You must give the original authors credit. (2) You may
% not use or redistribute this work for commercial purposes. (3) You may
% not alter, transform, or build upon this work. (4) For any reuse or
% distribution, you must make clear to others the license terms of this
% work. (5) Any of these conditions can be waived if you get permission
% from the authors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [AllareaofSR,AllmeanIntensityofSR,AllContrast,...
    AllAreaRatio,AllmuDiff,AllNucleiArea]=LCalInfo4Melanocytes(ROI_GC,ROI_bw,AllSP,shown)

cc=bwconncomp(ROI_bw);
stats=regionprops(cc,'BoundingBox');
imsize=size(ROI_bw);

for i=1:cc.NumObjects
    %     curSP=AllSP{i};
    [curSP_r,curSP_c]=ind2sub(imsize,AllSP{i});
    %%% the current object turn it to binary mask
    curbw4SP=poly2mask(curSP_c,curSP_r,imsize(1),imsize(2));
    if shown
        show(curbw4SP,1);
    end
    
    % find the support regions' ind
    curbw4SPInd=find(curbw4SP==1);
    curbwInd=cc.PixelIdxList{i};
    curbwIndDiff=setdiff(curbw4SPInd,curbwInd);
    
    ColorIncurSR=ROI_GC(curbwIndDiff);
       
    AllmeanIntensityofSR(i)=mean(ColorIncurSR);
    AllareaofSR(i)=length(curbwIndDiff);
    AllNucleiArea(i)=length(curbwInd);
    %% for cal the contrast in SR and nuclei, the higher the more likely the melanocytes
    ColorinNuclei=ROI_GC(curbwInd);
    AllContrast(i)=AllmeanIntensityofSR(i)/mean(ColorinNuclei);
    %% for cal the ratio of the area of SR and the area of the nuclei
    AllAreaRatio(i)=AllareaofSR(i)/length(curbwInd);
    %% the mean differences (mu_1-mu_2) of the whole SR+nuclei regions by the two GMM.
    % current bounding box
%     curBBx=round(stats(i).BoundingBox);
%     % current image
%     curIm=ROI_GC(curBBx(2):curBBx(2)+curBBx(4)-1,curBBx(1):curBBx(1)+curBBx(3)-1);
%     curbwIndSRPlusNuclei= union(curbw4SPInd,curbwInd);
%     ColorIncurSRPlusNuclei=ROI_GC(curbwIndSRPlusNuclei);
%     
%     muDiff=LAnalysis8MoG4RadialLine(ColorIncurSRPlusNuclei,curIm,0);
%     AllmuDiff(i)=muDiff;
AllmuDiff(i)=0;
end
end