%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is for detect the melanocyte with given nuclei mask
% implement two main methods: RLS and LDED
% Input:
%   -im_ConfLHR  epidermis image
%   -maskConfLHR epidermis mask
%   -maskAllCells_ATLRRS  the presegmented candidate regions
%   -TAreaRatio: Threshold for the area ratio
%   -TsmalNucleiArea: Threshold for the small area
%   -debug: show intermedia result or not

% Output:
%   -bwM     the binary mask for detected melanocytes
%
% (c) Edited by Cheng Lu,
% Deptment of Eletrical and Computer Engineering,
% University of Alberta, Canada.  12th Aug, 2011
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

function bwM=LDetectMelanocytes_RLS(im_ConfLHR,maskConfLHR,maskAllCells_ATLRRS,TAreaRatio,TsmalNucleiArea,debug)


if ~exist('debug','var')
    debug=0;
end

if ~exist('TAreaRatio','var')
    TAreaRatio=.6;
end

if ~exist('TsmalNucleiArea','var')
    TsmalNucleiArea=80;
end

if sum(maskAllCells_ATLRRS(:))~=0
    %% Initial Radial line scanning for the supporting point
    
    GC=im_ConfLHR(:,:,2);
    ROI_bw=maskAllCells_ATLRRS;
    ROI_GC=GC;
    
    cc=bwconncomp(ROI_bw);
    stats=regionprops(cc,'Centroid','Area');
    MeanIntenInEpi=LgetMeanColorInEpiArea(im_ConfLHR,maskConfLHR);
    AllSP=[];
    AllIrre=[];
    flagRegularization=1;
    
    for i=1:cc.NumObjects
        AllSP{i}=LfindOutterSPV3(stats(i).Centroid,stats(i).Area,...
            ROI_GC,ROI_bw,maskAllCells_ATLRRS,'GaussianBlur',MeanIntenInEpi,flagRegularization,cc.PixelIdxList{i},0);
        AllIrre(i)=LcalIrregularity(stats(i).Centroid,AllSP{i},size(ROI_bw));
    end
    
%     [AllSP_Area,AllSP_MeanIntensity,AllSP_Constrast]=LCalInfo4Melanocytes(ROI_GC,ROI_bw,AllSP,0);
    [AllSP_Area]=LCalInfo4Melanocytes_AreaofSRonly(ROI_GC,ROI_bw,AllSP,0);
    
    % AllRatio_Irre_SRArea=AllIrre(i)*100/AllSP_Area(i);% original Irre/SRArea ratio is good
    
    if debug
        % % plot out
        % figure(32);imshow(ROI_GC,'InitialMagnification','fit');hold on;
        LshowMaskCountouronIM(ROI_bw,ROI_GC,32);hold on;
        for i=1:length(AllSP)
            [curSubSP_r,curSubSP_c]=ind2sub(size(ROI_bw),AllSP{i});
            curSubSP_r=[curSubSP_r curSubSP_r(1)];
            curSubSP_c=[curSubSP_c curSubSP_c(1)];
            plot(curSubSP_c,curSubSP_r,'y','Linewidth',2);
            
            %     text(stats(i).Centroid(1),stats(i).Centroid(2),num2str(AllIrre(i)*100/AllSP_Area(i),'%.2f'),'color','y');
            
        end
        hold off;
    end
    %% smooth it  using morphological operations
    % AllSP_s=LGetSmoothBnd(ROI_GC,ROI_bw,AllSP,0);
    % for i=1:cc.NumObjects
    %        AllIrre(i)=LcalIrregularity(stats(i).Centroid,AllSP_s{i},size(ROI_bw));
    % end
    % %%% plot out
    % figure(33);imshow(ROI_GC,'InitialMagnification','fit');hold on;
    % for i=1:length(AllSP_s)
    %     [curSubSP_r,curSubSP_c]=ind2sub(size(ROI_bw),AllSP_s{i});
    %     curSubSP_r=[curSubSP_r curSubSP_r(1)];
    %     curSubSP_c=[curSubSP_c curSubSP_c(1)];
    %
    %     plot(curSubSP_c,curSubSP_r,'y');
    %          text(stats(i).Centroid(1),stats(i).Centroid(2),num2str(AllIrre(i)*100/stats(i).Area,'%.2f'),'color','y');
    %
    % end
    % hold off;
    %% resolve the overlap
    % AllNucleiArea=[stats.Area];
    
    AllSP_NoOverlap=LResolveOverlap4RSP(ROI_GC,ROI_bw,AllSP,AllIrre*100./AllSP_Area,0);
    
    if debug
        %%% plot out
        % figure(34);imshow(ROI_GC,'InitialMagnification','fit');hold on;
        LshowMaskCountouronIM(ROI_bw,ROI_GC,34);hold on;
        title('Overlap Resovled');
        for i=1:length(AllSP_NoOverlap)
            [curSubSP_r,curSubSP_c]=ind2sub(size(ROI_bw),AllSP_NoOverlap{i});
            curSubSP_r=[curSubSP_r curSubSP_r(1)];
            curSubSP_c=[curSubSP_c curSubSP_c(1)];
            
            plot(curSubSP_c,curSubSP_r,'y','Linewidth',2);
            %     text(stats(i).Centroid(1),stats(i).Centroid(2),num2str(AllIrre(i)*100/AllSP_Area(i),'%.2f'),'color','y');
            
        end
        hold off;
    end
    %% Analysis the info in the SP
    [AllSP_Area,AllSP_MeanIntensity,AllSP_Contrast,...
        AllAreaRatio,AllmuDiff,AllNucleiArea]=LCalInfo4Melanocytes(ROI_GC,ROI_bw,AllSP_NoOverlap,0);
    
    %%% the thrshold should be determined carefully
    
    if debug
        %%% plot out
        LshowMaskCountouronIM(ROI_bw,ROI_GC,37);hold on;
        for i=1:length(AllSP_NoOverlap)
            if 1%AllAreaRatio(i) > TAreaRatio
                %AllSP_Area(i)>TDiffarea%&&AllSP_MeanIntensity(i)>TMeanInten...
                % &&AllSP_Constrast(i)>TContrast&&(AllIrre(i)*100/stats(i).Area)<TIrre
                [curSubSP_r,curSubSP_c]=ind2sub(size(ROI_bw),AllSP_NoOverlap{i});
                curSubSP_r=[curSubSP_r curSubSP_r(1)];
                curSubSP_c=[curSubSP_c curSubSP_c(1)];
                plot(curSubSP_c,curSubSP_r,'b');
                text(stats(i).Centroid(1),stats(i).Centroid(2),num2str( AllAreaRatio(i),'%.2f'),'color','y');
            end
            
        end
        hold off;
    end
    %% Filtering out the true malenocytes
    
    cc=bwconncomp(ROI_bw);
    
    idx = find(AllAreaRatio > TAreaRatio& AllNucleiArea>TsmalNucleiArea);% & AllSP_Contrast>TContrast);
    bwM = ismember(labelmatrix(cc), idx);
    
    if debug
        LshowMaskCountouronIM(bwM,im_ConfLHR,38);
    end
else
    bwM=maskAllCells_ATLRRS;
end