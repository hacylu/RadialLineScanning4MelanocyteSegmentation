% RSP=radial line sppurt points
function AllSP_NoOverlap=LResolveOverlap4RSP(ROI_GC,ROI_bw,AllSP,AllIrre,shown)

cc=bwconncomp(ROI_bw,8);
stats=regionprops(cc,'Centroid','Area');
imsize=size(ROI_bw);
allCtemp=[stats.Centroid];
allC(:,1)=allCtemp(1:2:end);allC(:,2)=allCtemp(2:2:end);

for i=1:cc.NumObjects
    %     curSP=AllSP{i};
    curIrre=AllIrre(i);
    [curSP_r,curSP_c]=ind2sub(imsize,AllSP{i});
    %%% the current object turn it to binary mask
    curbw4SP=poly2mask(curSP_c,curSP_r,imsize(1),imsize(2));
    if shown
        show(curbw4SP,1);
    end
    %%% search the possible overlap mask of other nearby nuclei
    curC=stats(i).Centroid;
    curalldist=sqrt(sum((repmat(curC,size(allC,1),1)-allC).^2,2));
    
    % the searching range for detection of overlap
    Trang4DO=ceil(5*sqrt(stats(i).Area/pi));
    curLRInd=find(curalldist<Trang4DO&curalldist>0);% the current local region index
    if shown
        disp(sprintf('object %d has %d neighbour objects within range of %d pixels',i,length(curLRInd),Trang4DO));
    end
    % check all the local area regions if they are overlap with current analyis region
    if ~isempty(curLRInd)
        for ii=1:length(curLRInd)
            [curLRSP_r,curLRSP_c]=ind2sub(imsize,AllSP{curLRInd(ii)});
            curLRIrre=AllIrre(curLRInd(ii));
            %%% turn it to binary mask
            curbw4LRSP=poly2mask(curLRSP_c,curLRSP_r,imsize(1),imsize(2));
            if shown
                show(curbw4LRSP,2);
            end
            %%% check if has overlap with the curbw4SP
            bwOver=curbw4LRSP&curbw4SP;
            if shown
                show(bwOver,3);
            end
            if sum(bwOver(:))~=0 % in the case has overlap
                %%%% the one with larger nuclei and SP takes all the overlap
                %% !! decide which one get the overlap can be a better criterion !!
                %%%% use the shape irregularity as a criterion
                if curIrre<curLRIrre
%                 %%%% use the SP region's mean intensity as a criterion
%                 %%%% which one close to 255 is chosen
%                 curLRSPbwDiff=xor(curbw4LRSP,ismember(labelmatrix(cc),curLRInd(ii)));
%                 meancurLRSP=mean(ROI_GC(curLRSPbwDiff));%show(curSPbwDiff,1);
%                 curSPbwDiff=xor(curbw4SP,ismember(labelmatrix(cc),i));
%                 meancurSP=mean(ROI_GC(curSPbwDiff));
%                 if meancurSP>meancurLRSP
                    %                 %%%% use the SP region area as a criterion
                    %                 curLRSPAreaDiff=(sum(curbw4LRSP(:))-stats(curLRInd(ii)).Area);
                    %                 curSPAreaDiff=(sum(curbw4SP(:))-stats(i).Area);
                    %                 if curSPAreaDiff>curLRSPAreaDiff
                    
                    %%%% use the nuclei area as a criterion
%                     if stats(i).Area>stats(curLRInd(ii)).Area
                    
                    curbw4LRSP=xor(curbw4LRSP,bwOver);%show(curbw4LRSP,7);
                    %%% take the major result then track boundary
                    cctemp=bwconncomp(curbw4LRSP,8);
                    stemp=regionprops(cctemp,'Area');
                    [maxval,maxind]=max([stemp.Area]);
                    curbw4LRSP=ismember(labelmatrix(cctemp), maxind);
                    %show(curbw4SP,7);
                    if sum(curbw4LRSP(:))~=0
                    curLRBnd=bwboundaries(curbw4LRSP);
                    curLRBnd=curLRBnd{1};
                    curLRBndInd=sub2ind(imsize,curLRBnd(:,1), curLRBnd(:,2));
                    AllSP{curLRInd(ii)}=curLRBndInd';
                    else
                        continue;
                    end
                else
                    curbw4SP=xor(curbw4SP,bwOver);
                    if shown
                        show(curbw4SP,4);
                    end
                end
                
            else% in the case has no overlap
                continue;
            end
        end
        %%% store the changed boundary
        
        %%% take the major result then track boundary
        cctemp=bwconncomp(curbw4SP);
        if cctemp.NumObjects~=0
            stemp=regionprops(cctemp,'Area');
            [maxval,maxind]=max([stemp.Area]);
            curbw4SP=ismember(labelmatrix(cctemp), maxind);
            %show(curbw4SP,7);
            curBnd=bwboundaries(curbw4SP);
            
            curBnd=curBnd{1};
            curBndInd=sub2ind(imsize,curBnd(:,1), curBnd(:,2));
            AllSP_NoOverlap{i}=curBndInd';
        else
            AllSP_NoOverlap{i}=AllSP{i};
        end
        
    else
        AllSP_NoOverlap{i}=AllSP{i};
    end
    
end
end
