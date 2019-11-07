% addpath(genpath('C:\Nutstore\Nutstore\PathImAnalysis_Program\Program\RaidalLines\'));
% addpath(genpath('C:\Nutstore\Nutstore\PathImAnalysis_Program\Program\Miscellaneous\'));% for other dependencies

im=imread('your image here');
maskConfLHR=ones(size(im,1),size(im,2));
bw_nuclei=;%your nuclei mask here, binary map
TAreaRatio=.2% the ratio between the white region and the nuclei region
TsmalNucleiArea=60% the threshold for small nuclei, unit in pixel
debug=0; %show the imtermediate step or not
bwM=LDetectMelanocytes_RLS(im_ConfLHR,maskConfLHR,bw_nuclei,TAreaRatio,TsmalNucleiArea,debug);