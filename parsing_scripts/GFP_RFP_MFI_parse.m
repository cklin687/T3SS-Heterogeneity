clear;
close all;

%specify input directory 
directory = '~\example_data_and_outputs\example_outputs\';
outdirectory='~\example_data_and_outputs\outputs\twocolor\';
%For Macs, remember to change slashes to backslashes (/).

matfileobj=dir([directory '*.mat']);
[matfilenames{1:length(matfileobj)}]=matfileobj(:).name;
%'F:\Dropbox\Christina_data\PA14 WT Pt-sfGFP exsA-rbs-mTagRFPt NTA TL s1 col1_edited_lin\';

%Load in segmentation mat files and fluorescence images. I named filenames with the channel first (GFP-/RFP-).
numfiles=length(matfileobj);
for fnum=1:numfiles
matfilename=[directory matfilenames{fnum}];
load(matfilename)
imgfilestr=['GFP-' matfilenames{fnum}(1:end-7) '.tif'];
imgfilename=[directory imgfilestr];
rfpimgfilestr=['RFP-' matfilenames{fnum}(1:end-7) '.tif'];
rfpimgfilename=[directory rfpimgfilestr];

%tspan=length(cellListN);
tframes=find(cellListN>0);
%Frames I used for timepoints: [1,5,9,13,17,27,33]
outdirname=matfilename(1:end-4);   
mkdir(outdirname)
for t=1:length(tframes)
img=double(imread(imgfilename,tframes(t)));    
rfpimg=double(imread(rfpimgfilename,tframes(t)));    
    allIDlist=cellList.cellId;
numcells=max(cell2mat(allIDlist(~cellfun('isempty',allIDlist))));
cellIDs=cellList.cellId;
meshinfo=cellList.meshData;
  
wholemask=zeros(size(img));
gfpmfi=zeros(1,cellListN(tframes(t)));
rfpmfi=zeros(1,cellListN(tframes(t)));

for cellnum=1:cellListN(tframes(t)) 
    meshcoord=meshinfo{tframes(t)}{cellnum}.mesh;
    if meshcoord==0
        areavec(cellnum)=0;
        gfpmfi(cellnum)=NaN;
        rfpmfi(cellnum)=NaN;
    else
        
        borders=[flip(meshcoord(:,4)) flip(meshcoord(:,3)); meshcoord(:,2) meshcoord(:,1)];
   infrows=find(ismember(isinf(borders),[1 1],'rows')==1);
   if isempty(infrows)==0
       borders(infrows,:)=[];
   end
 %Defining segmentation masks
 %areavec(cellnum)=polyarea([meshcoord(:,1);flip(meshcoord(:,3))],[meshcoord(:,2);flip(meshcoord(:,4))]);
    bactmask=poly2mask((double(borders(:,2))),(double(borders(:,1))),size(img,1),size(img,2));
    gfpmfi(cellnum)=mean(img(bactmask==1)); %struct2array(regionprops(bactmask,img,'MeanIntensity'));
    rfpmfi(cellnum)=mean(rfpimg(bactmask==1)); %struct2array(regionprops(bactmask,img,'MeanIntensity'));

    wholemask=wholemask+cellnum*bactmask;
    end
end

%Defining background fluorescence values
if t==1
wholemaskBWbin=imdilate(wholemask>0,strel('disk',10));
backmask=imcomplement(wholemaskBWbin);
backval=mean(img(backmask>0));
rfpbackval=mean(rfpimg(backmask>0));
end

%Subtracting background fluorescence
gfpmfi=gfpmfi-backval;
rfpmfi=rfpmfi-rfpbackval;

%Saving GFP and RFP MFI values in new folders.
gfpoutfilename=[outdirectory '/twocolor/' matfilenames{fnum}(1:end-4) 'gfpmfi_only' num2str(tframes(t)) '.dat'];
save(gfpoutfilename,'gfpmfi','-ascii')

rfpoutfilename=[outdirectory '/twocolor/' matfilenames{fnum}(1:end-4) 'rfpmfi_only' num2str(tframes(t)) '.dat'];
save(rfpoutfilename, 'rfpmfi','-ascii')
end
end
