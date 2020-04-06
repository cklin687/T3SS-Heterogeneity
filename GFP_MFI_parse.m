clear


% specify matching oufti output(.mat) and fluoresence image file

matfilename='/Users/Christina/Documents/Ph.D 2015-2019/Lab Data/Microscopy/PA14 WT vs ladS01 vs PA01 Lineage Analysis/2.14.19 PA14 WT/Colony 45/PA14 WT Pt-sfGFP MinS_NTA s17 col45.mat';
imgfilename='/Users/Christina/Documents/Ph.D 2015-2019/Lab Data/Microscopy/PA14 WT vs ladS01 vs PA01 Lineage Analysis/2.14.19 PA14 WT/Colony 45/GFP-PA14 WT Pt-sfGFP MinS_NTA s17 col45.tif';

%load matfile

load(matfilename);
%number of timepoints in movie
tspan=length(cellListN);

%prompt desired frames for analysis, e.g., [1,5,9,13,17,27,33] (can also set 1:tspan)

tframes=input(['time frame? input number from 1 to ' num2str(tspan) '; if 1 frame just enter 1']);
%tframes=find(cellListN>0);  %can also set this


%make a directory for each image file in same place as mat file
outdirname=matfilename(1:end-4);   
mkdir(outdirname)


%loop through frames

for t=1:length(tframes)
%load raw fluoresence image in
img=double(imread(imgfilename,tframes(t)));    
%parse oufti file for segmentation

allIDlist=cellList.cellId;
numcells=max(cell2mat(allIDlist(~cellfun('isempty',allIDlist))));
cellIDs=cellList.cellId;
meshinfo=cellList.meshData;

%generate image mask based on oufti segmentation
wholemask=zeros(size(img));
%preallocate output matrix

gfpmfi=zeros(1,cellListN(tframes(t)));

%loop through each segmented cell


for cellnum=1:cellListN(tframes(t)) 
%parse oufti output for coordinates of cell
meshcoord=meshinfo{tframes(t)}{cellnum}.mesh;
%ignore empty coordinates which are output sometimes for some reason
    if meshcoord==0
        areavec(cellnum)=0;
        gfpmfi(cellnum)=NaN;
    else
       
   borders=[flip(meshcoord(:,4)) flip(meshcoord(:,3)); meshcoord(:,2) meshcoord(:,1)];
   infrows=find(ismember(borders,[Inf Inf],'rows')==1);
   
   if isempty(infrows)==0
       borders(infrows,:)=[];
   end
 %   areavec(cellnum)=polyarea([meshcoord(:,1);flip(meshcoord(:,3))],[meshcoord(:,2);flip(meshcoord(:,4))]);
%generate mask containing individual cell
 bactmask=poly2mask((double(borders(:,2))),(double(borders(:,1))),size(img,1),size(img,2));
 %calculate average fluoresence intensity of individual cell
    gfpmfi(cellnum)=struct2array(regionprops(bactmask,img,'ManIntensity'));
 %collate a whole-colony mask by accumulating individual cell mask together
    wholemask=wholemask+cellnum*bactmask;
    end
end
%dilate the mask containing all colonies and calculate the average
%intensity outside it to calculate the background fluoresence intensity
wholemaskBWbin=imdilate(wholemask>0,strel('disk',10));

backmask=imcomplement(wholemaskBWbin).*img;
%can alternatively replace this value with a manually input one
backval=mean(img(backmask>0));
gfpmfi=gfpmfi-backval;
%output a .dat file containing intensity values
outfilename=[outdirname '/gfpmfi_t' num2str(tframes(t)) 's17col45.dat'];
save(outfilename,'gfpmfi','-ascii')
end

%%
%fit a gaussian mixture model to the GFP intensity distribution (optional)
close all
mdl=fitgmdist(gfpmfi',2,'CovarianceType','diagonal');
binnum=ceil(sqrt(length(gfpmfi)));
xaxis=linspace(min(gfpmfi),max(gfpmfi),binnum);
figure(2)
xlabel('intensity (minus background)')
ylabel('Cumulative Probability Function')
title('CDF evalutated at each data point')
hold on
[f,x]=ecdf(gfpmfi);
cdfsum=cdf(mdl,x);
plot(x,cdfsum,'-g')
plot(x,f,'.k')
xax=linspace(min(x),max(x),100);
smallcomp=find(mdl.mu==min(mdl.mu));
largecomp=find(mdl.mu==max(mdl.mu));
cdf1=mdl.ComponentProportion(smallcomp)*normcdf(xax,mdl.mu(smallcomp),sqrt(mdl.Sigma(smallcomp)));

cdf2=mdl.ComponentProportion(largecomp)*normcdf(xax,mdl.mu(largecomp),sqrt(mdl.Sigma(largecomp)));
plot(xax,cdf1,'--b')
plot(xax,cdf2+mdl.ComponentProportion(smallcomp),'--r')
legend({'GMM fit','data','CDF1','CDF2+w1'})
figure(3)
[n, c]=ecdfhist(f,x,binnum);
hold on;
%bar(c,n)
plot(c,n,'ok')
nfit=ecdfhist(cdfsum,x,binnum);
plot(c,nfit,'-g')
xlabel('Intensity')
ylabel('probability density')
pdf1=ecdfhist(cdf1,xax,binnum);
pdf2=ecdfhist(cdf2,xax,binnum);
plot(c,pdf1,'--b')
plot(c,pdf2,'--r')
legstr0='EM Fit';
legstr1=['c_{1}=' num2str(mdl.ComponentProportion(smallcomp)) ' \mu_{1}=' num2str(mdl.mu(smallcomp)) ' \sigma_{1}=' num2str(sqrt(mdl.Sigma(smallcomp))) ];
legstr2=['c_{2}=' num2str(mdl.ComponentProportion(largecomp)) ' \mu_{2}=' num2str(mdl.mu(largecomp)) ' \sigma_{2}=' num2str(sqrt(mdl.Sigma(largecomp))) ];
legend({'combined data',legstr0,legstr1,legstr2})
title(['PDF evaluated with N=' num2str(binnum) 'bins'])
%%
