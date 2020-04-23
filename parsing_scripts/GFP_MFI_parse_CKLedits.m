clear
close all

%For Macs, remember to use forward slashes (/) in directory names.

% specify data source directory. To read files from the current working 
% directory, use '.'
directory = '.'; 

% specify output directory; in this case, a subdirectory of `directory` 
% named `GFP` e.g.,
outdirname = fullfile(directory, 'GFP');

% how are your `.mat` files named? This should be a regular expression
% containing at least `channel` and `image` tokens/capture groups. See
% documentation for the `regexp` function. 
% this particular pattern expects files like this:
% Ph-some_long_experiment_name.mat
MAT_FILE_PATTERN = '^(?<channel>\w+)-(?<image>.*)\.(?<ext>\w+)';

% how are your `.tif` image files named? should be a format specifier
% suitable for `sprintf`. `%1$s` is the channel, `%2$s` is the base image
% name. This particular pattern expects files like this:
% Ph-some_long_experiment_name.tif
TIF_FILE_PATTERN = '%1$s-%2$s.tif';

% how do you want your output files named? should be a format string
% suitable for `sprintf` with two placeholders: `%1$s` is the base image
% name, `%2$d` is the frame number
OUT_FILE_PATTERN = 'gfpmi_%1$s_%2$d.dat';

% which frames to export? if this set to the empty array `[]`, all frames
% with 1+ cells will be exported. 
frames_to_export = [];

%%

% create the output directory
mkdir(outdirname)

% find all `.mat` files in `directory`
% note: use of `fullfile` function normalizes for system-specific path
% separators (`/` for Mac/Linux, `\` for Windows) and accounts for trailing
% slashes in directory names, etc. 
matfileobj = dir(fullfile(directory,'*.mat'));
[matfilenames{1:length(matfileobj)}] = matfileobj(:).name;
numfiles=length(matfileobj);

% for each `.mat` file
for fnum=1:numfiles
    
    % check that this file matches the expected file name pattern, and 
    % extract the base image name. If this is not a Phase image, skip
    file_name_parts = regexp(matfilenames{fnum}, MAT_FILE_PATTERN, 'names');
    if (isempty(file_name_parts))
        fprintf('Skipping .mat file "%s" which does not fit the expected pattern\n', matfilenames{fnum})
        continue
    elseif (~strcmp(file_name_parts.channel,'Ph')) 
        fprintf('Skipping .mat file "%s" which does not appear to be a phase contrast image\n', matfilenames{fnum})
        continue
    else 
        fprintf('Processing .mat file "%s"...\n', matfilenames{fnum})
    end
    
    % build path to the GFP image
    %imgfilestr=['GFP-' file_name_parts.image '.tif'];
    matfilename = fullfile(directory, matfilenames{fnum});
    img_base_name = file_name_parts.image; 
    imgfilestr=sprintf(TIF_FILE_PATTERN, 'GFP', img_base_name);
    fprintf('- Loading GFP data from "%s"...\n',imgfilestr)
    imgfilename=fullfile(directory,imgfilestr);
    
    
    %load matfile 
    load(matfilename)
    
    % which frames should we export? 
    % number of timepoints in movie
    tspan=length(cellListN);
    if (isempty(frames_to_export)) 
        % export all frames that have >0 cells
        tframes = find(cellListN>0); 
    else 
        tframes = frames_to_export; 
    end
     
    %loop through frames
    for t=1:length(tframes)
        %load raw fluoresence image in
        img=double(imread(imgfilename,tframes(t)));
        
        %parse oufti file for segmentation
        allIDlist = cellList.cellId;
        numcells = max(cell2mat(allIDlist(~cellfun('isempty',allIDlist))));
        cellIDs = cellList.cellId;
        meshinfo = cellList.meshData;
        
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
                %areavec(cellnum)=0;
                gfpmfi(cellnum)=NaN;
            else
                
                borders=[flip(meshcoord(:,4)) flip(meshcoord(:,3)); meshcoord(:,2) meshcoord(:,1)];
                infrows=find(ismember(borders,[Inf Inf],'rows')==1);
                
                if isempty(infrows)==0
                    borders(infrows,:)=[];
                end
                %areavec(cellnum)=polyarea([meshcoord(:,1);flip(meshcoord(:,3))],[meshcoord(:,2);flip(meshcoord(:,4))]);
                
                %generate mask containing individual cell
                bactmask=poly2mask((double(borders(:,2))),(double(borders(:,1))),size(img,1),size(img,2));
                
                %calculate average fluoresence intensity of individual cell
                gfpmfi(cellnum)=struct2array(regionprops(bactmask,img,'MeanIntensity'));
                
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
        outfilename=fullfile(outdirname,sprintf(OUT_FILE_PATTERN, img_base_name, t)); 
        fprintf('- Writing output data to "%s"...\n',outfilename)
        %outfilename=fullfile(outdirname,['gfpmfi_' imgfilestr '.dat']);
        save(outfilename,'gfpmfi','-ascii')
    end
end

%%
%fit a Gaussian mixture model to the GFP intensity distribution (optional)
close all
mdl=fitgmdist(gfpmfi',2,'CovarianceType','diagonal');
binnum=ceil(sqrt(length(gfpmfi)));
xaxis=linspace(min(gfpmfi),max(gfpmfi),binnum);
figure(2)
xlabel('Mean Fluorescence Intensity') %minus background
ylabel('Cumulative Probability Function')
title('CDF evaluated at each data point')
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

%Plotting histograms of GFP MFIs
figure(3)
[n, c]=ecdfhist(f,x,binnum);
hold on;
%bar(c,n)
plot(c,n,'ok')
nfit=ecdfhist(cdfsum,x,binnum);
plot(c,nfit,'-g')
xlabel('Mean Fluorescence Intensity')
ylabel('Probability Density')
pdf1=ecdfhist(cdf1,xax,binnum);
pdf2=ecdfhist(cdf2,xax,binnum);

%Plotting Gaussian mixed model fits
plot(c,pdf1,'--b')
plot(c,pdf2,'--r')
legstr0='EM Fit';
legstr1=['c_{1}=' num2str(mdl.ComponentProportion(smallcomp)) ' \mu_{1}=' num2str(mdl.mu(smallcomp)) ' \sigma_{1}=' num2str(sqrt(mdl.Sigma(smallcomp))) ];
legstr2=['c_{2}=' num2str(mdl.ComponentProportion(largecomp)) ' \mu_{2}=' num2str(mdl.mu(largecomp)) ' \sigma_{2}=' num2str(sqrt(mdl.Sigma(largecomp))) ];
legend({'combined data',legstr0,legstr1,legstr2})
title(['PDF evaluated with N=' num2str(binnum) 'bins'])
%%
