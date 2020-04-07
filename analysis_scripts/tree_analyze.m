%this script parses oufti outputs containing lineage trees and performs
%change point analysis in two channels. 
clear
close all
directory='F:\Dropbox\Christina_data\PA14 WT Pt-sfGFP exsA-RFP in NTA\4.11.19 Lineage Tracings\';
xlsfileobj=dir([directory '*.csv']);
[xlsfilenames{1:length(xlsfileobj)}]=xlsfileobj(:).name;
numfiles=length(xlsfileobj);
%manually defined threshold, based on classification 
thresh=316.3844354;


Ts=[];
Gs=[];
Rs=[];
Ps=[];
tick=1;
Nleaves=0;
%
for fnum=1:numfiles
xlsfilename=[directory xlsfilenames{fnum}];

[num raw txt]=xlsread(xlsfilename);

times=num(:,1);
celllist=num(:,2);
gfplist=num(:,3);
rfplist=num(:,4);
ancestorlist=num(:,5);
leafstatus=num(:,6);
ontime=zeros(size(celllist));

cellids=unique(num(:,2));
for cellnum=1:length(cellids)
    indlist=find(num(:,2)==cellids(cellnum));
    timelist=num(indlist,1);
        gfpints=(num(indlist,3));
    rfpints=(num(indlist,4));
    finalint{fnum}(cellnum)=(gfpints(timelist==max(timelist)));
            finalrint{fnum}(cellnum)=(rfpints(timelist==max(timelist)));
        if max(timelist)<25 & min(timelist)>1
    lifetime{fnum}(cellnum)=length(indlist);

  oncat{fnum}(cellnum)=sum([finalint{fnum}>thresh])./length([finalint{fnum}]);
    flist{fnum}(cellnum)=fnum;
  cellIDlist{fnum}(cellnum)=cellids(cellnum);
  

        else 
            lifetime{fnum}(cellnum)=NaN;
  oncat{fnum}(cellnum)=NaN;
  flist{fnum}(cellnum)=NaN;
  cellIDlist{fnum}(cellnum)=NaN;
        end


end
leaflist=unique(celllist(leafstatus==1));
    Nleaves=Nleaves+length(leaflist);
for pathnum=1:length(leaflist)

    traj=leaflist(pathnum); 
    nextnode=mean(ancestorlist(celllist==leaflist(pathnum)));
    while isnan(nextnode)==0
traj=[traj  nextnode];
nextnode=mean(ancestorlist(celllist==nextnode));
    end
    trajinds=find(ismember(celllist,traj')==1);
     

gcoefs{fnum,pathnum}=        polyfit(times(trajinds),(log(gfplist(trajinds)./gfplist(trajinds(1)))),1)        ;
rcoefs{fnum,pathnum}=        polyfit(times(trajinds),(log(rfplist(trajinds)./rfplist(trajinds(1)))),1);

if(sum(isnan(gfplist(trajinds)))>0)
gipt{fnum,pathnum}=NaN;
ript{fnum,pathnum}=NaN;    

else
gipt{fnum,pathnum}=findchangepts(real(log(gfplist(trajinds)./gfplist(trajinds(1)))),'MaxNumChanges',1,'Statistic','linear'); %find gfpchangepoints
ript{fnum,pathnum}=findchangepts(real(log(rfplist(trajinds)./rfplist(trajinds(1)))),'MaxNumChanges',1,'Statistic','linear'); %find rfp changepoints
end
      


    ontime = find((gfplist(trajinds))>thresh); %use time 27 threshold
    if isempty(ontime)==0
    ontime(pathnum)=min(ontime);
    
    tbeforeon=times(trajinds(1:ontime(pathnum)))-ontime(pathnum);
    gfpbeforeon=gfplist(trajinds(1:ontime(pathnum)));
    rfpbeforeon=rfplist(trajinds(1:ontime(pathnum)));
    gfpbeforenorm=(gfpbeforeon-min(gfpbeforeon))./(max(gfpbeforeon)-min(gfpbeforeon));
    rfpbeforenorm=(rfpbeforeon-min(rfpbeforeon))./(max(rfpbeforeon)-min(rfpbeforeon));
    
    leafnums=tick.*ones(size(gfpbeforeon));
%    Ts=[Ts;times(trajinds(1:ontime(pathnum)))-ontime(pathnum)];

    Ts=[Ts;times(trajinds(1:ontime(pathnum)))];
Gs=[Gs;gfpbeforeon];
    Rs=[Rs;rfpbeforeon];
    Ps=[Ps;leafnums];
hold on



    end


tick=tick+1;


end
end
%%
all_lifetimes=cell2mat(lifetime);
all_ints=cell2mat(finalint);
all_rints=cell2mat(finalrint);
all_velo=cell2mat(velo);

allfnums=cell2mat(flist);
allcellnums=cell2mat(cellIDlist);


alloncat=cell2mat(oncat);
allgparams=cell2mat(reshape(gcoefs,[size(gcoefs,1)*size(gcoefs,2) 1]));
allrparams=cell2mat(reshape(rcoefs,[size(gcoefs,1)*size(gcoefs,2) 1]));

emptyIndex = cellfun('isempty', gipt);     % Find indices of empty cells
gipt(emptyIndex) = {NaN};                    % Fill empty cells with 0
emptyIndex = cellfun('isempty', ript);     % Find indices of empty cells
ript(emptyIndex) = {NaN};                    % Fill empty cells with 0

allript=cell2mat(ript);
allgipt=cell2mat(gipt);
figure()
histogram(allgipt-allript)
nanmean(allgipt(:)-allript(:))
