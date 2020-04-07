%this script classfies colonies based on given criteria output from control
%data and threshold_calculate.m



clear
%directory containing original mat files
directory ='F:\Dropbox\Christina_data\PA14 WT Pt-sfGFP in NTA\combined_4.4.19_4.14.19\'
matfileobj=dir([directory '*.mat']);
[matfilenames{1:length(matfileobj)}]=matfileobj(:).name;
numfiles=length(matfileobj)
%declare which frames should be analyzed
%tframes=[1,3,5,7,9,11,13,15,17,25];
tframes=25;
classinfo=xlsread('F:\Dropbox\Christina_scripts_for_github_upload\class_info.csv');
%loop through datfiles
for fnum=1:numfiles 
matfilename=[directory matfilenames{fnum}];  
for t=1:length(tframes)
%find output files based on original directory
outdirname=matfilename(1:end-4);   
outfilename=[outdirname '\gfpmfi_t' num2str(tframes(t)) '.dat'];

gfpmfi{t}=load(outfilename)
meanint(fnum,t)=nanmean(gfpmfi{t})
cellnum(fnum,t)=sum(isnan(gfpmfi{t})==0);
threshtmp=classinfo(find(classinfo(:,1)==tframes(t)),2);
if isempty(threshtmp)==1
threshtmp=classinfo(find(classinfo(:,1)==27),2);
end
thresh(t)=threshtmp;
numon=sum(sum(gfpmfi{t}>thresh(t)))
numoff=sum(sum(gfpmfi{t}<thresh(t)))
fon(fnum,t)=numon./(numon+numoff);
end
figure(10)
[F x]=ecdf(gfpmfi{end})

%plot CDF of each colony
plot(x,F,['-.'],'LineWidth',2)
hold on
if fnum==1
   plot(thresh(end)*ones(size(F)),F,'k','LineWidth',2)
hold on
end
xlabel('Intensity')
ylabel('Empirical CDF')

end

%%
%plot distribution of fraction on 
figure(1)
histogram(fon(:,end),'BinMethod','sqrt','FaceColor','k','FaceAlpha',0.5)
xlabel('fraction on')
ylabel('Number colonies')
ylim([0 40])

title('t=6 hrs')
%%
%plot number on vs time 
figure(2)
hold on
for i=1:size(fon,1)
plot(tframes, fon(i,:),'.-')
xlabel('frame')
ylabel('f_{on}')
end
ylim([-.1 1.1])

