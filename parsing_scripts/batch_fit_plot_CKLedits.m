clear
close all
%directory='C:\Users\dswle_000\Dropbox\Christina_data\12.29.2018 PA14 vs PA01 vs PA14 ladS01 Liquid Media Analysis\PA14 WT Pt-sfGFP exsA-rbs-mTagRFPt H12 MinS_NTA\';
%directory='C:\Users\dswle_000\Dropbox\Christina_data\01.01.2019 PA14 vs PA14 ladS01 and Pu-sfGFP vs Pt-sfGFP Liquid Media Analysis\PA14 WT PexoT-sfGFP H8 NTA\';
directory='/Users/Christina/Documents/Ph.D 2015-2019/Lab Data/Microscopy/2.1.19 PA14 PA01 WT ladS Y1 Y2 check TP #2/Copy_PA14 WT 12H/GFP/';
%day=dir([directory 'Expt*' ]);
numdays=1;%length(day);
%gfpints={};
%[daydirs{1:numdays}]=day(:).name;
%for daynum=1:numdays
 %   daydirectory=[directory daydirs{daynum} '\'];
    datobj=dir([directory '*.dat']);
   numfiles=length(datobj);
   
   
   daynum = 1;
   [datfilenames{1:numfiles}]=datobj(:).name;
    for fnum=1:numfiles
        filename=[directory datfilenames{fnum}];
        gfpints{daynum}{fnum}=load(filename);
    end
%end

for daynum=1%:numdays
   batch_ints{daynum}=cell2mat(gfpints{daynum});
   gfpmfi=[batch_ints{daynum}]';
   mdl=fitgmdist(gfpmfi,2,'CovarianceType','diagonal');
   figure(daynum)
   binnum=ceil(sqrt(length(gfpmfi)));
xaxis=linspace(min(gfpmfi),max(gfpmfi),binnum);
figure(daynum)
%subplot(1,2,1)
%xlabel('intensity (minus background)')
%ylabel('Cumulative Probability Function')
%title('CDF evaluated at each data point')
%hold on
[f, x]=ecdf(gfpmfi);
cdfsum=cdf(mdl,x);
%plot(x,cdfsum,'-g')
%plot(x,f,'.k')
xax=linspace(min(x),max(x),100);
smallcomp=find(mdl.mu==min(mdl.mu));
largecomp=find(mdl.mu==max(mdl.mu));
cdf1=mdl.ComponentProportion(smallcomp)*normcdf(xax,mdl.mu(smallcomp),sqrt(mdl.Sigma(smallcomp)));
cdf2=mdl.ComponentProportion(largecomp)*normcdf(xax,mdl.mu(largecomp),sqrt(mdl.Sigma(largecomp)));
%plot(xax,cdf1,'--b')
%plot(xax,cdf2+mdl.ComponentProportion(smallcomp),'--r')
%legend({'GMM fit','data','CDF1','CDF2+w1'})

%subplot(1,2,2)
[n c]=ecdfhist(f,x,binnum);
hold on;
bar(c,n,'g');
%plot(c,n,'.k')
nfit=ecdfhist(cdfsum,x,binnum);
%plot(c,nfit,'-g')
xlabel('GFP Mean Fluorescence Intensity (AU)')
ylabel('Probability Density')
pdf1=ecdfhist(cdf1,xax,binnum);
pdf2=ecdfhist(cdf2,xax,binnum);
plot(c,pdf1, 'LineWidth',2,'Color','b')
plot(c,pdf2,'LineWidth',2,'Color','r')
%legstr0='EM Fit';
%legstr1=['c_{1}=' num2str(mdl.ComponentProportion(smallcomp)) ' \mu_{1}=' num2str(mdl.mu(smallcomp)) ' \sigma_{1}=' num2str(sqrt(mdl.Sigma(smallcomp))) ];
%legstr2=['c_{2}=' num2str(mdl.ComponentProportion(largecomp)) ' \mu_{2}=' num2str(mdl.mu(largecomp)) ' \sigma_{2}=' num2str(sqrt(mdl.Sigma(largecomp))) ];
legstr1=['%T3SS-OFF=' num2str(mdl.ComponentProportion(smallcomp)) ' \mu_{1}=' num2str(mdl.mu(smallcomp)) ' \sigma_{1}=' num2str(sqrt(mdl.Sigma(smallcomp))) ];
legstr2=['%T3SS-ON=' num2str(mdl.ComponentProportion(largecomp)) ' \mu_{2}=' num2str(mdl.mu(largecomp)) ' \sigma_{2}=' num2str(sqrt(mdl.Sigma(largecomp))) ];
%legend({'Histogram of RFP MFI'})
legend({'Histogram of GFP MFI',legstr1,legstr2})
%title(['PDF evaluated with N=' num2str(binnum) 'bins'])
   
   coff(daynum)=mdl.ComponentProportion(smallcomp);
   con(daynum)=mdl.ComponentProportion(largecomp);
   muon(daynum)=mdl.mu(largecomp);
   muoff(daynum)=mdl.mu(smallcomp);
   sigon(daynum)=mdl.Sigma(largecomp);
   sigoff(daynum)=mdl.Sigma(smallcomp);
   
   
   
end
    
outfilename=['gfpmfi_PA14WT_12H.dat'];
save(outfilename,'gfpmfi','-ascii')
