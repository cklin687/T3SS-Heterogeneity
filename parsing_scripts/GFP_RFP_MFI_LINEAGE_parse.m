%This script parses .mat files containing lineage data and two fluorescent
%channels, saves the results as a .csv where each row is a timepoint and
%columns contain fluoresence and lineage information.
%and prints a plot of the lineage for each colony. It also contains two in-line
%functions, trimtreeplot and trimtreelayout, which are required for
%plotting the trees. 

clear
close all
%specify input and output directories
directory='~\example_data_and_outputs\example_outputs\';
outdirname='~\example_data_and_outputs\example_outputs\lineage\';  
%Remember for Macs that slashes need to be backslashes (/).

%Load Oufti output .mat files from directory
matfileobj=dir([directory '*.mat']);
[matfilenames{1:length(matfileobj)}]=matfileobj(:).name;
numfiles=length(matfileobj);

%loop through colony files 
for fnum=1:numfiles
matfilename=[directory matfilenames{fnum}];
imgfilestr=['GFP-' matfilenames{fnum}(1:end-7) '.tif'];
rfpimgfilestr=['RFP-' matfilenames{fnum}(1:end-7) '.tif'];
imgfilename=[directory imgfilestr];
rfpfilename=[directory rfpimgfilestr];

load(matfilename);
tspan=length(cellListN);

%Uncomment the "tframes=input()" command if you want to analyze specific frames of the lineage tree. 
%Can do a range (Ex: [1:27]) or specific frames (Ex: [1,5,9,13,17,27,33]).
%tframes=input(['time frame? input number from 1 to ' num2str(tspan) '; if 1 frame just enter 1'])
tframes=find(cellListN>0);

%loop through timepoints 
for t=1:length(tframes)
img=double(imread(imgfilename,tframes(t)));    
rfpimg=double(imread(rfpfilename,tframes(t)));    
    allIDlist=cellList.cellId;

%Load CellIDs and segmentation mesh data.
numcells=max(cell2mat(allIDlist(~cellfun('isempty',allIDlist))));
cellIDs=cellList.cellId{t};
meshinfo=cellList.meshData{t};

%Define GFP and RFP MFIs
wholemask=zeros(size(img));
gfpmfi=zeros(1,cellListN(tframes(t)));
rfpmfi=zeros(1,cellListN(tframes(t)));
ancestorlist=zeros(1,cellListN(tframes(t)));
leafstatus=zeros(1,cellListN(tframes(t)));

%loop through cells 
for cellnum=1:cellListN(tframes(t)) 
    meshcoord=meshinfo{cellnum}.mesh;
    if meshcoord==0
        areavec(cellnum)=0;
        gfpmfi(cellnum)=NaN;
        rfpmfi(cellnum)=NaN;
    else
        
        borders=[flip(meshcoord(:,4)) flip(meshcoord(:,3)); meshcoord(:,2) meshcoord(:,1)];
   infrows=find(ismember(borders,[Inf Inf],'rows')==1);
   if isempty(infrows)==0
       borders(infrows,:)=[];
   end

    bactmask=poly2mask((double(borders(:,2))),(double(borders(:,1))),size(img,1),size(img,2));
    gfpmfi(cellnum)=mean(img(bactmask==1));     
    rfpmfi(cellnum)=mean(rfpimg(bactmask==1));
    wholemask=wholemask+cellnum*bactmask;
    
    birthframe(cellIDs(cellnum))=meshinfo{cellnum}.birthframe;
    polarity(cellIDs(cellnum))=meshinfo{cellnum}.polarity;
    ancestors{cellIDs(cellnum)}=meshinfo{cellnum}.ancestors;
    descendants{cellIDs(cellnum)}=meshinfo{cellnum}.descendants;
    
    	if isempty(ancestors{cellIDs(cellnum)})==0
		ancestorlist(cellnum)=max(ancestors{cellIDs(cellnum)});
   	 else
        	ancestorlist(cellnum)=NaN;
    	end
    end
end

%Calculating background fluorescence values
wholemaskBWbin=imdilate(wholemask>0,strel('disk',10));
backmask=imcomplement(wholemaskBWbin).*img;
backval=mean(img(backmask>0));
rfpbackval=mean(rfpimg(backmask>0));

%Subtract background fluorescence values
gfpmfi=gfpmfi-backval;
rfpmfi=rfpmfi-rfpbackval;

gfpmficell{t}=gfpmfi;
rfpmficell{t}=rfpmfi;
celllists{t}=double([cellList.cellId{tframes(t)}]);
ancestorscell{t}=ancestorlist;
timelabels{t}=tframes(t)*ones(size(gfpmfi));
end

%Define nodes for lineage tree
nodes=zeros(1,numcells);
A=find(cellfun('isempty',ancestors)==0);
for i=1:length(A)
    nodes(A(i))=max(ancestors{A(i)});
end

%Functions for Trimtreeplot and Trimtreelayout defined below; may need MATLAB 2018 and above to work
figure(fnum); hold on
title(matfilenames{fnum})
trimtreeplot(nodes);
[x y h s]=trimtreelayout(nodes);
for i = 1 :numcells
text(x(i),y(i),num2str(i),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
end

%Save GFP/RFP MFIs and define variables needed
outfilename=[outdirname matfilenames{fnum} '_twocolor_lineage.csv'];
allts=cell2mat(timelabels)';
allcellids=cell2mat(celllists)';
allgfps=cell2mat(gfpmficell)';
allrfps=cell2mat(rfpmficell)';
alllastancestors=cell2mat(ancestorscell)';

[count cellnums]=hist(nodes,unique(nodes));

allIDlist(cellfun('isempty', allIDlist))=[];
cellIDlist=unique(cell2mat(allIDlist));

leafIDs=setdiff(cellIDlist,cellnums);
leafnum=length(leafIDs);
leafstatus=ismember(allcellids,leafIDs);

%Tint=NaN(leafnum,tspan);
output=table(allts,allcellids,allgfps,allrfps,alllastancestors,leafstatus,'VariableNames',{'tframe','Cell','GFPint','RFPint','Last_Ancestor','LeafStatus'});
writetable(output,outfilename);
fnum
end

%%
%Plot intensity of GFP along a branch of the tree
[count cellnum]=hist(newtree,unique(newtree))
allIDlist(cellfun('isempty', allIDlist))=[];
cellIDlist=unique(cell2mat(allIDlist));
leafIDs=setdiff(cellIDlist,cellnum);
leafnum=length(leafIDs);

figure(1)
hold on
for ii=1:leafnum% =find(leafIDs==677)% 1 :leafnum
    paths=cell2mat(new_ancestors(leafIDs(ii)));
    if isempty(paths)==1
        fullpath=[leafIDs(ii)];
    else
        fullpath=[paths leafIDs(ii)];
    end
   % tmp=%mean_ints_one(fullpath,:);

    gfpleaftrace=allgfps(find(ismember(allcellids,fullpath)==1))
    tleaftrace=allts(find(ismember(allcellids,fullpath)==1))
    figure(1)
    hold on
    plot(tleaftrace,gfpleaftrace,'--')
end

xlabel('Time (Frame)') 
%Need to multiply tleaftrace by 15min/frame if want to plot in minutes
ylabel('Mean Fluorescence Intensity (AU)')


%%
%Functions defined below.
function trimtreeplot(p,c,d)
% TRIMTREEPLOT A modified standard function TREEPLOT. Plots the
%   leaves differently from TREEPLOT. They appear in their
%   respective layer instead of the deepest layer so that the tree
%   appears as "trimmed". The format is the same as for TREEPLOT,
%   see below. 
%
% trimtreeplot(p [,c,d])
%
% TREEPLOT Plot picture of tree.
%   TREEPLOT(p) plots a picture of a tree given a vector of
%   parent pointers, with p(i) == 0 for a root. 
%
%   TREEPLOT(P,nodeSpec,edgeSpec) allows optional parameters nodeSpec
%   and edgeSpec to set the node or edge color, marker, and linestyle.
%   Use '' to omit one or both.
%
%   See also ETREE, TREELAYOUT, ETREEPLOT.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.8 $  $Date: 1997/11/21 23:44:59 $
%
%   Modified by RG 01-Nov-1 to use trimtreelayout instead of
%   TREELAYOUT 

% RG use another treelayout 
%[x,y,h]=treelayout(p);
[x,y,h]=trimtreelayout(p);

f = find(p~=0);
pp = p(f);
X = [x(f); x(pp); repmat(NaN,size(f))];
Y = [y(f); y(pp); repmat(NaN,size(f))];
X = X(:);
Y = Y(:);

if nargin == 1,
    n = length(p);
    if n < 500,
        plot (x, y, 'r.', X, Y, 'r-');
    else,
        plot (X, Y, 'r-');
    end;
else,
    [ignore, clen] = size(c);
    if nargin < 3, 
        if clen > 1, 
            d = [c(1:clen-1) '.']; 
        else,
            d = 'r.';
        end;
    end;
    [ignore, dlen] = size(d);
    if clen>0 & dlen>0
        plot (x, y, c, X, Y, d);
    elseif clen>0,
        plot (x, y, c);
    elseif dlen>0,
        plot (X, Y, d);
    else
    end;
end;
xlabel(['height = ' int2str(h)]);
axis([0 1 0 1]);

end


function [x,y,h,s] = trimtreelayout(parent,post)
%   TRIMTREELAYOUT A modified standard function TREELAYOUT. 
%     Produces different heights for the leaves. They appear in
%     their respective layer instead of the deepest layer. The
%     format is the same as for TREELAYOUT, see below.
%
%   [x,y,h,s] = trimtreelayout(parent,post)
% 
%
%   TREELAYOUT Lay out tree or forest.
%   [x,y,h,s] = treelayout(parent,post)
%       parent is the vector of parent pointers, with 0 for a root.
%       post is a postorder permutation on the tree nodes.
%       (If post is omitted we compute it here.)
%       x and y are vectors of coordinates in the unit square at which 
%       to lay out the nodes of the tree to make a nice picture.
%       Optionally, h is the height of the tree and s is the 
%       number of vertices in the top-level separator.
%
%   See also ETREE, TREEPLOT, ETREEPLOT, SYMBFACT.

%   Copyright 1984-2000 The MathWorks, Inc. 
%   $Revision: 5.11 $  $Date: 2000/06/08 20:18:46 $
% 
%   Modified by RG 01-Nov-1 to plot the leaves in the right layer
%

% This is based on the C code in sptrees.c by John Gilbert.
% Leaves are spaced evenly on the x axis, and internal
% nodes are centered over their descendant leaves with
% y coordinate proportional to height in the tree.

n = length(parent);

 pv = [];
 if (size(parent,1)>1), parent = parent(:)'; end
 if (nargin<2) & ~all(parent==0 | parent>(1:n))
     % This does not appear to be in the form generated by ETREE.
     if (any(parent>n) | any(parent<0) | any(parent~=floor(parent)) ...
	 | any(parent==[1:n]))
        error('Bad vector of parent pointers.');
     end
     [parent,pv] = fixparent(parent);
 end

if nargin < 2,

    % Create the adjacency matrix A of the given tree,
    % and get the postorder with another call to etree.

    j = find(parent);
    A = sparse (parent(j), j, 1, n, n);
    A = A + A' + speye(n,n);
    [ignore, post] = etree(A);
%    post
end;

% Add a dummy root node #n+1, and identify the leaves.

parent = rem(parent+n, n+1) + 1;  % change all 0s to n+1s
isaleaf = ones(1,n+1);
isaleaf(parent) = zeros(n,1);

% In postorder, compute heights and descendant leaf intervals.
% Space leaves evenly in x (in postorder).

xmin = n(1,ones(1,n+1)); % n+1 copies of n
xmax = zeros(1,n+1);
height = zeros(1,n+1);
nkids = zeros(1,n+1);
nleaves = 0;

for i = 1:n,
    node = post(i);
    if isaleaf(node),
        nleaves = nleaves+1;
        xmin(node) = nleaves;
        xmax(node) = nleaves;
    end;
    dad = parent(node);
    %RG
%    height(dad) = max (height(dad), height(node)+1);
    xmin(dad)   = min (xmin(dad),   xmin(node));
    xmax(dad)   = max (xmax(dad),   xmax(node));
    nkids(dad)  = nkids(dad)+1;
end;

% RG compute heights
% traverse the tree from the root downwards in a layer-manner
lay_ind = n+1;
par_ind = n+1;
while(1)
  lay_ind = find(ismember(parent, lay_ind));
  if isempty(lay_ind)
    break;
  end  
  par_ind = [par_ind lay_ind];
  height(par_ind) = height(par_ind)+1;
end
height(1:n) = height(1:n)-1;

% Compute coordinates, leaving a little space on all sides.

treeht = height(n+1) - 1;
deltax = 1/(nleaves+1);
deltay = 1/(treeht+2);
x = deltax * (xmin+xmax)/2;
y = deltay * (height+1);

% Omit the dummy node.

x = x(1:n);
y = y(1:n);

% Return the height and top separator size.


h = treeht;
s = n+1 - max(find(nkids~=1));

if ~isempty(pv)
   x(pv) = x;
   y(pv) = y;
end

% ----------------------------
function [a,pv] = fixparent(parent)
%FIXPARENT  Fix order of parent vector
%   [A,PV] = FIXPARENT(B) takes a vector of parent nodes for an
%   elimination tree, and re-orders it to produce an equivalent vector
%   A in which parent nodes are always higher-numbered than child
%   nodes.  If B is an elimination tree produced by the TREE
%   functions, this step will not be necessary.  PV is a
%   permutation vector, so that A = B(PV);

n = length(parent);
a = parent;
a(a==0) = n+1;
pv = 1:n;

niter = 0;
while(1)
   k = find(a<(1:n));
   if isempty(k), break; end
   k = k(1);
   j = a(k);
   
   % Put node k before its parent node j
   a  = [ a(1:j-1)  a(k)  a(j:k-1)  a(k+1:end)]; 
   pv = [pv(1:j-1) pv(k) pv(j:k-1) pv(k+1:end)]; 
   t = (a >= j & a < k);
   a(a==k) = j;
   a(t) = a(t) + 1;

   niter = niter+1;
   if (niter>n*(n-1)/2), error('Bad vector of parent pointers.'); end
end

a(a>n) = 0;

end
end
