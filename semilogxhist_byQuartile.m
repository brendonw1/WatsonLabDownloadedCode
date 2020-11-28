function [edges,counts,histidxs,patches] = semilogxhist_byQuartile(val,M)
% semilogxhist - generate histogram with M bars and log-scale x axis
% Downloaded from http://stackoverflow.com/questions/6812899/how-to-plot-hist-with-log-scale
% Modified by Brendon Watson 6/2014
%

if nargin<2; 
    M=min(30,sqrt(length(val))); 
end

vmin=min(val(val>0));%minimum non-zero value, so logs don't get messed up
vmax=max(val);
edges=vmin*(vmax/vmin).^([0:M]/M);
[counts,histidxs] = histc(val,edges); 

zerovals = histidxs==0;
histidxs(zerovals) = 1;
counts(1) = counts(1)+sum(zerovals);

if size(counts,2)==1
    counts=counts';
end


% % x=edges(sort([1:M 1:M]));
% % y=[0 count(sort([1:M-1 1:M-1])) 0];
% % outline only: semilogx(x, y, '-');
% % plot(x, y, '-'); 
% % fill(x, y, 'b');
QuadColors = GetQuadColors;

for a = 1:M
    thisq = ceil(a/M/.25);
    if a == 1;
        patches(a) = patch([0 0 edges(a) edges(a)],[0 counts(a) counts(a) 0],QuadColors(thisq,:),'EdgeColor','k');
    else
        patches(a) = patch([edges(a-1) edges(a-1) edges(a) edges(a)],[0 counts(a) counts(a) 0],QuadColors(thisq,:),'EdgeColor','k');
    end        
end
set(gca,'XScale','log');

edges = edges(2:end);%get rid of extra bin added at beginning
counts = counts(2:end);% "

