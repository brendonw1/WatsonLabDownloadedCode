function [centers,counts,minedge,histidxs,patches] = semilogxhist(val,M,plotting)
% semilogxhist - generate histogram with M bars and log-scale x axis.
% Plotting is a boolean to determine whether to plot patches or not.
% Default = 1(plottig on)
% Downloaded from http://stackoverflow.com/questions/6812899/how-to-plot-hist-with-log-scale
% Modified by Brendon Watson 6/2014
%

if nargin<2; 
    M=min(30,sqrt(length(val))); 
end
if ~exist('plotting','var')
    plotting = 1;
end

vmin=min(val(val>0));
vmax=max(val(val<Inf));
edges=vmin*(vmax/vmin).^(linspace(eps,M,M)/M);
minedge = edges(1);
[counts,histidxs] = histc(val,edges); 

if size(counts,2)==1
    counts=counts';
end

patches = [];
if plotting
    for a = 1:M
        if a == 1;
            te = edges(1)/(edges(2)/edges(1));
            patches(a) = patch([te te edges(a) edges(a)],[0 counts(a) counts(a) 0],'b','EdgeColor','k');
        else
            patches(a) = patch([edges(a-1) edges(a-1) edges(a) edges(a)],[0 counts(a) counts(a) 0],'b','EdgeColor','k');
        end        
    end
    % set(gca,'XScale','log');
end

centers = mean([edges(1:end-1);edges(2:end)],1);%get rid of extra bin added at beginning
counts = counts(2:end);% "

quartilebounds = 1;

