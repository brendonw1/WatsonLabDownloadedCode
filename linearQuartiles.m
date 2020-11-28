function quartileidxs = semilogQuartiles(val)
% semilogxhist - generate histogram with M bars and log-scale x axis
% Downloaded from http://stackoverflow.com/questions/6812899/how-to-plot-hist-with-log-scale
% Modified by Brendon Watson 6/2014
%


M = 4;%num bins... quartiles
% vmin=min(val);
vmax=max(val);
% 
% edges=vmin*(vmax/vmin).^([0:M]/M);%log arraying of edges
edges = linspace(0,vmax,5);
[counts,quartileidxs] = histc(val,edges); 
