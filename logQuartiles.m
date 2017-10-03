function quartileidxs = logQuartiles(val)
% semilogxhist - generate histogram with M bars and log-scale x axis
% Downloaded from http://stackoverflow.com/questions/6812899/how-to-plot-hist-with-log-scale
% Modified by Brendon Watson 6/2014
%

M = 4;%num bins... quartiles
vmin=min(val(val>0));%minimum non-zero value, so logs don't get messed up
vmax=max(val);

edges=vmin*(vmax/vmin).^([0:M]/M);%log arraying of edges
[counts,quartileidxs] = histc(val,edges); 

quartileidxs(quartileidxs==0) = 1;%all zeros go into bin 1

% if size(counts,2)==1
%     counts=counts';
% end

% % x=edges(sort([1:M 1:M]));
% % y=[0 count(sort([1:M-1 1:M-1])) 0];
% % outline only: semilogx(x, y, '-');
% % plot(x, y, '-'); 
% % fill(x, y, 'b');


% for a = 1:M
%     if a == 1;
%         patches(a) = patch([0 0 edges(a) edges(a)],[0 counts(a) counts(a) 0],'b','EdgeColor','k');
%     else
%         patches(a) = patch([edges(a-1) edges(a-1) edges(a) edges(a)],[0 counts(a) counts(a) 0],'b','EdgeColor','k');
%     end        
% end
% set(gca,'XScale','log');

% edges = edges(1:end-1);%get rid of extra bin added at beginning
% counts = counts(1:end-1);% "

