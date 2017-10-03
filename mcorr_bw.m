function varargout= mcorr_bw(varargin)
%MCORR_BW Multi-plot of all correlations between columns of a matrix
% Variant of mcorr, documentation is below.  Input is now a series of
% vectors rather than a matrix, and resultantly I'm able to plot the names
% of each input vector next to the related axes.  I also plot the linear
% fit line.
% Brendon Watson 2014
%
% MCORR (X) plots correlations between all possible combinations of the
% columns of array X, in a single figure. If the first argument is the 
% name of a file in the current directory, mcorr reads 
% it (including variable names in the first row) as X. Otherwise it is
% assumed that the first argument is an array, in which case consecutive 
% numbers will be used as variable names. If there is a second (numeric) 
% argument, mcorr will plot only the columns indicated in the second argument.
% The frequency distribution of each variable is plotted in the main
% diagonal with HIST using n/20 bins.
% It is unpractical to try to plot more than, say, 8 variables, since each
% individual plot becomes too small.
%
% OUTPUT= MCORR (X,'sig',ALPHA) calculates the correlation coefficient between
% each pair of columns and, if the correlation is sgnificant at the ALPHA
% level, points are plotted in red and the column numbers, corr. coefficient
% and p-value are returned in the OUTPUT array. Requires Statistical Toolbox
% 
% OUTPUT= MCORR (X,...,'corr',CORTYPE) specifies the type of correlation,
% allowed values are 'Pearson' (default), 'Spearman' and 'Kendall'
%
% EXAMPLES:
% mcorr ('myfile')
% mcorr ('myfile',[1:5])
% mcorr (X,[3:6])
% output= mcorr (X,'sig',0.05)
%
% Last modified: Jun. 2012
% BW: plotlogmode can be linear(default), logx, logy or loglog


if nargin < 1 
    error ('Need at least 1 argument: name of file or array')
end
if ~exist('plotlogmode','var')
    plotlogmode = 'linear';
end

color= [0 0 0];
output= [];
alpha= 0.05;
plotit= 'on';

datatoplot = {};
varnames = {};

for a = 1:length(varargin);
    if isstr(varargin{a})
        switch lower(varargin{a})
            case 'plotlogmode'%these must be at end
                plotlogmode = varargin{a+1};
        end
    else
        datatoplot{end+1} = varargin{a};
        varnames{end+1} = inputname(a);
    end
end


% corrtype= 'Pearson';

% % Its a file => Read. Else is a matrix
% if ischar(X) && ~isempty(dir(X)) 
% 	[A,varnames]= tblread(X);
% else
% 	A= X;
%  	varnames= num2cell(1:size(A,2));
% end

if length(datatoplot) == 1;
    A = datatoplot{1};
else
    for a = 1:length(datatoplot);
        A(a,:) = datatoplot{a}(:);
    end
end
A = A';
[nrow,ncol]= size(A);

% Function
for j= 1:ncol
    for k= 1:j
        if strcmp(plotit,'on')
            subplot(ncol,ncol,(j-1)*ncol+k);
            % Correlation plot or histogram,plotlogmode
            if j ~= k
                ScatterWStats(A(:,k),A(:,j),plotlogmode,'min',plotlogmode);
            else
                switch lower(plotlogmode)
                    case {'linear','logy'}  
                        hist(A(:,j),max([nrow/10,20])); % N/20 bins
                    case {'logx','loglog'}  
                        semilogxhist(A(:,j),max([nrow/10,10])); % N/20 bins
                        set(gca,'XScale','log')
                end
                axis tight
            end
            if k == 1 ylabel(varnames(j),'FontSize',7); end
            if j == ncol xlabel(varnames(k),'FontSize',7); end
        end


    end
end

varargout{1}= output;





