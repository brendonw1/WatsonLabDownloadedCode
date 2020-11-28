function varargout= mcorr_bw_ForCellIDs(CellIDs,varargin)
%MCORR_BW Multi-plot of all correlations between columns of a matrix
% Variant of mcorr, documentation is below.  Input is now a series of
% vectors rather than a matrix, and resultantly I'm able to plot the names
% of each input vector next to the related axes.  I also plot the linear
% fit line.
%  In the upper triangle of the figure now are log plots of percent changes
%  for each cell across each pair of condions (
%
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


if nargin < 2 error ('Need at least 2 argument: name of file or array'); end
color= [0 0 0];
output= [];
alpha= 0.05;
plotit= 'on';
corrtype= 'Pearson';

% % Its a file => Read. Else is a matrix
% if ischar(X) && ~isempty(dir(X)) 
% 	[A,varnames]= tblread(X);
% else
% 	A= X;
%  	varnames= num2cell(1:size(A,2));
% end

for a = 1:length(varargin);
    varnames{a} = inputname(a+1);
    A(a,:) = varargin{a}(:);
end

allE = CellIDs.EAll;
allI = CellIDs.IAll;
isE = zeros(1,size(A,2));
isE(allE) = 1;
dotcolors = [1 0 0 ; 0 1 0];

% % Optional arguments
% for j= 1:2:length(varargin)
% 	string= lower(varargin{j});
% 	switch string(1:min(3,length(string)))
% 		case 'col'
% 			columns= varargin{j+1};
% 			A= A(:,columns);
% 			varnames= varnames(varargin{j+1});
% 		case 'sig'
% 			pearson= 1;
% 			alpha= varargin{j+1};
% 		case 'plo'
% 			plotit= varargin{j+1};
% 		case 'cor'
% 			corrtype= varargin{j+1};
% 		otherwise
% 			error('MCORR Unknown parameter name');
% 	end
% end

A = A';
[nrow,ncol]= size(A);


% Function
for j= 1:ncol
    for k= 1:ncol
        
        if j<k%in upper part of matrix... do a log plot
            subplot(ncol,ncol,(j-1)*ncol+k);
            hold on
            for ix = 1:size(A,1)
                dotcolor = dotcolors(1+isE(ix),:);
                plot(ix,log10(A(ix,k)./A(ix,j)),'.','color',dotcolor)
            end
            plot([0 size(A,1)],[0 0],'k')
            axis tight
        else
            % Correlation
            if j ~= k
                [yfit,r,pval]= RegressAndFindR2(A(:,k),A(:,j));
                if pval < alpha
                    linecolor= [1 0 0];
                    output= [output;k,j,r,pval];
                else
                    linecolor= [.5 .5 .5]; 
                end
            end

            if strcmp(plotit,'on')
                subplot(ncol,ncol,(j-1)*ncol+k);
                hold on

                % Correlation plot or histogram
                if j ~= k
                    for ix = 1:size(A,1)
                        dotcolor = dotcolors(1+isE(ix),:);
                        plot(A(ix,k),A(ix,j),'.','MarkerSize',5,'MarkerEdgeColor',dotcolor);
                    end
                else
                    hist(A(:,j),nrow/20); % N/20 bins
                end

                % Colors and fonts
                h = findobj(gca,'Type','patch');
                set(h,'FaceColor',[.2 .4 1],'EdgeColor','none')
                set (gca,'FontSize',6);
                if k == 1 ylabel(varnames(j),'FontSize',7); end
                if j == ncol xlabel(varnames(k),'FontSize',7); end

                if j ~= k
                    hold on;%plot overlaid best fit line
                    xl = get(gca,'xlim');
                    yl = get(gca,'ylim');
                    plot(A(:,k),yfit,'color',linecolor)
                    text(xl(1)+0.1*diff(xl),yl(1)+0.8*diff(yl),['p=',num2str(pval)])
                end
            end

        end
    end
end

varargout{1}= output;





