function rscstr = getrscstr(EmbeddedFonts,NeededFonts)
%GETRSCSTR   Form DSC headers DSC v3 compliant resource tags for the fonts
% '%%DocumentNeededResources' & '%%DocumentSuppliedResources:'

% Copyright 2013 Takeshi Ikuma
% History:
% rev. - : (04-25-2013) original release

rscstr = '%%DocumentNeededResources:';
if ~isempty(NeededFonts)
   rscstr = sprintf('%s font',rscstr);
   N = 127;
   for n = 1:numel(NeededFonts)
      rscstr_new = sprintf('%s %s',rscstr,NeededFonts{n});
      if numel(rscstr_new)>N % must change line
         N = numel(rscstr)+127;
         rscstr = sprintf('%s\n%%%%+ font %s',rscstr,NeededFonts{n});
      else
         rscstr = rscstr_new;
      end
   end
end

rscstr = sprintf('%s\n%%%%DocumentSuppliedResources:',rscstr);
if ~isempty(EmbeddedFonts)
   rscstr = sprintf('%s font',rscstr);
   N = 127;
   for n = 1:numel(EmbeddedFonts)
      rscstr_new = sprintf('%s %s',rscstr,EmbeddedFonts{n});
      if numel(rscstr_new)>N % must change line
         N = numel(rscstr)+127;
         rscstr = sprintf('%s\n%%%%+ font %s',rscstr,EmbeddedFonts{n});
      else
         rscstr = rscstr_new;
      end
   end
end
rscstr = sprintf('%s\n',rscstr);
