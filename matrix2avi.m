function matrix2avi(mtrx,varargin);

%Creates a avifile from a 3-D matrix and specified parameters
%
%mtrx is the 3-D matrix
%outfile is the *.avi file that you specify (in nothing is specified, a pop-up window appears
%varargin takes in parameters and values for the avifile in the format
%
% matrix2avi(mtrx,'Parameter1',Value_of_parameter1,'Parameter2',Value_of_parameter 2, ...)
%
%Valid parameters are
%
%   file:   Save into avi file (default is temp.avi in the current directory)
%
%   map :   The colormap to be used (default is jet)
%
%   fps :   Frames per second (default is 25)
%
%   tm_dm:  Time's Dimension in the 3-D matrix (default is 3)
%
%   chi:    Maximum value to clip to in the avi
%
%   clo:    Minimum value to clip to in the avi
%
%   lp:     Number of times, the matrix is to be played and recorded
%
%   trnsps: Set to 1 if transpose of the image is necessary

%Suresh E Joel, 19 Feb 2003.


if nargin<2, disp('Too few arguements'); return; end;
if size(size(mtrx),2)~=3, disp('Matrix has to be 3-D'); return; end;
video=struct( ...
        'mtrx',mtrx, ...
        'file','temp.avi', ...
        'map',jet, ...
        'fps',25, ...
        'tm_dm',3, ...
        'chi',max(max(max(mtrx))), ...
        'clo',min(min(min(mtrx))), ...
        'lp',1, ...
        'trnsps',0 ...
        );
close(gcf);
clear mtrx outfile; %To free up memory

%set variable parameters
i=1;
if nargin>2,
    while i<length(varargin),      
        if isfield(video,varargin{i}),
            if strcmp(varargin{i},'map'),
                figure;
                video=setfield(video,'map',colormap(varargin{i+1}));
                close(gcf);
            else
                video=setfield(video,varargin{i},varargin{i+1});
            end;
            i=i+2;
        else
            disp(strcat('Wrong parameter passed',varargin{i}));
            return;
        end;
    end;
    switch(video.tm_dm)
    case 1,
        video.mtrx=permute(video.mtrx,[2 3 1]);
    case 2,
        video.mtrx=permute(video.mtrx,[3 1 2]);
    otherwise
    end;
end;

%Verify default file
if(strcmp(video.file,'temp.avi')),
    y=questdlg('Save in temp.avi?');
    if(strcmp(y,'No'))
        [filename pathname]=uiputfile('*.avi','Save file');
        file=strcat(pathname,filename);
        if(~strcmp(file(length(file-3):length(file)),'.avi'))
            file=strcat(file,'.avi');
        end
    elseif (strcmp(y,'Cancel'))
        disp('No file written');
        return;
    end;
end;

xn=size(video.mtrx,1);
yn=size(video.mtrx,2);

%Open aviobject and write to file
av=avifile(video.file,'compression','none','colormap',video.map,'FPS',video.fps);
for k=1:video.lp,
    hw=waitbar(0,'Wait till all frames are recorded');
    for i=1:size(video.mtrx,3),
        b=reshape(video.mtrx(:,:,i),xn,yn);
        if video.trnsps, b=b'; end;
        c=find(b>video.chi); if ~isempty(c),for k=1:length(c),b(c(k))=video.chi; end;end;
        c=find(b<video.clo); if ~isempty(c),for k=1:length(c),b(c(k))=video.clo; end;end;
        b=fix(((b+abs(video.clo))/(video.chi-video.clo))*(length(video.map)-1)+1);
        av=addframe(av,im2frame(b,video.map));
        waitbar(i/size(video.mtrx,3),hw);
    end;
    close(hw);    
end;
av=close(av);