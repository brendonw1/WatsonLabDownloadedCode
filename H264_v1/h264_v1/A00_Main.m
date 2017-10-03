%% Clean all and Close all previuse compution
clc
clear all;
close all;
system_dependent('DirChangeHandleWarn', 'Never');
addpath(genpath('.'));
%% Encoding Video Paramiters
Frame_start = 1;            % I frame
Frame_end = 10;             % Following P frames
NumberOfFrames=Frame_start+Frame_end;
%% Read video and save it in 00OriginalVideo.mat
A01Read('xylophone.mpg',NumberOfFrames); % You have to do it once, then comment it
load 00OriginalVideo.mat;
%% Option: playing & resizing 
A01Play(mov,Obj);
A01Resize(mov,Frame_start,Frame_end,128,128) 
%% Encoding Inputs/paramiters
load 01ResizedFrames.mat;   % 1. Input video sequance (VideoSeq_Input)
QP = 27;                    % 2. Quality Praramter, QP values

%% Encoding Setup
global h w              % Hight and Width
block_size = 16;        % Macroblock size for P frames
ext = 0;                % switch for extended ME
[h,w,N] = size(VideoSeq_Input); % Extract input video sequance diminsion

%% Encoding Outputs intialization
Frames_PSNR = zeros(N,1);       % Initialize PSNR for each frame
Frames_BR =  zeros(N,1);        % Initialize BitRate for each frame
VideoSeq_Rec = zeros(h,w,N);    % Initialize Reconstructed video sequance
bitstream = '';                 % Initialize output bitstream
%% Encoding Processes
%% 1. Saving the header: Encoding Inputs/paramiters (1x byte each)
[bits] = header(h,w,QP,Frame_start,Frame_end);
bitstream = bits;
%% 2. Encode I-Frame
disp(['Encoding I-Frame: ',num2str(Frame_start)]);      % Dispaly
bitstream = [bitstream '1111'];     % Appending I-Frame header ('1111')
Seq(:,:,1) = double(VideoSeq_Input(:,:,1)); % Extracting I-Frame
[Seq_r(:,:,1),bits] = encode_i_frame(Seq(:,:,1),QP);    % Encoding
Frames_Rec(:,:,1)=Seq(:,:,1);   % Storing Rec. Frame
bitstream = [bitstream bits];   % Appending I-Frame bitstream ('11100...')
X(:,:,1) = Seq_r(:,:,1);        % X[1]: Refrence Frame (for P-Frames)
%% 3. Encoding P-Frams
for K = 2:Frame_end
    k = K-1;
    disp(['Encoding P Frame: ',num2str(K)]);            % Dispaly
    bitstream = [bitstream '0000'];     % Appending P-Frame header ('0000')
    Seq(:,:,2) = double(VideoSeq_Input(:,:,K)); % Extracting P-Frame
    X(:,:,2) = Seq(:,:,2); % X[2]: Rnter coded Frame (P-Frames)
    [Seq_r(:,:,2),bits] = encode_p_frame(X,QP,ext,block_size); % Encoding
    Frames_Rec(:,:,K)=Seq_r(:,:,2); % Storing Rec. Frame
    bitstream = [bitstream bits]; % Appending P-Frame bitstream ('1100...')
    X(:,:,1) = Seq_r(:,:,2); % X[1]: Refrence Frame (for the next P-Frames)
end

%% Storage
save('02BitStream.mat','bitstream');                % 1. Bitstream
save('02ReconstructedVideo.mat','Frames_Rec');      % 2. Reconstracted Frames


%% H.264 Decoding
clc
clear all;
close all;
%% Decoding Bitsteam
load 02BitStream.mat
% load 02ReconstructedVideo.mat
global h w QP block_size
idx = 1;
block_size = 16;
%---------------------------------------------------------
%% Decode header
[h,w,QP,Frame_start,Frame_end,m] = dec_header(bitstream);
idx = idx + m - 1;
N = 1 + (Frame_end - Frame_start);
%% Decode I-Frame
if (bitstream(idx:idx+3)=='1111')
   disp(['Decoding I Frame: ',num2str(Frame_start)])
   idx = idx + 4;
   [Ceq_r(:,:,1),idx]=decode_i_frame(idx,bitstream);
   Frames_Dec(:,:,1)=Ceq_r(:,:,1);
end
%% Decode the following P-Frames
for k = 2:N       
    if (bitstream(idx:idx+3)=='0000')
        disp(['Decoding P Frame: ', num2str(k)])
        idx = idx + 4;
        [Ceq_r(:,:,k),idx]= decode_p_frame(idx,bitstream,Ceq_r(:,:,k-1));
        Frames_Dec(:,:,k)=Ceq_r(:,:,k);
    end  
end
% End the decoding
%% Option: Save the decoded video and paly it
save('03DecodedVideo.mat','Frames_Dec');
A01SaveAndPlay('HaveDecoded.avi',Frames_Dec);