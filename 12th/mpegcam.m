function mpegcam()
clear all;
close all;
N = 8;% block size is N x N pixels
N2 = 2*N;
W = 16; % search window size isWxW pixels
quantizer_scale = 4; % used only for the base layer

[chosenfile, chosenpath] = uigetfile('*.avi', 'Select a video');
  if ~ischar(chosenfile)
    return;   %user canceled dialog
  end
  filename = fullfile(chosenpath, chosenfile);

frames = VideoReader(filename); 
vidHeight=frames.Height;
vidWidth=frames.Width;

F=read(frames,[1 10]);
M = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),'colormap',[]);
k = 1;
figure;
subplot(2,3,1)
%hold on, title('original video');
for k=1:10
M(k).cdata = F(:,:,:,k);
image(M(k).cdata);
title('RGB Video Before comp.');
pause(1/frames.FrameRate);
end


[X,Y,Z] = size(M(1).cdata);
if mod(X,8)~=0
    Height = floor(X/8)*8;
else
    Height = X;
end
if mod(Y,8)~=0
   Width = floor(Y/8)*8;
else
   Width = Y;
end
Depth = Z;
clear X Y Z
%
if Depth == 3
   A = rgb2ycbcr(M(1).cdata);% Convert RGB to YCbCr & retain only Y
   y_ref = A(:,:,1);
else
   A = M(1).cdata;
   y_ref = A;
end
% pad the reference frame left & right and top & bottom
y_ref = double(padarray(y_ref,[W/2 W/2],'replicate'));
% arrays to store SNR and PSNR values
Base_snr = zeros(1,k-2); Enhanced_snr = zeros(1,k-2);
Base_psnr = zeros(1,k-2); Enhanced_psnr = zeros(1,k-2);
% Encode the monochrome video using MPC
subplot(2,3,3)
for f = 2:k-1
   if Depth == 3
      B = rgb2ycbcr(M(f).cdata);
      y_current = B(:,:,1);
   else
      y_current = M(f).cdata;
   end
   y_current = double(padarray(y_current,[W/2 W/2],'replicate'));
   for r = N:N:Height
      rblk = floor(r/N);
      for c = N:N:Width
         cblk = floor(c/N);
         D = 1.0e+10;% initial city block distance
         for u = -N:N
            for v = -N:N
             d=y_current(r+1:r+N,c+1:c+N)-y_ref(r+u+1:r+u+N,c+v+1:c+v+N);
             d = sum(abs(d(:)));% city block distancebetween pixels
             if d < D
                D = d;
                x1 = v; y1 = u; % motion vector
             end
            end
          end
% MC compensated difference coding
        temp = y_current(r+1:r+N,c+1:c+N)-y_ref(r+1+y1:r+y1+N,c+1+x1:c+x1+N);
        TemP = dct2(temp); % DCT of difference
        s = sign(TemP); % extract the coefficient sign
        TemP1 = s .* round(abs(TemP)/(16*quantizer_scale))*(16*quantizer_scale); % quantize/dequantize DCT
        temp = idct2(TemP1); % IDCT
        Base(r-N+1:r,c-N+1:c) = y_ref(r+1+y1:r+y1+N,c+1+x1:c+x1+N)+temp; % reconstructed block - base quality
        delta_DCT = TemP - TemP1; % incremental DCT
        s1 = sign(delta_DCT); % extract the sign of incremental DCT
        delta_DCT = s1 .* round(abs(delta_DCT)/4)*4;
        temp1 = idct2(TemP1 + delta_DCT);
        Enhanced(r-N+1:r,c-N+1:c) = y_ref(r+1+y1:r+y1+N,c+1+x1:c+x1+N) +temp1;
      end
  end
% Calculate the respective SNRs and PSNRs
    Base_snr(f-1) = 20*log10(std2(y_current(N+1:Height+N,N+1:Width+N))/std2(y_current(N+1:Height+N,N+1:Width+N)-Base));
    Enhanced_snr(f-1) = 20*log10(std2(y_current(N+1:Height+N,N+1:Width+N))/std2(y_current(N+1:Height+N,N+1:Width+N)-Enhanced));
    Base_psnr(f-1) = 20*log10(255/std2(y_current(N+1:Height+N,N+1:Width+N)-Base));
    Enhanced_psnr(f-1) = 20*log10(255/std2(y_current(N+1:Height+N,N+1:Width+N)-Enhanced));
% replace previous frames by the currently reconstructed frames
    y_ref = Base;
    y_ref = double(padarray(y_ref,[W/2 W/2],'replicate'));
    y_ref2 = Enhanced;
    y_ref2 = double(padarray(y_ref,[W/2 W/2],'replicate'));
image(y_ref);
title('Monochrome video After Decomp Base');
subplot(2,3,6)
image(y_ref2);
title('Monochrome video After Decomp Enhanced');
pause(1/frames.FrameRate);
end
FNO = int16(1:k);
subplot(2,3,4),plot(FNO(2:end-1),Base_snr,'k*','LineWidth',1), hold on
%figure,plot(FNO(2:end-1),Base_snr,'k*','LineWidth',1), hold on
plot(FNO(2:end-1),Enhanced_snr,'kd','LineWidth',2), title('SNR (dB)')
% axis([M(2) M(end) min(Base_snr)-2 max(Enhanced_snr)+2]) % for Rhinos sequence
legend('Base Quality','Enhanced Quality','Best')
xlabel('Frame #'), ylabel('SNR (dB)'), hold off

subplot(2,3,5),plot(FNO(2:end-1),Base_psnr,'k*','LineWidth',1), hold on
% figure,plot(FNO(2:end-1),Base_psnr,'k*','LineWidth',1), hold on
plot(FNO(2:end-1),Enhanced_psnr,'kd','LineWidth',2), title('PSNR (dB)')
% axis([F(2) F(end) min(Base_psnr)-2 max(Enhanced_psnr)+2]) % for Rhinos sequence
legend('Base Quality','Enhanced Quality','Best')
xlabel('Frame #'), ylabel('PSNR (dB)'), hold off
