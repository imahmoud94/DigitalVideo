video = VideoReader('3Docean.mov');

wVideo = VideoWriter('img\\16x16\\16x16myVideo.avi');
wVideo.FrameRate =1;
wcVideo = VideoWriter('compressed\\16x16\\16x16CompressedMyVideo.avi');
wcVideo.FrameRate =1;

open(wVideo);
open(wcVideo);

v = read(video,[1 10]);

for i = 1:10
    frame = v(:,:,:,i);
    
    imwrite(frame,strcat('.\\img\\img',int2str(i),'.png'));
    writeVideo(wVideo,frame);
    cframe = compress(frame, 16);
    imwrite(cframe,strcat('.\\compressed\\16x16\\16x16compressedImg',int2str(i),'.png'));
    writeVideo(wcVideo,cframe);
end

close(wVideo);
close(wcVideo);