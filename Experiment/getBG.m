%% https://github.com/harshjn/OpticalTweezers/

addOpen='Y:\smita\optical tweezer expts\2um-circle\';    %Path of video file
addSave='Y:\smita\optical tweezer expts\Analysis\';

TIMES=[0,100];  %Starting and Ending time intervals

%%

filename='4'
v=VideoReader(strcat(addOpen,filename,'.avi'));
Fs=v.FrameRate; 


%%
v.CurrentTime=100
a1 = readFrame(v);
imshow(a1(461:761,589:892))

v.CurrentTime=120
a2 = readFrame(v);
imshow(a2(461:761,589:892))

bg1 = a1(612:761,589:892)
bg2 = a2(461:611,589:892)
imshow(bg2)
bg=[bg2;bg1]

imshow(bg)
