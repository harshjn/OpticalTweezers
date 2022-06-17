%% https://github.com/harshjn/OpticalTweezers/


clear all
close all
addOpen='Y:\smita\optical tweezer expts\2um-circle\';    %Path of video file
addSave='Y:\smita\optical tweezer expts\Analysis\';
filenames=['4'];   % 
load('bg.mat')  % Background generated via getBG.m
%%
filename='4';
v=VideoReader(strcat(addOpen,filename,'.avi'));
Fs=v.FrameRate;  

%%
startTime=0;
endTime=222;
len=round((endTime-startTime)*Fs)+1;
cents=zeros(len,2);
rads=zeros(len,1);
time=zeros(len,1);
v.CurrentTime=startTime;
%%
k=2;
l=zeros(length(time),1);
v.CurrentTime = 1

%%

% while(v.hasFrame)
while(v.CurrentTime<endTime)
    A_=readFrame(v);
    A=A_(461:761,589:892);
%     A=A-bg./2;
%     A=rgb2gray(A);
    A2=A-bg;
%     A=rgb2gray(A);
    A3=imadjust(A2);
%     A3=~A3;
%     imshow(A3)

%
%     A=im2bw(A,0);
    
    
    try
        [centers, radii] = imfindcircles(A3,[13 19],'ObjectPolarity','bright','Sensitivity',0.95)
        cents(k,1)=centers(1);
        cents(k,2)=centers(2);
        if mod(k,100)==0 
            sprintf('%d',k)
        end
        rads(k)=radii;
%         
%         [centsx,centsy,radii]=radialcenter(im2double(A));
%         rads(k)=radii;
%         sprintf('%d',k)
%         cents(k,1)=centsx;
%         cents(k,2)=centsy;
                
    catch causeexception
        l(k)=k;
        sprintf('%d',k)
        cents(k,1:2)=cents(k-1,1:2); %Hiding the errors for continuity
        rads(k)=rads(k-1);
    end
    time(k)=v.CurrentTime;
    k=k+1;
    
end
l2=find(l);


%%


for i = 1:3000
viscircles(cents(i,:),rads(i))
pause(0.01)
end
