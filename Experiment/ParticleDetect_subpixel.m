%% https://github.com/harshjn/OpticalTweezers/

addOpen='F:\Data\';    %Path of video file
addSave='F:\Analysis\';

filenames=['0.6particle2'];   % 
TIMES=[0,528];  %Starting and Ending time intervals

%%
load('bg.mat')  % Background generated via getBG.m
bg=bg12;
A1=readFrame(v);
bg1 = A1(:,157:312,:);
%%
savefile=1;


for i=1:1:size(filenames,1)
filename=filenames(i,:);
v=VideoReader(strcat(addOpen,filename,'.avi'));
Fs=v.FrameRate; 

%%
j=i;  % Select your interested time interval from TIMES
startTime=TIMES(j,1);
endTime=TIMES(j,2);

len=round((endTime-startTime)*Fs)+1;
cents=zeros(len,2);
rads=zeros(len,1);
time=zeros(len,1);
v.CurrentTime=startTime;
%%
k=1;
l=zeros(length(time),1);
%%

% while(v.hasFrame)
while(v.CurrentTime<endTime)
    
    A=readFrame(v);
%     A=A-bg./2;
%     A=rgb2gray(A);
    A=A-bg;
    A=rgb2gray(A);

    A=imadjust(A);
    A=~A;

%
%     A=im2bw(A,0);
    
    
    try
        [centers, radii] = imfindcircles(A,[13 19],'ObjectPolarity','bright','Sensitivity',0.93);
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

%% Removing the errors by iterating over errors with varying sensitivities or other parameters

% m=[];
% for p=1:length(l2)
%     k=l2(p);
%     v.CurrentTime=time(k-1);
%     A=readFrame(v);
%     A=imadjust(A);
% 
% %     A=rgb2gray(A);
% %     A=A-bg;
% %     A=A(337:737,478:887); %Location of the circular trap
% 
% %     A=A-bg;
% 
% %     A=im2bw(A,0);
%     
%     try
% %         [centers, radii] = imfindcircles(A,[14 25],'ObjectPolarity','dark','Sensitivity',0.92);
%         cents(k,1)=centers(1);
%         cents(k,2)=centers(2);
%         rads(k)=radii;
% %         sprintf('%d',k)
%     catch causeexception
%         m(end+1)=k;
%         sprintf('%d',k)
%         cents(k,1:2)=cents(k-1,1:2);
%         rads(k)=rads(k-1);
%     end
% end
% l2=m;
%%
cents=cents(1:end-1,:);

%% Angular plot
%Finding Centers from scatter Plot
% hold on
% scatter(CenterX,CenterY)
% hold off


CenterX=163.5;
CenterY=130.82;
[Angles,rhos]=cart2pol((cents(:,1)-CenterX),(cents(:,2)-CenterY));
% thetas=Angles(rhos<200);
thetas=Angles;

    
%% histogram of angular position
fig1=figure();
HistBW=0.01;
ThetaHist=histogram(thetas,'binWidth',HistBW);

ProbValues=ThetaHist.Values./sum(ThetaHist.Values);
bar(ThetaHist.BinEdges(1:end-1),ProbValues);

titl=strcat('histogram of angle distribution',filename);
title(titl)
xlabel('angle')
ylabel('Probability')
set(gca, 'YScale', 'log')

if savefile==1
    saveas(fig1,strcat(addSave,titl,'.fig'))

    close(fig1)
end

%% plot of the unwrapped angle
fig2=figure();

theta=unwrap(thetas);
plot(theta)
titl=strcat('AngleEvolution',filename);
title(titl);
xlabel('angle')
ylabel('Probability')
if savefile==1
    saveas(fig2,strcat(addSave,titl,'.png'))

    close(fig2)
end
%% Plot LDF
figure();
HistBW=0.01;

ThetaHist=histogram(thetas,'binWidth',HistBW);
ProbX=ThetaHist.BinEdges(1:end-1);

ProbY=-log(ThetaHist.Values./sum(ThetaHist.Values));
bar(ProbX,ProbY)


%% Where's a particular angle?
scatter(cents(:,1),cents(:,2),'.')
hold on


%%
a=find((theta<-0*pi) & (theta>-50*pi));
startT=min(a)
endT=max(a)
figure();
scatter(cents(startT:endT,1),cents(startT:endT,2),'.')

theta2=theta(startT:endT);
fig2=figure();
HistBW=0.01;
ThetaHist=histogram(theta2,'binWidth',HistBW);

titl=strcat('histogram of angle distribution',filename);
title(titl)
xlabel('angle')
ylabel('Probability')
set(gca, 'YScale', 'log')

if savefile==1
    saveas(fig1,strcat(addSave,titl,'.png'))

    close(fig1)
end



%% Save the relevant angle
if savefile==1    
    savefilename=strcat(addSave,'data',filename,'.mat');
    save(savefilename,'theta','cents','rads','CenterX','CenterY','thetas','filename','l2','startTime','endTime','Fs')
end

end
