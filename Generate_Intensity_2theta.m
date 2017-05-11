%% Get the XRD Pattern
% Run through StructureFactor with the same parameters except the
% Lattice.Normal to obtain multiple 2theta angles.
% diffraction geometry of the (102) peak (noncoplanar)
% It needs the Lattice and Probe to get information for StructureFactor,
% Threshold to limit the peaks to only the important ones (it is
% intensity), Resolution to define how close two peaks can be together before
% they are interpreted as a single one from the same family of planes and
% span/IndexMax for what planes should be analyzed 
% 
%
%
% Last updated 5-8-2017 Cosmin Popescu

function Table=Generate_Intensity_2theta(Lattice, Probe,Threshold, Resolution, hkl,FigNum)
%% Get the XRD Pattern
% Run through StructureFactor with the same parameters except the
% Lattice.Normal to obtain multiple 2theta angles.
% diffraction geometry of the (102) peak (noncoplanar)
% Last updated 4-10-2017 Cosmin Popescu

addpath(genpath('C:\Users\Cosmin\Desktop\Grand Diffraction Master'))
addpath(genpath('C:\Users\Cosmin\Desktop\Cr2AlC'))
addpath(genpath('C:\Users\Cosmin\Desktop\Diffraction-master\StructureLibrary'))
addpath(genpath('C:\Users\Cosmin\Desktop\Diffraction-master\TestScripts'))
%addpath(genpath('C:\Users\Cosmin\Desktop\Diffraction-master'))

% % Load your material
% load Cr2AlC.mat
% 
% % DEFINE YOUR X-RAYS
% Probe.Type = 'x-ray';
% Probe.Energy = 8047; % [eV] % define either Energy or lambda
% Probe.Polarization = 's'; % s (perpendicular) or p (parallel)
% 
% Threshold=1;
% Resolution=0.1;
% IndexMax=9;
% FigNum=[];

%% Loop over different hkl values individually. 
% There will be repetitions
% like 001 and 002
% The software goes through all combinations given by user of Miller
% indeces and returns the resulting information from Structure Factor. This
% information is saved into DataOfIntensity which is used as file of raw
% data. It contains among other things: complex values for the intensity
% and the Bragg angle. This are removed further in the script.
% DataOfIntensity is the first storage point/variable but the main
% operations are done with PlotIntensity which is a matrix type variable as
% well.
TypeOfFile=input('What kind of file extension do you want?\nOptions: nothing (press Enter), txt, xlsx, xls, dat, csv \n','s');

DataOfIntensity=zeros(10,4);
countofdata=1;


for h=-hkl:hkl
    for k=-hkl:hkl
        for l=-hkl:hkl
            if h==0 && k==0 && l==0
                % This is to eliminate the 000
            else
            Lattice.Reflection=[h k l];
            DataOfIntensity(countofdata,1)=h;
            DataOfIntensity(countofdata,2)=k;
            DataOfIntensity(countofdata,3)=l;
            [Result,Lattice,Probe]= StructureFactor(Lattice,Probe);
            DataOfIntensity(countofdata,4)=Result.Intensity;
            DataOfIntensity(countofdata,5)=2*Result.BraggAngle;
            countofdata=countofdata+1;
            end
        end
    end
end

%% Eliminate the values with intensity below a threshold.
% PlotIntensity is created by going through all raw data and removing all
% the BraggAn
countofplot=1;
PlotIntensity=zeros(1,5);
for i=1:countofdata-1
    if(DataOfIntensity(i,4)>Threshold && isreal(DataOfIntensity(i,5)))
        PlotIntensity(countofplot,:)=DataOfIntensity(i,:);
        countofplot=countofplot+1;
    end
end
%% Sort everything by the Bragg angle.
sorted=false;
while sorted==false
    n=1;
   while n<countofplot-1
        if PlotIntensity(n,5)>PlotIntensity(n+1,5)
            
            store=PlotIntensity(n+1,:);
            PlotIntensity(n+1,:)=PlotIntensity(n,:);
            PlotIntensity(n,:)=store;
            n=1;
        else 
            n=n+1;
        end
   
   end
   sorted=true;
    
end
    

%% Make Plot

if nargin>5
    figure(FigNum)
    
    plot(PlotIntensity(:,5),PlotIntensity(:,4),'ob','linewidth',2);
    xlabel('2\theta');
    ylabel('Intensity');
    Title=strcat(Lattice.Symbol,'-',Lattice.Type);
    title(Title);
    axis([0 180 0 1.1*max(PlotIntensity(:,4))]);
    hold on;
end

%% Add the Miller indeces to the plot
xseparation=-2;
yseparation=max(PlotIntensity(:,4))*0.04;
% This is a way to separate multiple combinations of hkl that give the same
% theta. 
PlotIntensity(:,6)=ones(countofplot-1,1);
for i=2:countofplot-1
    if PlotIntensity(i,5)==PlotIntensity(i-1,5)
        PlotIntensity(i,6)=PlotIntensity(i-1,6)+1;
    end
end
% Write text at a separation from the point. If multiple point in the same
% location, stack the hkl in order defined by previous for.
if nargin >5
    for i=1:countofplot-1
        h=PlotIntensity(i,1);
        k=PlotIntensity(i,2);
        l=PlotIntensity(i,3);
        text( PlotIntensity(i,5)+xseparation,PlotIntensity(i,4)+yseparation*PlotIntensity(i,6),strcat([num2str(h), num2str(k), num2str(l)]))
    end
end
%% Intensities on the same angle are added up

PlotIntensity(:,7)=zeros(countofplot-1,1);
PlotIntensity(:,8)=PlotIntensity(:,4);
[length,~]=size(PlotIntensity);
for i=1:length
    n=i;
    currentposition=i;
    if(n<length-1)
        while(PlotIntensity(n+1,6)>1)
            PlotIntensity(currentposition,8)=PlotIntensity(currentposition,8)+PlotIntensity(n+1,8);
            PlotIntensity(n+1,8)=0;
            if(n+1<length)
            n=n+1;
            else
                break
            end
        end
    else
        if(PlotIntensity(n,6)==1)
            PlotIntensity(n,4)=PlotIntensity(n,4);
        end
    end
    
end

%% Generating XRD plot with resolution based on range.
if nargin>5
    figure(FigNum+1)
else
    figure
end
countAngle=1;
XRD_plot=zeros(ceil(180/Resolution+length),5);
indexPlot=1;
index=1;

for angle=0:Resolution:180
    if (angle<PlotIntensity(countAngle,5))&& (angle+Resolution>PlotIntensity(countAngle,5)&&countAngle<length-1)
        Table.h=PlotIntensity(countAngle,1); %h
        Table.k=PlotIntensity(countAngle,2); %k
        Table.l=PlotIntensity(countAngle,3); %l
        Table.TwoTheta=PlotIntensity(countAngle,5);
        Table.Intensity=PlotIntensity(countAngle,8); %Piggyback on data processing.
        Array(index)=Table;
        index=index+1;
        XRD_plot(indexPlot,1)=PlotIntensity(countAngle,5); %2theta
        XRD_plot(indexPlot,2)=PlotIntensity(countAngle,8); %Intensity of the actual point
        indexPlot=indexPlot+1;
        if countAngle<length-1
            while(PlotIntensity(countAngle+1,6)>1&&(countAngle<(length-2)))
                countAngle=countAngle+1;
            end
        end
        countAngle=countAngle+1;
    else
        XRD_plot(indexPlot,1)=angle; %Looping through the angles.
        XRD_plot(indexPlot,2)=0; %No intensity
        indexPlot=indexPlot+1;
        
    end       
end
% subplot(2,1,1);
hold on;
plot(XRD_plot(:,1),XRD_plot(:,2)*100/max(XRD_plot(:,2)),'-b','linewidth',1.5);
xlabel('2\theta');
ylabel('Intensity percent I/I_M_a_x');
Title=strcat(Lattice.Symbol,' - ',Lattice.Type,' 8048.3 eV');
title(Title);
axis([0 180 0 105]);
hold on;
maximum=max(XRD_plot(:,2));
PlotIntensity(:,9)=100*PlotIntensity(:,8)/maximum;
legend('MATLAB');
%% A new plot similar to actual XRD is generated.
% Xtheta=zeros(5,180+length);
% count=1;
% % Make 0 in intensity everything that is not in the data.
% for i=1:180+length
%     if count<=length
%         if((PlotIntensity(count,5)-i+count)<Resolution)
%             Xtheta(1,i)=PlotIntensity(count,5);
%             Xtheta(2,i)=PlotIntensity(count,8);
%             Xtheta(3,i)=PlotIntensity(count,1);
%             Xtheta(4,i)=PlotIntensity(count,2);
%             Xtheta(5,i)=PlotIntensity(count,3);
%             count=count+1;
%         else
%             Xtheta(1,i)=i-count;
%         end
%     else
%         Xtheta(1,i)=i-count;
%         
%     end
% end
% Xtheta=Xtheta';
% figure;
% plot(Xtheta(:,1),Xtheta(:,2),'-b');
% xlabel('2\theta');
% ylabel('Intensity');
% title(Lattice.Symbol);
% axis([0 180 0 1.1*max(Xtheta(:,2))]);
% hold on;
%% Make name for the text file
% Dealt with this below. See comment after section below.
% Mostly old code here.
%     Filename=strcat(Lattice.Symbol,'_',Lattice.Type,'.txt');
%     MAX=max(PlotIntensity(:,8));
%     PlotIntensity(:,8)=PlotIntensity(:,8)/MAX;

%% Make table with Miller indeces.
Table=struct2table(Array);
time=datestr(datetime('now'));
time(15)='_';
time(18)='_';
time(12)='_';
NameOFile=strcat(Lattice.Symbol,'_',Lattice.Type,'_',num2str(Probe.Energy),'eV','_',time,'.');

if size(TypeOfFile,2)==4
    if TypeOfFile=='xlsx'
        writetable(Table,strcat(NameOFile,'xlsx'));
    end
elseif size(TypeOfFile,2)==3
    if TypeOfFile== 'txt'
        writetable(Table,strcat(NameOFile,'txt'));
    elseif TypeOfFile=='dat'
        writetable(Table,strcat(NameOFile,'dat'));
    elseif TypeOfFile == 'csv'
        writetable(Table,strcat(NameOFile,'csv'));
    elseif TypeOfFile == 'xls'
        writetable(Table,strcat(NameOFile,'xls'));
    end
else %default do nothing
end




% Former text file generation on the basis of the data generated. Updated
% version of data output above in table generation with different
% extension.
% fid=fopen(Filename,'a');
% fprintf(fid,'\r\n');
% fprintf(fid,'Start of a new data set \r\n');
% fprintf(fid,strcat('The energy is  ',num2str(Probe.Energy)));
% fprintf(fid,'\r\n %s %s %s %s %s\r\n','h','k','l','2theta','Intensity');
% for i=1:length
%     if(PlotIntensity(i,6)==1)
%         fprintf (fid,'%d %d %d %f %f\r\n',PlotIntensity(i,1),PlotIntensity(i,2),PlotIntensity(i,3),PlotIntensity(i,5),PlotIntensity(i,8));
%         
% %     else
% %         fprintf(fid,'%d %d %d\r\n',PlotIntensity(i,1),PlotIntensity(i,2),PlotIntensity(i,3));
%         
%     end
% end
% fclose(fid);

