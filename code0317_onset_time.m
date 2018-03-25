% code for finding miscible fingering elongate time
% 0317
% by Rui Luo

clear
clc
close all

particleSize = 125;  % unit: um
gapThickness = 1.397; % unit: mm
FlowRate = 150;      % unit: ml/min

MainDirectory = 'C:\Users\lr546\Desktop\';

% phiInitial = [0.14,0.17,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35];
phiInitial = [0.31];
timeFrame = linspace(90,600,511);
% timeFrame = [540];
[~,phiIndexTotal] = size(phiInitial);
for phiIndex = 1:phiIndexTotal
      DataDirectory = [MainDirectory,num2str(particleSize),'particle ',...
      num2str(gapThickness),'gap\phi',num2str(phiInitial(phiIndex)*100)];
      concentrationDirectory = [DataDirectory,'\concentration data'];
      indexTimeTotal = length(timeFrame);
      
      for indexTime = 1:indexTimeTotal
            fileDirectory = [concentrationDirectory,'\',num2str(timeFrame(indexTime)),'.csv'];
            concentrationData = csvread(fileDirectory,1,0);
            rho = concentrationData(:,1);
            concentrationVar = concentrationData(:,4);
            dataLength = length(rho);
            rho(dataLength-4:dataLength)=[];
            rho(1:3)=[];
            concentrationVar(dataLength-4:dataLength)=[];
            concentrationVar(1:3)=[];
            dataLength = length(rho);
            threshold = mean(concentrationVar(ceil(dataLength*0.2):ceil(dataLength*0.65)));
            a = find(concentrationVar>=threshold*1.3);
            a(a<ceil(dataLength*0.8)) = [];
            rhoFinger = rho(a);
            if isempty(rhoFinger)== 0 
            fingerLength(indexTime) = max(rhoFinger(:))-min(rhoFinger(:));
            else
            fingerLength(indexTime) =0;     
            end
            
%             figure;plot(rho,concentrationVar,'*-')
%             figure;plot(rho(a),concentrationVar(a),'^-')

      end
       figure; plot(timeFrame/30,smooth(fingerLength),'*-')
      if isempty(fingerLength) == 0
      fingerOnsetTime(phiIndex) = timeFrame(find(fingerLength>0,1));
      end
end
figure; plot(phiInitial,fingerOnsetTime,'*-')