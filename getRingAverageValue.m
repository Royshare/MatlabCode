% this function is used for get ring average value.
% By Rui Luo 2017/12/27

function ringAverageData = getRingAverageValue(dataPolar,ringWidth)
   rPosition = dataPolar(:,1);
   thetaPosition = dataPolar(:,2);
   physicalValue = dataPolar(:,3);
   ringStart = min(dataPolar(:,1));
   ringEnd = max(dataPolar(:,1));
   ringRangePixel = ringStart:ringWidth:ringEnd;
   for ringNum = 1:length(ringRangePixel)-1
        indexLogical = logical(rPosition > ringRangePixel(ringNum) &...
                                rPosition <= ringRangePixel(ringNum+1));
        indexInRing = find(indexLogical);
        valueInRing = [];
        valueInRing = physicalValue(indexInRing);
        valueRingAverage(ringNum) = mean(valueInRing);
        valueRingVar(ringNum) = var(valueInRing);
        ringLocation(ringNum) = (ringRangePixel(ringNum)+ringRangePixel(ringNum+1))/2;
   end
   ringAverageData = [ringLocation',valueRingAverage',valueRingVar'];
   
%    figure
%    plot(ringLocation',valueRingAverage','*-')
   
end 