% this function is used for calculating mobility ratio.
% By Rui Luo 2017/12/29

function mobilityRatioRingAverage = getMobilityRatio(valueRingAverage)
    viscosity = valueRingAverage(:,2);
    ringPosition = valueRingAverage(:,1);
    [ringNum,columnValue] = size(valueRingAverage);
    mobilityRatioRingAverage=zeros(ringNum,columnValue);
    for iNum =1:ringNum-1
        mobilityRatioRingAverage(iNum,1) = ringPosition(iNum);
        mobilityRatioRingAverage(iNum,2) = viscosity(iNum+1)/viscosity(iNum);
    end
%     figure
%     plot(mobilityRatioRingAverage(:,1),mobilityRatioRingAverage(:,2));
end