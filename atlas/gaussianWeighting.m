function [weights] = gaussianWeighting(ages, bdLow, bdHigh, sig, center)
% generate Gaussain weights based on ages

weights = zeros(size(ages));
for iI = 1:length(ages)
    weights(iI) = exp( -(ages(iI) - center)^2 / (2*sig^2) );
    if ages(iI) < bdLow + sig
        weights(iI) = weights(iI) + exp( -(2*bdLow - ages(iI) - center)^2 / (2*sig^2) );
    elseif ages(iI) > bdHigh - sig
        weights(iI) = weights(iI) + exp( -(2*bdHigh - ages(iI) - center)^2 / (2*sig^2) );
    end
end
