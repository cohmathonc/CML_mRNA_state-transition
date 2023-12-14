% Concordance Index calculation function
function cIndex = calculateConcordance(observedTime, predictedTime, censoringStatus)
    [~, sortIdx] = sort(predictedTime, 'descend');
    observedTime = observedTime(sortIdx);
    censoringStatus = censoringStatus(sortIdx);

    concordantPairs = 0;
    discordantPairs = 0;

    for i = 1:length(observedTime)-1
        for j = i+1:length(observedTime)
            if censoringStatus(i) == 0 && censoringStatus(j) == 0
                if observedTime(i) < observedTime(j)
                    concordantPairs = concordantPairs + 1;
                elseif observedTime(i) > observedTime(j)
                    discordantPairs = discordantPairs + 1;
                end
            end
        end
    end

    totalPairs = concordantPairs + discordantPairs;
    cIndex = concordantPairs / totalPairs;
end