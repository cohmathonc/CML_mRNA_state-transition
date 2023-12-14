function cIndexValue = calculateCIndex(predicted1, predicted2)
    n = length(predicted1);
    concordant = 0;
    discordant = 0;

    for i = 1:n
        for j = i + 1:n
            if (predicted1(i) < predicted1(j) && predicted2(i) < predicted2(j)) || ...
               (predicted1(i) > predicted1(j) && predicted2(i) > predicted2(j))
                concordant = concordant + 1;
            elseif (predicted1(i) < predicted1(j) && predicted2(i) > predicted2(j)) || ...
                   (predicted1(i) > predicted1(j) && predicted2(i) < predicted2(j))
                discordant = discordant + 1;
            end
        end
    end

    cIndexValue = concordant / (concordant + discordant);
end
