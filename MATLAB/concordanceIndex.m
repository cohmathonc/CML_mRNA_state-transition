function cIndex = concordanceIndex(times, censored, predicted_score)
%% example
% times = 1:5000;
% censored = floor(rand(1, 5000)+ 0.5);
% for i = 1:5000
% predicted_score(i) = rand();%i + i*(rand-0.5);
% end
if length(times) ~= length(censored) || length(times) ~= length(predicted_score) || length(censored) ~= length(predicted_score)
   disp('DIMENSION ERROR') 
end
%% code
%throw out NaNs
nanIndex = isnan(predicted_score);
times(nanIndex) = [];
censored(nanIndex) = [];
predicted_score(nanIndex) = [];
%get pairs
all_pairs = nchoosek(1:length(times), 2);
%throw out useless pairs
throwOut =  ( censored(all_pairs(:,1)) & censored(all_pairs(:,2)) ) |...
( (times(all_pairs(:,1)) < times(all_pairs(:,2))) & censored(all_pairs(:,1)) ) |...%patient 1 is censored but has smaller follow up than event
( (times(all_pairs(:,2)) < times(all_pairs(:,1))) & censored(all_pairs(:,2)) ) ;%patient 2 is censored but has smaller follow up than event
all_pairs(throwOut, :) = [];
%count tied pairs before deleting
tied_pair_ind = predicted_score(all_pairs(:,1)) == predicted_score(all_pairs(:,2));
tied_pairs = sum(tied_pair_ind) ;
all_pairs(tied_pair_ind, :) = [];
number_pairs = length(all_pairs) + tied_pairs;
time_concordance = times(all_pairs(:,1)) > times(all_pairs(:,2));
score_concordance = predicted_score(all_pairs(:,1)) > predicted_score(all_pairs(:,2)) ;
concordant_pairs = sum( time_concordance == score_concordance);
cIndex = (concordant_pairs + tied_pairs/2)/number_pairs;
end