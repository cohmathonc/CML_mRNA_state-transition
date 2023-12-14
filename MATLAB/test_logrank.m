% Sample survival data (time to event and group indicator)
time_to_event_group1 = [10, 12, 14, 16, 18];
time_to_event_group2 = [8, 11, 15, 17, 20];
event_status_group1 = ones(size(time_to_event_group1)); % 1 indicates an event
event_status_group2 = ones(size(time_to_event_group2));

% Combine data and group indicators
time_to_event = [time_to_event_group1, time_to_event_group2];
event_status = [event_status_group1, event_status_group2];
group = [ones(size(time_to_event_group1)), 2 * ones(size(time_to_event_group2))];

% Sort data by time-to-event
[time_to_event_sorted, sort_indices] = sort(time_to_event);
event_status_sorted = event_status(sort_indices);
group_sorted = group(sort_indices);

% Compute Kaplan-Meier survival curves
[km_curve_group1, km_curve_group2] = calculate_kaplan_meier(time_to_event_sorted, event_status_sorted, group_sorted);

% Compute Log-Rank Test
test_statistic = calculate_logrank_test(km_curve_group1, km_curve_group2);

% Calculate p-value
p_value = chi2cdf(test_statistic, 1, 'upper');

% Display the results
fprintf('Log-Rank Test Results:\n');
fprintf('Test Statistic: %f\n', test_statistic);
fprintf('P-value: %f\n', p_value);

% Interpret the results
if p_value < 0.05
    fprintf('There is a significant difference between the survival curves.\n');
else
    fprintf('There is no significant difference between the survival curves.\n');
end
