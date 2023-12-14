% Generate example survival data
rng(123); % Set seed for reproducibility
n = 100;
time = exprnd(1/0.02, n, 1);
status = randi([0, 1], n, 1);
covariate1 = randn(n, 1);
covariate2 = randn(n, 1);

data = table(time, status, covariate1, covariate2, 'VariableNames', {'Time', 'Status', 'Covariate1', 'Covariate2'});

% Split the data into training and testing sets (adjust the split ratio)
trainFraction = 0.8;
trainIdx = randperm(n, round(trainFraction * n));
trainData = data(trainIdx, :);
testData = data(setdiff(1:n, trainIdx), :);

% Fit Cox Proportional Hazards Model on training data
formula = 'Time ~ Covariate1 + Covariate2';
% Fit Cox Proportional Hazards Model on training data
covariates = table2array(trainData(:, {'Covariate1', 'Covariate2'}));
coxModel = coxphfit(covariates,  'Censoring', ~trainData.Status);

% Display the summary of the Cox model
disp(coxModel)
coxModel = coxphfit(trainData, 'Survival', 'on', 'Censoring', ~trainData.Status, 'Covariates', formula);

% Display the summary of the Cox model
disp(coxModel)

% Predict survival probabilities on testing data
predictedSurvival = coxModel.predict(testData);

% Combine observed and predicted survival data
observedAndPredicted = table(testData.Time, testData.Status, predictedSurvival, ...
    'VariableNames', {'Time', 'Status', 'PredictedSurvival'});

% Display the first few rows of observed and predicted survival data
disp(head(observedAndPredicted))

% Plot observed and predicted survival curves for a subset of individuals
figure;
plotSurvival(observedAndPredicted, 'GroupVariable', 'Status', 'TimeVariable', 'Time');
legend('Observed: Status = 0', 'Observed: Status = 1', 'Predicted: Status = 0', 'Predicted: Status = 1');
title('Observed and Predicted Survival Curves');
