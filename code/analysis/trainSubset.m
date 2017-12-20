function err = trainModel(dataPath, outPath, featureRows, labelRow)
  % Expected to load data matrix that is 0-1 normalized
  load(dataPath);

  rng('shuffle')

  subset = 0.10

  numFolds = 10;

  featureRows = str2num(featureRows);
  labelRow = str2num(labelRow);

  numTuples = size(data, 2);
  shuffle = randperm(numTuples);

  data = data(:, shuffle);

  % perform subset
  numTuples = floor(size(data, 2) * subset);
  data = data(:, 1:numTuples);
  trainingSetSize = floor(numTuples*((numFolds-1)/numFolds));

  trainingSet = data(:, 1:trainingSetSize);
  validationSet = data(:, trainingSetSize+1:numTuples);

  disp('Training Set Size')
  size(trainingSet)
  disp('Validation Set Size')
  size(validationSet)

  trainingFeatures = trainingSet(featureRows,:)';
  trainingLabels = trainingSet(labelRow, :)';


  model = svmtrain(trainingLabels,...
                   trainingFeatures,...
                   '-s 3 -t 2 -h 0');
  disp('TRAINING SET:')
  svmpredict(trainingLabels,...
             trainingFeatures,...
             model);

  disp('trainingFeatures')
  size(trainingFeatures)
  disp('trainingLabels')
  size(trainingLabels)

  validationFeatures = validationSet(featureRows, :)';
  validationLables = validationSet(labelRow, :)';
  disp('Validation SET:')
  svmpredict(validationLables,...
             validationFeatures,...
             model);


  save(outPath, 'model');

  err = 0;
end
