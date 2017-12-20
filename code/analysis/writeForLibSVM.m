function err = trainModel(dataPath, trainingPath, validationPath, featureRows, labelRow)
  % Expected to load data matrix that is 0-1 normalized
  load(dataPath);

  rng('shuffle')

  numFolds = 10;

  featureRows = str2num(featureRows);
  labelRow = str2num(labelRow);

  numTuples = size(data, 2);
  trainingSetSize = floor(numTuples*((numFolds-1)/numFolds));
  shuffle = randperm(numTuples);

  data = data(:, shuffle);

%  for fold = 1:numFolds

%    disp('Fold:')
%    disp(fold)

  trainingSet = data(:, 1:trainingSetSize);
  validationSet = data(:, trainingSetSize+1:numTuples);

  disp('Training Set Size')
  size(trainingSet)
  disp('Validation Set Size')
  size(validationSet)

  trainingFeatures = trainingSet(featureRows,:)';
  trainingLabels = trainingSet(labelRow, :)';


%  model = svmtrain(trainingLabels,...
%                   trainingFeatures,...
%                   '-s 3 -t 2 -h 0');
%  disp('TRAINING SET:')
%  svmpredict(trainingLabels,...
%             trainingFeatures,...
%             model);
%
  validationFeatures = validationSet(featureRows, :)';
  validationLables = validationSet(labelRow, :)';
%  disp('Validation SET:')
%  svmpredict(validationLables,...
%             validationFeatures,...
%             model);
%
%    data = [data(:, trainingSetSize+1:numTuples) data(:, 1:trainingSetSize)];

%  end

%  save(outPath, 'model');

  libsvmwrite(trainingPath, trainingLabels, sparse(trainingFeatures))
  libsvmwrite(validationPath, validationLables, sparse(validationFeatures))

  err = 0;
end
