function out = createDataMatrix(dirPath, outPath)

PossibleDays = [1:365];
Days = [];

for day = PossibleDays
    if exist(getDayPath(dirPath, day), 'file') == 2
        Days = [Days ; day];
    end
end

globalLoad(dirPath, Days(1));
global Tact
numDays = size(Days, 1);
pointsPerDay = size(Tact,2);

exampleVec = getObservation(1, 1);

numFeatures = size(exampleVec, 1);

% each col is an observation, each row is a feature
% although we'll transform this when we're done
data = zeros(numFeatures, numDays * pointsPerDay);

for dayIdx = 1:numDays
    day = Days(dayIdx);
    startObIdx = (dayIdx - 1) * pointsPerDay + 1;
    endObIdx = startObIdx + pointsPerDay - 1;
    disp(day);
    globalLoad(dirPath, day);

    for pointIdx = 1:pointsPerDay
        data(:, (dayIdx - 1)*pointsPerDay + pointIdx) = ...
            getObservation(dayIdx, pointIdx);
    end

end

mins = min(data')'
maxes = max(data')'
diffs = maxes - mins;

data = (data - mins) ./ diffs;

data = data';
mins = mins';
maxes = maxes';

save(outPath, 'data', 'mins', 'maxes')

out = 0

end

function ret = getDayPath(rootPath, day)
    ret = fullfile(rootPath, strcat(int2str(day), '.mat'));
end

function globalLoad(rootPath, day)
    global Rain
    global L
    global mu
    global I
    global bbeta
    global WindVelocity
    global ForceWind
    global Tact
    global VehicleMass
    global WheelMass
    global SuspensionStiffness
    global WheelStiffness
    global SuspensionDamping
    global WheelDamping
    global fv
    global V
    global VehicleFrequencyData
    global RoadProfile
    global Time
    load(getDayPath(rootPath, day));
end

function vec = getObservation(dayIdx, pointIdx)
    global Rain
    global L
    global mu
    global I
    global bbeta
    global WindVelocity
    global ForceWind
    global Tact
    global VehicleMass
    global WheelMass
    global SuspensionStiffness
    global WheelStiffness
    global SuspensionDamping
    global WheelDamping
    global fv
    global V
    global VehicleFrequencyData
    global RoadProfile
    global Time

    vec = [
        % Booleans
        Rain(dayIdx)                                            % 1

        % Scalars
        % L
        mu(dayIdx)                                              % 2
        %I
        %bbeta

        % Vectors
        % WindVelocity(Day, :)
        %ForceWind(dayIdx, :)
        Tact(dayIdx, pointIdx)                                  % 3
        VehicleMass(dayIdx, pointIdx)                           % 4
        WheelMass(dayIdx, pointIdx)                             % 5
        SuspensionStiffness(dayIdx, pointIdx)                   % 6
        WheelStiffness(dayIdx, pointIdx)                        % 7
        SuspensionDamping(dayIdx, pointIdx)                     % 8
        WheelDamping(dayIdx, pointIdx)                          % 9
        fv(dayIdx, pointIdx) % warning, this is still complex   % 10
        V(dayIdx, pointIdx)                                     % 11

        % Vector of Cells
        RoadProfile{dayIdx, pointIdx} % results in a vector of observations

        max(Time{dayIdx, pointIdx})                             %
   ];
end
