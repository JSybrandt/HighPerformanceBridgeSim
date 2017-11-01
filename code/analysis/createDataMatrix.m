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
    global AllFrequencyData
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
    global AllFrequencyData
    global Time

    vec = [
        % Booleans
        Rain(dayIdx)

        % Scalars
        % L
        mu(dayIdx)
        %I
        %bbeta

        % Vectors
        % WindVelocity(Day, :)
        %ForceWind(dayIdx, :)
        Tact(dayIdx, pointIdx)
        VehicleMass(dayIdx, pointIdx)
        WheelMass(dayIdx, pointIdx)
        SuspensionStiffness(dayIdx, pointIdx)
        WheelStiffness(dayIdx, pointIdx)
        SuspensionDamping(dayIdx, pointIdx)
        WheelDamping(dayIdx, pointIdx)
        fv(dayIdx, pointIdx) % warning, this is still complex
        V(dayIdx, pointIdx)

        % Vector of Cells
        min(AllFrequencyData{dayIdx, pointIdx}(1,:))
        max(AllFrequencyData{dayIdx, pointIdx}(1,:))
        mean(AllFrequencyData{dayIdx, pointIdx}(1,:))
        min(AllFrequencyData{dayIdx, pointIdx}(2,:))
        max(AllFrequencyData{dayIdx, pointIdx}(2,:))
        mean(AllFrequencyData{dayIdx, pointIdx}(2,:))
        max(Time{dayIdx, pointIdx})
   ];
end
