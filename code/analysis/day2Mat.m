function res = day2Mat(day, path)

    NUM_PEAKS = 10;
    NUM_FEATURES = 2 * NUM_PEAKS + 6;


    filePath = strcat(path, '/', num2str(day), '.mat');
    load(filePath, 'VehicleFrequencyData', ...
                   'VehicleFrequencyAmpData', ...
                   'BridgeFrequencyData', ...
                   'V', ...
                   'DamageClass', ...
                   'VehicleMass');
                   % 'Tact', ...



    cars = size(VehicleFrequencyData, 2);

    res = zeros(cars, NUM_FEATURES);

    for car = 1:cars
        damageClass = DamageClass(day, car);
        mass = VehicleMass(day, car);
        speed = V(day, car);
        %temp = Tact(day, car);
        freq = VehicleFrequencyData{day, car};
        amp = VehicleFrequencyAmpData{day, car};
        bridgeFreq = BridgeFrequencyData{day, car}(1, :);
        minBF = min(bridgeFreq);
        maxBF = max(bridgeFreq);
        [peakAmp, peakFreq] = findpeaks(amp, freq);
        peaks = [peakFreq' peakAmp'];
        keyPeaks = peaks(1:NUM_PEAKS, :);

        res(car, :) = [keyPeaks(:)', ... % 2*NUM_PEAKS points
                       mass, ...
                       speed, ...  % after this, no more training feat
                       day, ...
                       minBF, ...
                       maxBF, ...
                       damageClass];
        %              temp, ...
    end

    disp(size(res))
end


