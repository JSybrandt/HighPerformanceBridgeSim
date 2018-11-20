function res = day2Mat(path, use_other)
    disp(path)

    NUM_PEAKS = 10;
    NUM_FEATURES = 2 * NUM_PEAKS + 6;
    warning('off','all')
    load(path, 'Monitor_Vehicle_Frequency_Data', ...
               'Monitor_Vehicle_Frequency_Amp_Data', ...
               'V', ...
               'Tact', ...
               'Rain', ...
               'DamageClass', ...
               'MonitorVehicleMass',...
               'Day',...
               'Multiple_Vehicles', ...
               'n',...
               'Order_of_Vehicles',...
               'Other_Vehicle_Frequency_Data',...
               'Other_Vehicle_Frequency_Amp_Data',...
               'SecondVehicleMass', ...
               'DayDamage1', ...
               'DayDamage2', ...
               'DayDamage3', ...
               'DayDamage4', ...
               'DayDamage5', ...
               'Damage_Case');

    function dc = getDamageClass(day)
      % if Damage_Case == 2
        if day >=  DayDamage5
          dc = 5;
        elseif day >= DayDamage4
          dc = 4;
        elseif day >= DayDamage3
          dc = 3;
        elseif day >= DayDamage2
          dc = 2;
        elseif day >= DayDamage1
          dc = 1;
        else
          dc = 0;
        end
      % else
        % dc = DamageClass(Day, car);
      % end
    end


    function data = getMonitorData(car)
        freq = Monitor_Vehicle_Frequency_Data{Day, car};
        amp = Monitor_Vehicle_Frequency_Amp_Data{Day, car};
        [peakAmp, peakFreq] = findpeaks(amp, freq);
        peaks = [peakFreq' peakAmp'];

        avalible_peaks = size(peaks, 1);
        select_peaks = min(NUM_PEAKS, avalible_peaks);
        keyPeaks = peaks(1:select_peaks, :);

        % We hit the case where there aren't enough peaks
        if select_peaks < NUM_PEAKS
            keyPeaks = [keyPeaks; zeros(NUM_PEAKS - select_peaks, 2)];
        end

        mass = MonitorVehicleMass(Day, car);
        speed = V(Day, car, 1);

        % these don't exist for different day types
        if exist('Tact')
            temp = Tact(Day, car);
        else
            temp = 0;
        end

        if exist('Rain')
            rain = Rain(Day);
        else
            rain = 0;
        end

        % damageClass = DamageClass(Day, car);
        damageClass = getDamageClass(Day);

        data = [keyPeaks(:)', ...
                mass, ...
                speed, ...
                temp, ...
                rain, ...
                Day, ...
                damageClass];
    end

    function data = getOtherData(car)
        freq = Other_Vehicle_Frequency_Data{Day, car};
        amp = Other_Vehicle_Frequency_Amp_Data{Day, car};
        [peakAmp, peakFreq] = findpeaks(amp, freq);
        peaks = [peakFreq' peakAmp'];

        avalible_peaks = size(peaks, 1);
        select_peaks = min(NUM_PEAKS, avalible_peaks);
        keyPeaks = peaks(1:select_peaks, :);

        % We hit the case where there aren't enough peaks
        if select_peaks < NUM_PEAKS
            keyPeaks = [keyPeaks; zeros(NUM_PEAKS - select_peaks, 2)];
        end

        mass = SecondVehicleMass(Day, car);
        speed = V(Day, car, 2);

        % these don't exist for different day types
        if exist('Tact')
            temp = Tact(Day, car);
        else
            temp = 0;
        end

        if exist('Rain')
            rain = Rain(Day);
        else
            rain = 0;
        end

        damageClass = getDamageClass(Day);

        data = [keyPeaks(:)', ...
                mass, ...
                speed, ...
                temp, ...
                rain, ...
                Day, ...
                damageClass];
    end

    num_cars = n;

    if Multiple_Vehicles == 1
        res = zeros(num_cars * 2, NUM_FEATURES);
    else
        res = zeros(num_cars, NUM_FEATURES);
    end


    curr_idx = 1;
    for car = 1:num_cars
        if Multiple_Vehicles
            for v = Order_of_Vehicles{Day, car}
                if v == 1
                    res(curr_idx,:) = getMonitorData(car);
                    curr_idx = curr_idx + 1;
                elseif use_other
                    res(curr_idx,:) = getOtherData(car);
                    curr_idx = curr_idx + 1;
                end
            end
        else
            res(curr_idx,:) = getMonitorData(car);
            curr_idx = curr_idx + 1;
        end
    end

    % return only used rows
    res = res(1:curr_idx-1,:);
end


