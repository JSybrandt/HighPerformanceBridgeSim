% This function is NOT designed to be called from command line
function result = day2Mat(input_path, record_both_cars)
  % In order to see progress, output the path
  disp(input_path)

  % Some constants to help the later calculations
  NUM_PEAKS = 10;
  NUM_FEATURES = 2 * NUM_PEAKS + 6;

  % Turn off load warnings, makes the progress clearer while running
  warning('off','all')

  % Load the day's data. By requesting each variable here, we get a clearer
  % error if something went wrong.
  load(input_path, 'Monitor_Vehicle_Frequency_Data', ...
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

  % This function is responsible for parsing a single car
  % Data that is not shared between primary and secondary cars is passed in
  function result = getVehicleData(car_idx, car_freqs, car_amps, car_masses)
    freq = car_freqs{Day, car_idx};
    amp = car_amps{Day, car_idx};
    [amp_peaks, freq_peaks] = findpeaks(amp, freq);
    peaks = [freq_peaks' amp_peaks'];

    avalible_peaks = size(peaks, 1);
    num_avalible_peaks = min(NUM_PEAKS, avalible_peaks);
    select_peaks = peaks(1:num_avalible_peaks, :);

    % We hit the case where there aren't enough peaks, pad with zero
    if num_avalible_peaks < NUM_PEAKS
      select_peaks = [select_peaks; zeros(NUM_PEAKS - num_avalible_peaks, 2)];
    end

    mass = car_masses(Day, car_idx);
    speed = V(Day, car_idx, 1);

    % Make sure to check that temperature effects actually exist
    if exist('Tact')
      temp = Tact(Day, car_idx);
    else
      temp = 0;
    end

    if exist('Rain')
      rain = Rain(Day);
    else
      rain = 0;
    end

    dc = DamageClass(Day, car_idx);

    % Warning, order matters, must equal the input order in python scripts
    result = [select_peaks(:)', ...
              mass, ...
              speed, ...
              temp, ...
              rain, ...
              Day, ...
              dc];
  end % getVehicleData

  % pre-allocate result, may over-allocate if Multiple_Vehicles
  if Multiple_Vehicles == 1
    result = zeros(num_cars * 2, NUM_FEATURES);
  else
    result = zeros(num_cars, NUM_FEATURES);
  else


  % track rows in result
  curr_idx = 1;
  for car_idx = 1:num_cars
    if Multiple_Vehicles
      for v = Order_of_Vehicles{Day, car_idx}
        if v == 1
          result(curr_idx,:) = getVehicleData(...
                                      car_idx, ...
                                      Monitor_Vehicle_Frequency_Data, ...
                                      Monitor_Vehicle_Frequency_Amp_Data, ...
                                      MonitorVehicleMass);
          curr_idx = curr_idx + 1;
        elseif record_both_cars
          result(curr_idx,:) = getVehicleData(...
                                      car_idx, ...
                                      Other_Vehicle_Frequency_Data, ...
                                      Other_Vehicle_Frequency_Amp_Data, ...
                                      SecondVehicleMass);
          curr_idx = curr_idx + 1;
        end % if
      end % for
    else % single vehicle
      result(curr_idx,:) = getVehicleData(...
                              car_idx, ...
                              Monitor_Vehicle_Frequency_Data, ...
                              Monitor_Vehicle_Frequency_Amp_Data, ...
                              MonitorVehicleMass);
      curr_idx = curr_idx + 1;
    end
  end % for each car

  % return only used rows
  result = result(1:curr_idx-1,:);

end % end function


