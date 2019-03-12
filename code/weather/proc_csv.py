#!/usr/bin/env python3

import csv
from pathlib import Path
from pprint import pprint
from datetime import datetime, timedelta, date
import numpy as np
from scipy.io import savemat
import sys
from scipy.interpolate import interp1d

STATION_CODE="USW00094846"

def daterange(start_date, end_date):
  for n in range(int ((end_date - start_date).days)):
    yield start_date + timedelta(n)

hourly_avg_csv = Path("../../data/chicago_hourly_averages.csv")

assert hourly_avg_csv.is_file()

# year1 is a dummy year
start = date(year=1, month=1, day=1)
end = date(year=2, month=1, day=1)
date2hour_mean = {d:np.full((24,), np.nan) for d in daterange(start, end)}

with open(hourly_avg_csv) as f:
  for row in csv.DictReader(f):
    if row["STATION"] != STATION_CODE:
      continue
    try:
      temp = float(row["HLY-TEMP-NORMAL"])
    except:
      continue

    dt = datetime.strptime(row["DATE"], "%m-%dT%H:%M:%S")
    dt = dt.replace(year=1)
    if np.isnan(date2hour_mean[dt.date()][dt.hour]):
      date2hour_mean[dt.date()][dt.hour] = temp

highs_and_lows = Path("../../data/chicago_10-11.csv")
date2min = {}
date2max = {}
with open(highs_and_lows) as f:
  for row in csv.DictReader(f):
    if row["STATION"] != STATION_CODE:
      continue
    try:
      temp_min = float(row["TMIN"])
      temp_max = float(row["TMAX"])
    except:
      continue
    dt = datetime.strptime(row["DATE"], "%Y-%m-%d").date()
    date2min[dt] = temp_min
    date2max[dt] = temp_max

# two full years
start = date(year=2010, month=1, day=1)
end = date(year=2012, month=1, day=1)

num_days = len(date2min)
result = np.empty((num_days, 24))
# these are just to give us an initial go
yesterday_high = 13
yesterday_high_hour = -6
for day_idx, today in enumerate(daterange(start, end)):
  tomorrow = today + timedelta(1)

  # loop at end
  if tomorrow == end:
    tomorrow = start

  today_pattern = date2hour_mean[today.replace(year=1)]

  today_high_hour = np.argmax(today_pattern)
  today_low_hour = np.argmin(today_pattern)
  today_high = date2max[today]
  today_low = date2min[today]

  tomorrow_pattern = date2hour_mean[tomorrow.replace(year=1)]
  tomorrow_low = date2min[tomorrow]
  tomorrow_low_hour = np.argmin(tomorrow_pattern) + 24

  times = [yesterday_high_hour,
           today_low_hour,
           today_high_hour,
           tomorrow_low_hour]
  temps = [yesterday_high,
           today_low,
           today_high,
           tomorrow_low]
  inter = interp1d(times, temps)
  result[day_idx, :] = inter(range(24))
  yesterday_high = today_high
  yesterday_high_hour = today_high_hour-24



savemat("../../data/chicago_linear_interp_10-11.mat",
        {"temperature": result})



