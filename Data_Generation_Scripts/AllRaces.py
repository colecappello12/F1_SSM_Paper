import fastf1 as f1
import pandas as pd
import os

# Enable cache (REQUIRED for performance)
#f1.Cache.enable_cache("fastf1_cache")

SEASON = 2025
DRIVER = "HAM"  # Lewis Hamilton
SAVE_INDIVIDUAL = False  # Set to False if you only want combined file

all_laps_list = []

# Get full event schedule minus testing
schedule = f1.get_event_schedule(SEASON)
schedule = schedule[schedule['Session5DateUtc'].notna()]

# Filter to official race sessions only
#race_events = schedule[schedule['EventFormat'].isin(['conventional', 'sprint'])]

for _, event in schedule.iterrows():

    event_name = event['EventName']
    print(f"Loading: {event_name}")

    try:
        session = f1.get_session(SEASON, event_name, "R")
        session.load()

        # Pick Lewis Hamilton laps
        laps = session.laps.pick_drivers(DRIVER).copy()

        # Drop in-laps and out-laps
        laps = laps[~laps['PitInTime'].notna()]
        laps = laps[~laps['PitOutTime'].notna()]

        # Convert LapTime to seconds
        laps = laps[laps['LapTime'].notna()]
        laps['LapTimeSeconds'] = laps['LapTime'].dt.total_seconds()

        # Add race identifier
        laps['Race'] = event_name
        laps['Year'] = SEASON
        laps['RoundNumber'] = event['RoundNumber']

        all_laps_list.append(laps)

        # Optional: Save each race separately
        if SAVE_INDIVIDUAL:
            filename = f"{event_name.replace(' ', '_')}_HAM_2025.csv"
            laps.to_csv(filename, index=False)

    except Exception as e:
        print(f"Skipping {event_name}: {e}")

# Combine all races into one dataframe
all_laps = pd.concat(all_laps_list, ignore_index=True)

# Save combined file
all_laps.to_csv("Hamilton_2025_All_Races.csv", index=False)

print("Done.")
