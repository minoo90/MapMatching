# MapMatching: A method for matching vehicles GPS data on the digital roads network map

The map matching uses both digital road network (DRN) topology and spatial relations. It matches GPS tracks to roads. In this work, OpenStreetMap digital road network of City Seattle extracted and used for DRN. Trajectory dataset is a vehicle’s GPS data that has latitude, longitude, time, speed and heading. The sampling rate is every second of its travel during a day. 
Three main functions are used in this work, Initial road and location detection (IMP), follow up a vehicle alongside the current road (SMP1), detect and locate the vehicle on a road when it cross the intersection(SMP2). 
Different sets of fuzzy logic rules is needed for each aforementioned function to gain an optimized speed for the algorithm. Therefore, corresponding zero-order Sugeno Fuzzy Inference System (FIS) is developed for each function. 
The code is implemented in Rstudio 3.4.
