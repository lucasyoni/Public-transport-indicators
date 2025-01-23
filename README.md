[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/P3VlJKWK)
# Final Assignment

### Status

Once you are finished with the final assignment, edit this readme and add "x" to the correct box:

* [x] Submitted

* [ ] I'm still working on my final assignment. 

# Public Transport Network Service Level Indicators

The code in this assignment uses the following python packages:
- pathlib
- pandas
- geopandas
- matplotlib
- shapely
- contextily
- folium

*Note: Don't upload large files into GitHub! If you are using large input files, provide downloading instructions and perhaps a small sample of the data in this repository for demonstrating your workflow.*

## Topic: 

I am personally interested in urban geography and the geography of urban mobility. I was aware of GTFS (General Transit Feed Specifications) data before, but I had never had the chance or know-how to use it properly. The goal of this assignment was to code a rudimentary program that calculates descriptive indicators for public transport systems. The indicators reveal information on the spatial distribution and service level of public transport stops, which are like access points for travelers to take advantage of public transport. In this assignment, I analyzed the public transport networks of the cities of Helsinki, Stockholm and Copenhagen.

The analysis produces the following materials which can already be found in the "Final_materials" folder:
- a static map of rush hour service level in each city's public transport service area
- an interactive and a static map of rush hour service level per administrative/statistical region in each city's public transport service area
- a static map of stop density per square kilometer in each sity's public transport service area

The indicators allow the viewer to, at a glance, get a quick understanding of the breadth of each city's structure, the breadth of its public transport network, and the consistency of service during rush hour. 

In thhe future, I could continue this analysis by incorporating population data, which is likely available at a European scale. Being more familiar with Finnish sources of open geographical data, finding the different datasets for each country was surprisingly time consuming. Incorporating population data would reveal information not only about each publictransport system's geographicalk reach, but about how well it serves city-dwellers in each area. 

### Structure of this repository:

The code itself can be found in 2 different files:
- `pubtrans_indicators_workspace.ipynb`: a Jupyter notebook with markdown annotations
- - Each function is in its own cell. In the notebook all the functions are placed first. The last two cells contain code used for preparing the data and running the analysis. 
- `pubtrans_indicators.py`: a Python script that can be imported and used elsewhere. The end of the script includes an example of how the functions in the script can be run to perform the analysis. 

All data necessary for this analysis is found in the "data" folder. The data folder contains 3 city-specific folders, which each contain similarly named data. Fortunately there are consistent GTFS data naming conventions, so this code could be used for other public transport authorities (with some tweaks). 

THe results of the analysis (plots) can already be found in the "Final_materials" folder. 

### Input data:

#### GTFS Data

The input data is GTFS data from public transport authorities in 3 Nordic cities, Helsinki, Stockholm and Copenhagen (https://gtfs.org/documentation/schedule/reference/). Public transport authorities publish data about the services they provide, including schedules, stop and station locations, and routes. 

In the **Helsinki** metropolitan area, the local public transport authority Helsinki Region Transport (HSL) provides open access to their GTFS data on their own platform (https://www.hsl.fi/en/hsl/open-data). The municipal geometries are produced by Statistics Finland and access through a WFS interface. 

**Note: The Helsinki data is very extensive and the "stop_times.txt" file exceeds GitHub file size limits by a long shot.** I wrote this code locally on my computer in VSCode and the analysis is designed accordingly. I have uploaded an abbreviated version of the 2024 stop time data with the rush hour rows and necessary columns selected. I also only included every 10th row of the original set of data. However, this results in an inaccurate analysis. I would strongly recommend using the whole data from HSL's site. The results of the analysis performed with the complete data can be found in the "Final_materials" folder as well as in the Jupyter notebook. 

The **Stockholm** region GTFS data was retrieved from the Transit Feed Open mobility data platform. Transit feed is no longer in use and has been replaced by the Mobility Database catalog. However, the local transport authority SL (Stockholms Lokaltrafic) has not uploaded their GTFS data to the platform in an easily accessible way. For that reason, the data from Stockholm used in this analysis is from 2020. (https://transitfeeds.com/p/storstockholms-lokaltrafik/1086) 

The **Copenhagen** GTFS was retrieved using Rejseplanen Labs API (https://labs.rejseplanen.dk/hc/). I had to send a request to have an account made, so that I could access the data. The GTFS data is attached in the data folder in this repository. 

#### Geometries
The service area is determined based on geometries from administrative regions in the three metropolitan areaas. The administrative areas are available the Pan-European EuroRegionalMap by Open Maps for Europe. (https://www.mapsforeurope.org/datasets/euro-regional-map). However, in order to be able to perform spatial calculations, I went with geopackages from each country's own platforms which are set in each country's local CRS, which uses meters as units of measurement. In retrospect, this was definitely not the easiest way to go about it, but it only occurred to me to look for the Europe-wide administrative division once the code was already completed. If I were to do this assignment again, I would definitely use the European map data, and look into ways to automate the conversion of each dataset to its local meter-based CRS, which would also make applying this same code to different European cities more feasible. 

Finnish administrative regions (municipalities) were provided by STatistics Finland's WFS interface (https://geo.stat.fi/geoserver/wfs). HSL's GTFS data includes stops and routes from outside the Helsinki metropolitan area, so I manually selected the municipalities that could be considered more integral parts of the Helsinki region service area (even a few of those outside of HSL's jurisdiction). 

Swedish districts are from the Swedish Landsurveying Institute (Lantmäteriet) and were accessed from their website (https://www.lantmateriet.se/sv/geodata/vara-produkter/produktlista/administrativ-indelning-nedladdning-vektor/). The geopackage contains many small areas that cover the entirety of Sweden. I selected the areas from East and Central Sweden and included them in the data folder of this repository. 

Danish municipalities were sourced from Planinfo, the Danish Planning and Rural Development Agency's datashare platform's WFS interface (https://www.planinfo.dk/plandatadk/vejledning-til-plandatadk/plandatadks-webservices/wfs-web-feature-service). I again selected the municipalities considered part of the Copenhagen metro area manually. The Danish GTFS contains transit feed data from the entire country, and without otherwise being so familiar with Denmark's regional geography, it would be difficult to determine what data to include. 

### Analysis steps:
In this section I will describe the stages of the analysis performed in the code written for this assignment. The Jupyter notebook is also annotated with markdown cells to ease the understanding of the workflow. 

The code first goes through the GTFS data imported for each city. This analysis requires data about the following features:
- stops
- stop times
- trips
- calendar dates
- scheduling exceptions

* You will notice mention of a routes and line geometries in the code. I initially planned on incorporating the route geometries into the analysis, but as the code progressed, they became sidelined. Since the route data was very large and not used in the analysis (I initially considered using it for visualization), I decided to leave it out of the final submission. The lines of code that refer to it have been commented out (#). 

The code is set up to iterate over a list `cities`, which contains the names of the cities, whose data has been collected. As the data from each city is processed it is added to a dictionary, `cities_data`, where it can later be accessed by iterating over the same `cities` list. 

The data in this dictionary which is used for analysis are:
- stop locations
- all departures made from stops within the service area on weekdays** during rush hour (7-9 am and 4-6 pm)
- a grid of a resolution of the user's choosing
- the administrative region geometries that form the service area

** The Stockholm GTFS data did not contain data about weekday service, so it was not possible to filter departures by day of the week. 

Once the dictionary has been built, stop density and average rush hour departure amounts are calculated for each grid cell. To standardize the presentation of both stop density and rush hour departure amounts, I have classified them with predetermined categories. Stop densities are divided into 9 categories (from <1 to >100 stops / sqkm). Due to the fact that the GTFS data from each city did not have consistent dating (i.e. when the data was from), it was not really possible to filter the stop time and trip data to only represent a certain part of the year. In Helsinki, for example, public transport schedules are divided into winter and summer periods. Since I did not take potential seasonal variations into account, it was especially necessary to generalize the rush hour calculations. Rush hour departures are divided into quantiles. The values in the rush hour data represent "levels of service" where 1 equates to unfrequent or nonexistent rush hour service and 5 refers to frequent, reliable rush hour service. 

### Results:

The plots demostrate that public transport service, both in terms of the density of stops and the frequency of rush hour departures, is strongest in the urban core of all 3 cities and dissipates linearily along motorways to peripheral centers and suburbs. In both Helsinki and Stockholm there were strong similarities between the distribution of stops and the level of rush hour service. Areas with high density of stops also had high quality rush hour service. The data from Copenhagen showed quite a stark difference between stop density and rush hour service level. Almost all of the Copenhagen metropolitan area showed a high density of public transport spots, which might lead one to believe that the entire region is easily accessible by public transport. However, areas of high quality rush hour service are confined to small parts the absolute urban core and smaller centers along suburban rail routes. Most of the area is, based on the GTFS data analysis, rendered practically inaccessible during rush hour. IT can be implied that this fault in service level does not encourage Copenhageners to use public transport as a more sustainable alternative to automobility for their commute. Of course, this analysis alone cannot provide a full picture of the true state of things in the Copenhagen metro area. 

### References:

### Use of AI:

I used ChatGPT for small parts of the code presented in this assignment as well as for finding errors (print debugging was still more effective). As the code got longer, sifting it through ChatGPT proved an easy way to spot redundancies, typos, and other mistakes. The only parts of the code that are more heavily AI-assisted are the plotting functions, in which I wanted the composition of the different map, table and text elements to be more precise. Otherwise I referenced documentation pages, forums, and previous assignments from this course. 

### Literature

Urban Mobility Indicators for Walking and Public Transport, European Commission. https://ec.europa.eu/futurium/en/system/files/ged/convenient-access-to-public-transport.pdf 

Eboli, L. & Mazzulla, G. (2012) Performace indicators for an objective measure of public transport service quality. European Transport, 51:3. https://www.openstarts.units.it/server/api/core/bitstreams/19aa2bb3-6cb7-43db-9ae3-a81b6a660a9c/content 

### Feedback



Great job! You work is well structured -  you have taken clear steps to solve the problem and you have presented your results well, both verbally and visually. Your maps look nice. The documentation could be be more detailed.



### Points

| Item (max)                            | Points |
| ------------------------------------- | ------ |
| Reading data and prep.                | 10      |
| Data analysis                         | 10      |
| Visualization                         | 10      |
| Documentation                         | 8      |
| Coding style                          | 10      |
| Other merits (e.g., exceptional work) | 8      |
| **Total**                             | 56      |

## Course grade

- Final exercise: 56/60

- Weekly exercises: 38/40

- Total points: 94/100

- Grader: Kamyar
  
  

**Preliminary final grading thresholds:**

* 90% of possible points, or more → 5

* 80% of possible points, or more → 4

* 70% of possible points, or more → 3

* 60% of possible points, or more → 2

* 50% of possible points, or more → 1
