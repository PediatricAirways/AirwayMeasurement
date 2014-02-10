Pediatric Airway Atlas

How to use:
1. use script "create_directory" to make the directories for outputs
2. change the input and output directories in following two files:
    1) buildAirwayAtlasOnXSectionalArea.m 
    2) buildAirwayAtlasOnHydraulicDiameter.m
3. Call 1) to build atlas for cross sectional area or 2) for hydraulic diameter.

Atlas score:
The scores are saved to your_result_directory/Age/boxplotPart4/diagnosis_data/posInAtlas.mat
There are six columns in the file, score | age | score | age | score | age.
The point wise boxplot corresponds to the first two columns, 
the function boxplot corresponds to the next two columns, and
the weighted functional boxplot corresponds to the last two columns.
