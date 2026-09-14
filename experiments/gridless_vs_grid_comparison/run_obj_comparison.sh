#!/bin/bash

# Run the shell script from experiments/gridless_vs_grid_comparison
# Move up to the parent directory
cd ../..

for s in {1..13}
do
   echo "--------------------------------------"
   echo "Starting Setting $s at $(date)"
   
   # Reset the timer
   START_TIME=$SECONDS
   
   # Run the R script
   Rscript experiments/gridless_vs_grid_comparison/time_and_objvalues.R $s 100
   
   # Calculate elapsed time
   ELAPSED_TIME=$(($SECONDS - $START_TIME))
   
   echo "Finished Setting $s"
   echo "Elapsed time: $(($ELAPSED_TIME / 60)) minutes and $(($ELAPSED_TIME % 60)) seconds."
done
