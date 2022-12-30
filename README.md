# Differential Expression Analysis Tools

Plot HeatMap

To run the above script, you can use the following command in the terminal:
``` 
Rscript plot_HeatMap.R --count_file count_file.txt --test_column_idx 1,2 --control_column_idx 3,4 --out_prefix HeatMap
```

This will run the script and generate a heatmap using the count file count_file.txt, with the test condition columns being indices 1 and 2, and the control condition columns being indices 3 and 4. The output heatmap will be saved as a PNG file with the filename HeatMap.png.

Note that you will need to have the required libraries (ggplot2, plotly, and argparse) installed in your R environment in order for the script to run properly. You can install these libraries by running the following commands in the R console:
```
install.packages("ggplot2")
install.packages("plotly")
install.packages("argparse")
```
