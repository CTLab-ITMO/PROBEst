# Exec:
# bash grid_search.sh

# Required:
# R > 4.2 with libs:
# shiny
# ggplot2
# plotly
# dplyr

#scripts/generator/
python test_parameters.py \
    -p data/test/general/param_grid_light.json

#R -e 'shiny::runApp("scripts/generator/grid_viz.R")'

Rscript scripts/generator/grid_viz.R data/param_search/param_stats.csv