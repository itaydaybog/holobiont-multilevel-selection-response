Readme 

Author: Itay Daybog
itay.daybog@mail.huji.ac.il

The following Code and Data are supplied:

##############
# MODEL CODE #
##############

# mls_general_code_original.py #
VVD's original file.
Defines several helper functions used in the rest of the code

# MLS_static_fast_original.py #
VVD's original file.
Implementation of two species (helper and non-helper) microbiome Multilevel selection model

# MLS_static_fast_new.py #
An edited file.
Similar to VVD's original file, with the exception of the relation between the helper frequency and the fitness, now expressed as a step function.

###############
# Figure CODE #
###############

# MLS_figure_3_original.py #
VVD's original file.
Minor edit- the 3D graph was altered from a scatterplot to a surface plot to fit our interpolated discrete parameter space. 
Run to recreate figure 3 described in VVD's Manuscript. Code tries to load the existing data files from Data_Paper folder and uses it to make the figure. If the data file cannot be found, or if the model parameters have changed, the model is rerun. (warning: can take several hours to days!)

# MLS_figure_3_new_single_point.py #
An edited file.
Runs our altered model on a single point in the parameter space, designated by the global parameters "x_axis_parameter_ratio", "y_axis_parameter_ratio".
the model's run was split to single points in order to run all desired points on a cluster in parallel in order to shorten the long running time.

# build_3d_graph.py #
A novel file.
Uses all the model created data from the "Step_Data" folder, containing all the data created by the model run on the single points in the parameter space by the above file "MLS_figure_3_new_single_point.py", and displays it in a single figure as a whole parameter space.
2D heatmap can be left discrete or be interpolated, by changing the global boolean parameter "interpolate_2d_graph".

##############
# Data Files #
##############

# Original Data #
Data outputted by VVD's original model, by running "MLS_figure_3_original.py".

# Step_Data #
Data outputted by our edited model, by running "MLS_figure_3_new_single_point.py" on the following sets of points as grids- (X: [-2, -1, 0, 1, 2] on Y: [-2, -0.5, 1, 2.5, 4]) and (X: [-1.5, -0.5, 0.5, 1.5] on Y: [-1.25, 0.25, 1.75, 3.25]).

# Step_figures # & # Original_Figures #
Empty folders to contain figures created by our edited model and the original model respectively.



