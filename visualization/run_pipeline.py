# Test run of VTK-logger by MW - Adapted from Tomek

# ### INSTRUCTIONS ###
"""
1. Set the directories with which the script should work:
"""
# Root directory of the eGFRD algorithm files:
# Use the full path (e.g. NOT a tilde)
egfrd_directory ='/storage1/wehrens/myfork_egfrd/'
# Output directory of the VTK logger 
# (This was defined when calling the logger in your
# script.)
simulation_data_directory = egfrd_directory + \
    'samples/memory/VTK_out/'
"""

2. Start ParaView. Make sure you give Paraview the PYTHONPATH 
   information:

    This means, in bash:
    $ export PYTHONPATH=/usr/lib/paraview/
    $ paraview

3. In ParaView, select tools > python shell > run script. Then
   select this script. If everything went well, you've now 
   imported the data correctly in ParaView.

"""



# paraview_directory = '/usr/lib/paraview/'

# Make sure python knows where to find eGFRD files
import sys
sys.path.append(egfrd_directory)
# sys.path.append(paraview_directory)


# import the pipeline script - which converts the VTKlogger data to data
# Paraview can understand.
from visualization import pipeline

# Create pipeline.
p = pipeline.Pipeline(egfrd_directory, simulation_data_directory)

# Script constants
p.particle_radius_scale_factor = 1
p.helix_radius_scale_factor = 1
p.cylinder_radius_scale_factor = 1

# Go!
p.rebuild()
