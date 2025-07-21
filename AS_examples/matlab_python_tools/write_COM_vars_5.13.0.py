# trace generated using paraview version 5.13.0-RC1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'EnSight Reader'
ngacase = EnSightReader(registrationName='nga.case', CaseFileName='/mnt/RRG2_Backup/Andrew_backup/aerobreakup/AIAAWG_shot3/2D/500cpd_perturbed_vs_unperturbed/500cpd_PLICNET_perturbed_7-7-2025/ensight/ShockDroplet/nga.case')

# Properties modified on ngacase
ngacase.CellArrays = ['Lrho', 'VOF']

UpdatePipeline(time=0.0, proxy=ngacase)

# create a new 'Append Location Attributes'
appendLocationAttributes1 = AppendLocationAttributes(registrationName='AppendLocationAttributes1', Input=ngacase)

# Properties modified on appendLocationAttributes1
appendLocationAttributes1.AppendPointLocations = 0

UpdatePipeline(time=0.0, proxy=appendLocationAttributes1)

# create a new 'Cell Size'
cellSize1 = CellSize(registrationName='CellSize1', Input=appendLocationAttributes1)

# Properties modified on cellSize1
cellSize1.ComputeVertexCount = 0
cellSize1.ComputeLength = 0

UpdatePipeline(time=0.0, proxy=cellSize1)

# save data
SaveData('/mnt/RRG2_Backup/Andrew_backup/aerobreakup/AIAAWG_shot3/2D/500cpd_perturbed_vs_unperturbed/500cpd_PLICNET_perturbed_7-7-2025/trajectory_processing/data/COM_vars.csv', proxy=cellSize1, WriteTimeSteps=1,
    WriteTimeStepsSeparately=1,
    ChooseArraysToWrite=1,
    CellDataArrays=['Area', 'CellCenters', 'Lrho', 'VOF', 'Volume'],
    Precision=10,
    FieldAssociation='Cell Data',
    AddTime=1)