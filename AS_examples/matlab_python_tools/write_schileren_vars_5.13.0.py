# trace generated using paraview version 5.13.0-RC1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'EnSight Reader'
ngacase = EnSightReader(registrationName='nga.case', CaseFileName='/mnt/RRG2_Backup/Andrew_backup/aerobreakup/AIAAWG_shot3/2D/100cpd_PLICNET_longrun_6-5-2025/ensight/ShockDroplet/nga.case')

# Properties modified on ngacase
ngacase.CellArrays = ['Density', 'VOF']

UpdatePipeline(time=0.0, proxy=ngacase)

# save data
SaveData('/mnt/RRG2_Backup/Andrew_backup/aerobreakup/AIAAWG_shot3/2D/100cpd_PLICNET_longrun_6-5-2025/temp2/schlieren_data.csv', proxy=ngacase, WriteTimeSteps=1,
    WriteTimeStepsSeparately=1,
    ChooseArraysToWrite=1,
    CellDataArrays=['Density', 'VOF'],
    Precision=10,
    FieldAssociation='Cell Data',
    AddTime=1)