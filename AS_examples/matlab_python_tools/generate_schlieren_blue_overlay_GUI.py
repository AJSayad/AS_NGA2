# trace generated using paraview version 5.13.0-RC1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

### user defined variables
case_path = '/mnt/RRG2_Backup/Andrew_backup/aerobreakup/AIAAWG_shot3/2D/100cpd_PLICNET_longrun_6-5-2025/ensight/ShockDroplet/nga.case'
image_name = '/mnt/RRG2_Backup/Andrew_backup/aerobreakup/AIAAWG_shot3/2D/100cpd_PLICNET_longrun_6-5-2025/schlieren_images/schlieren.png'
nt  = 148      # number of timesteps (data files)
d0  = 0.00187  # drop diameter
Lx  = 12*d0    # uniform domain length
Ly  = 4*d0     # uniform absolute domain height (top at Ly/2, bottom at -Ly/2)
tc2 = 1.47e-5  # characteristic breakup time

# STEP 1: Set background to white, hide xyz axis
# find settings proxy
generalSettings = GetSettingsProxy('GeneralSettings')

# find settings proxy
iOSettings = GetSettingsProxy('IOSettings')

# find settings proxy
renderViewInteractionSettings = GetSettingsProxy('RenderViewInteractionSettings')

# find settings proxy
renderViewSettings = GetSettingsProxy('RenderViewSettings')

# find settings proxy
representedArrayListSettings = GetSettingsProxy('RepresentedArrayListSettings')

# find settings proxy
colorPalette = GetSettingsProxy('ColorPalette')

# Properties modified on colorPalette
colorPalette.Background = [1.0, 1.0, 1.0]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# get the material library
materialLibrary1 = GetMaterialLibrary()

#STEP 2: load in simulation data (AS note: create a variable here for the ensight path)
# create a new 'EnSight Reader'
ngacase = EnSightReader(registrationName='nga.case', CaseFileName=case_path)

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on ngacase
ngacase.CellArrays = ['Density', 'VOF']

# show data in view
ngacaseDisplay = Show(ngacase, renderView1, 'UniformGridRepresentation')

# trace defaults for the display properties.
ngacaseDisplay.Representation = 'Outline'

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

# update the view to ensure updated data information
renderView1.Update()

#STEP 3: clip domain so that we only have the uniform region (AS note: create variables for domain using d0 and Lx,Ly)
# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=ngacase)

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [Lx, 0.0, 9.333333764516283e-06]

# show data in view
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'VOF'
vOFLUT = GetColorTransferFunction('VOF')

# get opacity transfer function/opacity map for 'VOF'
vOFPWF = GetOpacityTransferFunction('VOF')

# get 2D transfer function for 'VOF'
vOFTF2D = GetTransferFunction2D('VOF')

# create a new 'Clip'
clip2 = Clip(registrationName='Clip2', Input=clip1)

# Properties modified on clip2.ClipType
clip2.ClipType.Origin = [0.01122, Ly/2, 9.333333764516283e-06]
clip2.ClipType.Normal = [0.0, 1.0, 0.0]

# show data in view
clip2Display = Show(clip2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip2Display.Representation = 'Surface'

# hide data in view
Hide(clip1, renderView1)

# show color bar/color legend
clip2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Clip'
clip3 = Clip(registrationName='Clip3', Input=clip2)

# Properties modified on clip3.ClipType
clip3.ClipType.Origin = [0.01122, -Ly/2, 9.333333764516283e-06]
clip3.ClipType.Normal = [0.0, -1.0, 0.0]

# show data in view
clip3Display = Show(clip3, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip3Display.Representation = 'Surface'

# hide data in view
Hide(clip2, renderView1)

# show color bar/color legend
clip3Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

#STEP 4: liquid phase scaling factor
# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=clip3)

# Properties modified on calculator1
calculator1.ResultArrayName = 'Kliq'
calculator1.Function = '5'

# show data in view
calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'

# hide data in view
Hide(clip3, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'Kliq'
kliqLUT = GetColorTransferFunction('Kliq')

# get opacity transfer function/opacity map for 'Kliq'
kliqPWF = GetOpacityTransferFunction('Kliq')

# get 2D transfer function for 'Kliq'
kliqTF2D = GetTransferFunction2D('Kliq')

#STEP 6: gas phase scaling factor
# create a new 'Calculator'
calculator2 = Calculator(registrationName='Calculator2', Input=calculator1)

# Properties modified on calculator2
calculator2.ResultArrayName = 'Kgas'
calculator2.Function = '15'

# show data in view
calculator2Display = Show(calculator2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
calculator2Display.Representation = 'Surface'

# hide data in view
Hide(calculator1, renderView1)

# show color bar/color legend
calculator2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'Kgas'
kgasLUT = GetColorTransferFunction('Kgas')

# get opacity transfer function/opacity map for 'Kgas'
kgasPWF = GetOpacityTransferFunction('Kgas')

# get 2D transfer function for 'Kgas'
kgasTF2D = GetTransferFunction2D('Kgas')

#STEP 6: define beta function for schlieren
# create a new 'Calculator'
calculator3 = Calculator(registrationName='Calculator3', Input=calculator2)

# Properties modified on calculator3
calculator3.ResultArrayName = 'beta'
calculator3.Function = 'VOF*Kliq+Kgas*(1-VOF)'

# show data in view
calculator3Display = Show(calculator3, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
calculator3Display.Representation = 'Surface'

# hide data in view
Hide(calculator2, renderView1)

# show color bar/color legend
calculator3Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'beta'
betaLUT = GetColorTransferFunction('beta')

# get opacity transfer function/opacity map for 'beta'
betaPWF = GetOpacityTransferFunction('beta')

# get 2D transfer function for 'beta'
betaTF2D = GetTransferFunction2D('beta')

# set active source
SetActiveSource(calculator2)

#STEP 7: set annotation of gas phase scaling constant
# create a new 'Python Annotation'
pythonAnnotation1 = PythonAnnotation(registrationName='PythonAnnotation1', Input=calculator2)

# Properties modified on pythonAnnotation1
pythonAnnotation1.Expression = '"Kgas = "+ str(inputs[0].CellData[\'Kgas\'][0])'

# show data in view
pythonAnnotation1Display = Show(pythonAnnotation1, renderView1, 'TextSourceRepresentation')

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on pythonAnnotation1Display
pythonAnnotation1Display.WindowLocation = 'Any Location'

# Properties modified on pythonAnnotation1Display
pythonAnnotation1Display.Position = [0.65, 0.959558599695586]

# Properties modified on pythonAnnotation1Display
pythonAnnotation1Display.Position = [0.65, 0.25]

# Properties modified on pythonAnnotation1Display
pythonAnnotation1Display.Color = [0.0, 0.0, 0.0]

# Properties modified on pythonAnnotation1Display
pythonAnnotation1Display.Bold = 1

# set active source
SetActiveSource(calculator3)

#STEP 7: compute mixture density gradient
# create a new 'Gradient'
gradient1 = Gradient(registrationName='Gradient1', Input=calculator3)

# Properties modified on gradient1
gradient1.ScalarArray = ['CELLS', 'Density']
gradient1.ResultArrayName = 'grad_rho'

# show data in view
gradient1Display = Show(gradient1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
gradient1Display.Representation = 'Surface'

# hide data in view
Hide(calculator3, renderView1)

# show color bar/color legend
gradient1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

#STEP 8: compute magnitude of mixture density gradient
# create a new 'Calculator'
calculator4 = Calculator(registrationName='Calculator4', Input=gradient1)

# Properties modified on calculator4
calculator4.ResultArrayName = 'mag_grad_rho'
calculator4.Function = 'sqrt((grad_rho_X^2)+(grad_rho_Y^2)+(grad_rho_Z^2))'

# show data in view
calculator4Display = Show(calculator4, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
calculator4Display.Representation = 'Surface'

# hide data in view
Hide(gradient1, renderView1)

# show color bar/color legend
calculator4Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'mag_grad_rho'
mag_grad_rhoLUT = GetColorTransferFunction('mag_grad_rho')

# get opacity transfer function/opacity map for 'mag_grad_rho'
mag_grad_rhoPWF = GetOpacityTransferFunction('mag_grad_rho')

# get 2D transfer function for 'mag_grad_rho'
mag_grad_rhoTF2D = GetTransferFunction2D('mag_grad_rho')

#STEP 9: set constant scaling factor C
# create a new 'Calculator'
calculator5 = Calculator(registrationName='Calculator5', Input=calculator4)

# Properties modified on calculator5
calculator5.ResultArrayName = 'C'
calculator5.Function = '2.5e5'

# show data in view
calculator5Display = Show(calculator5, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
calculator5Display.Representation = 'Surface'

# hide data in view
Hide(calculator4, renderView1)

# show color bar/color legend
calculator5Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'C'
cLUT = GetColorTransferFunction('C')

# get opacity transfer function/opacity map for 'C'
cPWF = GetOpacityTransferFunction('C')

# get 2D transfer function for 'C'
cTF2D = GetTransferFunction2D('C')

#STEP 10: compute schlieren field
# create a new 'Calculator'
calculator6 = Calculator(registrationName='Calculator6', Input=calculator5)

# Properties modified on calculator6
calculator6.ResultArrayName = 'Schlieren'
calculator6.Function = 'exp(-beta*(abs(mag_grad_rho)/C))'

# show data in view
calculator6Display = Show(calculator6, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
calculator6Display.Representation = 'Surface'

# hide data in view
Hide(calculator5, renderView1)

# show color bar/color legend
calculator6Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'Schlieren'
schlierenLUT = GetColorTransferFunction('Schlieren')

# get opacity transfer function/opacity map for 'Schlieren'
schlierenPWF = GetOpacityTransferFunction('Schlieren')

# get 2D transfer function for 'Schlieren'
schlierenTF2D = GetTransferFunction2D('Schlieren')

#STEP 11: set properties for schlieren colormap and colorbar
# Properties modified on schlierenLUT
schlierenLUT.EnableOpacityMapping = 1

# Properties modified on schlierenPWF
schlierenPWF.Points = [0.0, 0.0892857164144516, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# Properties modified on schlierenPWF
schlierenPWF.Points = [0.0, 0.2901785969734192, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# Properties modified on schlierenPWF
schlierenPWF.Points = [0.0, 0.6116071939468384, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# Properties modified on schlierenPWF
schlierenPWF.Points = [0.0, 0.8392857313156128, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# Properties modified on schlierenPWF
schlierenPWF.Points = [0.0, 0.9464285969734192, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# Properties modified on schlierenPWF
schlierenPWF.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# Properties modified on schlierenPWF
schlierenPWF.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 0.8973214626312256, 0.5, 0.0]

# Properties modified on schlierenPWF
schlierenPWF.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 0.7455357313156128, 0.5, 0.0]

# Properties modified on schlierenPWF
schlierenPWF.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 0.504464328289032, 0.5, 0.0]

# Properties modified on schlierenPWF
schlierenPWF.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 0.1428571492433548, 0.5, 0.0]

# Properties modified on schlierenPWF
schlierenPWF.Points = [0.0, 1.0, 0.5, 0.0, 1.0, 0.0, 0.5, 0.0]

# get color legend/bar for schlierenLUT in view renderView1
schlierenLUTColorBar = GetScalarBar(schlierenLUT, renderView1)

# Properties modified on schlierenLUTColorBar
schlierenLUTColorBar.AutoOrient = 0
schlierenLUTColorBar.Orientation = 'Horizontal'
schlierenLUTColorBar.WindowLocation = 'Any Location'
schlierenLUTColorBar.Position = [0.35, 0.82]
schlierenLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
schlierenLUTColorBar.TitleBold = 1
schlierenLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
schlierenLUTColorBar.DrawScalarBarOutline = 1
schlierenLUTColorBar.ScalarBarOutlineColor = [0.0, 0.0, 0.0]
schlierenLUTColorBar.RangeLabelFormat = '%.1f'

# set active source
SetActiveSource(calculator5)

#STEP 12: create annotation for C constant
# create a new 'Python Annotation'
pythonAnnotation2 = PythonAnnotation(registrationName='PythonAnnotation2', Input=calculator5)

# Properties modified on pythonAnnotation2
pythonAnnotation2.Expression = '"C = "+ str(inputs[0].CellData[\'C\'][0])'

# show data in view
pythonAnnotation2Display = Show(pythonAnnotation2, renderView1, 'TextSourceRepresentation')

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on pythonAnnotation2Display
pythonAnnotation2Display.WindowLocation = 'Any Location'

# Properties modified on pythonAnnotation2Display
pythonAnnotation2Display.Position = [0.2, 0.959558599695586]

# Properties modified on pythonAnnotation2Display
pythonAnnotation2Display.Position = [0.2, 0.25]

# Properties modified on pythonAnnotation2Display
pythonAnnotation2Display.Color = [0.0, 0.0, 0.0]

# Properties modified on pythonAnnotation2Display
pythonAnnotation2Display.Bold = 1

# set active source
SetActiveSource(calculator6)

#STEP 13: create threshold for isolating the liquid phase
# create a new 'Threshold'
threshold1 = Threshold(registrationName='Threshold1', Input=calculator6)

# Properties modified on threshold1
threshold1.Scalars = ['CELLS', 'VOF']
threshold1.LowerThreshold = 0.95

# show data in view
threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'

# hide data in view
Hide(calculator6, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(calculator6)

# show data in view
calculator6Display = Show(calculator6, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
calculator6Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(threshold1)

#STEP 14: create a separate colormap for liqiud schlieren and set properties (blue overlay)
# change use separate color map
threshold1Display.UseSeparateColorMap = 1

# set scalar coloring using an separate color/opacity maps
ColorBy(threshold1Display, ('CELLS', 'Schlieren'), True)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(schlierenLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
threshold1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# get separate color transfer function/color map for 'Schlieren'
separate_threshold1Display_SchlierenLUT = GetColorTransferFunction('Schlieren', threshold1Display, separate=True)

# get separate opacity transfer function/opacity map for 'Schlieren'
separate_threshold1Display_SchlierenPWF = GetOpacityTransferFunction('Schlieren', threshold1Display, separate=True)

# get separate 2D transfer function for 'Schlieren'
separate_threshold1Display_SchlierenTF2D = GetTransferFunction2D('Schlieren', threshold1Display, separate=True)

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.EnableOpacityMapping = 1

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
separate_threshold1Display_SchlierenLUT.ApplyPreset('Blues', True)

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.031373, 0.188235, 0.419608, 0.12549, 0.031757, 0.318139, 0.612149, 0.1882355, 0.080969, 0.38113, 0.661361, 0.2509805, 0.130427, 0.444152, 0.710327, 0.3137255, 0.195386, 0.509112, 0.743791, 0.3764705, 0.260715, 0.573841, 0.777209, 0.4392155, 0.341423, 0.628958, 0.808704, 0.501960785, 0.422745, 0.684075, 0.839892, 0.564706, 0.523137, 0.739193, 0.861546, 0.627451, 0.622684, 0.793464, 0.883429, 0.690196, 0.701423, 0.826928, 0.910988, 0.7529410000000001, 0.778685, 0.8603, 0.937993, 0.8156865, 0.825928, 0.891795, 0.953741, 0.8784315, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.031373, 0.188235, 0.419608, 0.1882355, 0.080969, 0.38113, 0.661361, 0.2509805, 0.130427, 0.444152, 0.710327, 0.3137255, 0.195386, 0.509112, 0.743791, 0.3764705, 0.260715, 0.573841, 0.777209, 0.4392155, 0.341423, 0.628958, 0.808704, 0.501960785, 0.422745, 0.684075, 0.839892, 0.564706, 0.523137, 0.739193, 0.861546, 0.627451, 0.622684, 0.793464, 0.883429, 0.690196, 0.701423, 0.826928, 0.910988, 0.7529410000000001, 0.778685, 0.8603, 0.937993, 0.8156865, 0.825928, 0.891795, 0.953741, 0.8784315, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.031373, 0.188235, 0.419608, 0.3137255, 0.195386, 0.509112, 0.743791, 0.3764705, 0.260715, 0.573841, 0.777209, 0.4392155, 0.341423, 0.628958, 0.808704, 0.501960785, 0.422745, 0.684075, 0.839892, 0.564706, 0.523137, 0.739193, 0.861546, 0.627451, 0.622684, 0.793464, 0.883429, 0.690196, 0.701423, 0.826928, 0.910988, 0.7529410000000001, 0.778685, 0.8603, 0.937993, 0.8156865, 0.825928, 0.891795, 0.953741, 0.8784315, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.031373, 0.188235, 0.419608, 0.3764705, 0.260715, 0.573841, 0.777209, 0.4392155, 0.341423, 0.628958, 0.808704, 0.501960785, 0.422745, 0.684075, 0.839892, 0.564706, 0.523137, 0.739193, 0.861546, 0.627451, 0.622684, 0.793464, 0.883429, 0.690196, 0.701423, 0.826928, 0.910988, 0.7529410000000001, 0.778685, 0.8603, 0.937993, 0.8156865, 0.825928, 0.891795, 0.953741, 0.8784315, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.031373, 0.188235, 0.419608, 0.4392155, 0.341423, 0.628958, 0.808704, 0.501960785, 0.422745, 0.684075, 0.839892, 0.564706, 0.523137, 0.739193, 0.861546, 0.627451, 0.622684, 0.793464, 0.883429, 0.690196, 0.701423, 0.826928, 0.910988, 0.7529410000000001, 0.778685, 0.8603, 0.937993, 0.8156865, 0.825928, 0.891795, 0.953741, 0.8784315, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.031373, 0.188235, 0.419608, 0.564706, 0.523137, 0.739193, 0.861546, 0.627451, 0.622684, 0.793464, 0.883429, 0.690196, 0.701423, 0.826928, 0.910988, 0.7529410000000001, 0.778685, 0.8603, 0.937993, 0.8156865, 0.825928, 0.891795, 0.953741, 0.8784315, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.031373, 0.188235, 0.419608, 0.627451, 0.622684, 0.793464, 0.883429, 0.690196, 0.701423, 0.826928, 0.910988, 0.7529410000000001, 0.778685, 0.8603, 0.937993, 0.8156865, 0.825928, 0.891795, 0.953741, 0.8784315, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.031373, 0.188235, 0.419608, 0.690196, 0.701423, 0.826928, 0.910988, 0.7529410000000001, 0.778685, 0.8603, 0.937993, 0.8156865, 0.825928, 0.891795, 0.953741, 0.8784315, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.031373, 0.188235, 0.419608, 0.7529410000000001, 0.778685, 0.8603, 0.937993, 0.8156865, 0.825928, 0.891795, 0.953741, 0.8784315, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.031373, 0.188235, 0.419608, 0.8156865, 0.825928, 0.891795, 0.953741, 0.8784315, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.031373, 0.188235, 0.419608, 0.8784315, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.0, 0.188235, 0.419608, 0.8784315, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.0, 0.3, 0.419608, 0.8784315, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.0, 0.3, 1.0, 0.8784315, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.0, 0.3, 1.0, 0.33, 0.87328, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.0, 0.3, 1.0, 0.33, 0.0, 0.923291, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.0, 0.3, 1.0, 0.33, 0.0, 0.25, 0.969489, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.0, 0.3, 1.0, 0.33, 0.0, 0.25, 0.75, 0.9411765000000001, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.0, 0.3, 1.0, 0.33, 0.0, 0.25, 0.75, 0.66, 0.922491, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.0, 0.3, 1.0, 0.33, 0.0, 0.25, 0.75, 0.66, 0.0, 0.954787, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.0, 0.3, 1.0, 0.33, 0.0, 0.25, 0.75, 0.66, 0.0, 0.5, 0.985236, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.0, 0.3, 1.0, 0.33, 0.0, 0.25, 0.75, 0.66, 0.0, 0.5, 0.85, 1.0, 0.968627, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.0, 0.3, 1.0, 0.33, 0.0, 0.25, 0.75, 0.66, 0.0, 0.5, 0.85, 1.0, 0.7, 0.984314, 1.0]

# Properties modified on separate_threshold1Display_SchlierenLUT
separate_threshold1Display_SchlierenLUT.RGBPoints = [7.431936397955121e-248, 0.0, 0.3, 1.0, 0.33, 0.0, 0.25, 0.75, 0.66, 0.0, 0.5, 0.85, 1.0, 0.7, 0.9, 1.0]

# get color legend/bar for separate_threshold1Display_SchlierenLUT in view renderView1
separate_threshold1Display_SchlierenLUTColorBar = GetScalarBar(separate_threshold1Display_SchlierenLUT, renderView1)

# Properties modified on separate_threshold1Display_SchlierenLUTColorBar
separate_threshold1Display_SchlierenLUTColorBar.AutoOrient = 0
separate_threshold1Display_SchlierenLUTColorBar.Orientation = 'Horizontal'
separate_threshold1Display_SchlierenLUTColorBar.WindowLocation = 'Any Location'
separate_threshold1Display_SchlierenLUTColorBar.Position = [0.35, 0.8]
separate_threshold1Display_SchlierenLUTColorBar.Title = ''
separate_threshold1Display_SchlierenLUTColorBar.DrawScalarBarOutline = 1
separate_threshold1Display_SchlierenLUTColorBar.ScalarBarOutlineColor = [0.0, 0.0, 0.0]
separate_threshold1Display_SchlierenLUTColorBar.DrawTickMarks = 0
separate_threshold1Display_SchlierenLUTColorBar.DrawTickLabels = 0

# Properties modified on separate_threshold1Display_SchlierenLUTColorBar
separate_threshold1Display_SchlierenLUTColorBar.Position = [0.35, 0.79]

# Properties modified on separate_threshold1Display_SchlierenLUTColorBar
separate_threshold1Display_SchlierenLUTColorBar.AddRangeLabels = 0

# Properties modified on separate_threshold1Display_SchlierenLUTColorBar
separate_threshold1Display_SchlierenLUTColorBar.Position = [0.35, 0.795]

#STEP 15: create annotation for liquid phase scaling factor
# create a new 'Python Annotation'
pythonAnnotation3 = PythonAnnotation(registrationName='PythonAnnotation3', Input=threshold1)

# Properties modified on pythonAnnotation3
pythonAnnotation3.Expression = '"Kliq = " + str(inputs[0].CellData[\'Kliq\'][0])'

# show data in view
pythonAnnotation3Display = Show(pythonAnnotation3, renderView1, 'TextSourceRepresentation')

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on pythonAnnotation3Display
pythonAnnotation3Display.WindowLocation = 'Any Location'

# Properties modified on pythonAnnotation3Display
pythonAnnotation3Display.Position = [0.45, 0.959558599695586]

# Properties modified on pythonAnnotation3Display
pythonAnnotation3Display.Position = [0.45, 0.25]

# Properties modified on pythonAnnotation3Display
pythonAnnotation3Display.Color = [0.0, 0.0, 0.0]

# Properties modified on pythonAnnotation3Display
pythonAnnotation3Display.Bold = 1

# set active source
SetActiveSource(calculator6)

#STEP 16: create annotation for displaying time (AS note: create variable for scaling here 1/tc2)
# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(registrationName='AnnotateTimeFilter1', Input=calculator6)

# Properties modified on annotateTimeFilter1
annotateTimeFilter1.Format = 'Tau: {time:.2f}'
#annotateTimeFilter1.Scale = 67900.0
annotateTimeFilter1.Scale = 1/tc2

# show data in view
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1, 'TextSourceRepresentation')

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on annotateTimeFilter1Display
annotateTimeFilter1Display.WindowLocation = 'Any Location'

# Properties modified on annotateTimeFilter1Display
annotateTimeFilter1Display.Position = [0.45, 0.959558599695586]

# Properties modified on annotateTimeFilter1Display
annotateTimeFilter1Display.Position = [0.45, 0.74]

# Properties modified on annotateTimeFilter1Display
annotateTimeFilter1Display.Color = [0.0, 0.0, 0.0]

# Properties modified on annotateTimeFilter1Display
annotateTimeFilter1Display.Bold = 1

# set active source
SetActiveSource(calculator6)

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 1

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 0

# Properties modified on calculator6Display.DataAxesGrid
calculator6Display.DataAxesGrid.GridAxesVisibility = 1

# set active source
SetActiveSource(annotateTimeFilter1)

# toggle interactive widget visibility (only when running from the GUI)
ShowInteractiveWidgets(proxy=annotateTimeFilter1Display)

# toggle interactive widget visibility (only when running from the GUI)
HideInteractiveWidgets(proxy=annotateTimeFilter1Display)

# STEP 17: setup our axes grid (AS note: these should be normailzed by d0 in the datascale section)
# set active source
SetActiveSource(calculator6)

# Properties modified on calculator6Display.DataAxesGrid
calculator6Display.DataAxesGrid.GridAxesVisibility = 0

# Properties modified on calculator6Display.DataAxesGrid
calculator6Display.DataAxesGrid.GridAxesVisibility = 1

# hide data in view
Hide(ngacase, renderView1)

# Properties modified on calculator6Display.DataAxesGrid
calculator6Display.DataAxesGrid.GridAxesVisibility = 0

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.Visibility = 1

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.ZTitle = ''
renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XTitleBold = 1
renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YTitleBold = 1
renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.AxesToLabel = 3
renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.XLabelBold = 1
renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
renderView1.AxesGrid.YLabelBold = 1
renderView1.AxesGrid.DataScale = [0.00187, 0.00187, 1.0]
renderView1.AxesGrid.DataBoundsScaleFactor = 1.0

renderView1.ResetActiveCameraToNegativeZ()

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

# Properties modified on renderView1.AxesGrid
renderView1.AxesGrid.XTitle = 'X/d0'
renderView1.AxesGrid.YTitle = 'Y/d0'

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(808, 657)

# current camera placement for renderView1
renderView1.CameraPosition = [0.010330819214808142, -0.0011567720492031718, 0.03806556889379909]
renderView1.CameraFocalPoint = [0.010330819214808142, -0.0011567720492031718, -0.007630149236573997]
renderView1.CameraParallelScale = 0.011826922131777108

#STEP 18: save images (AS note: create a variable for image file name, image file path, and the number of timesteps(framewindow))
# save animation
SaveAnimation(filename=image_name, viewOrLayout=renderView1, location=16, ImageResolution=[808, 657],
    FrameWindow=[0, nt])

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(808, 657)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [0.010330819214808142, -0.0011567720492031718, 0.03806556889379909]
renderView1.CameraFocalPoint = [0.010330819214808142, -0.0011567720492031718, -0.007630149236573997]
renderView1.CameraParallelScale = 0.011826922131777108


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://kitware.github.io/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------
