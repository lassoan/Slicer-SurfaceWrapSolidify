import os
import vtk, qt, ctk, slicer
import logging
from SegmentEditorEffects import *
import numpy as np
import math
import vtkSegmentationCorePython

class SegmentEditorEffect(AbstractScriptedSegmentEditorEffect):
  """This effect uses shrinkwrap, raycasting, remesh, and solidifying algorithms to filter the surface from the input segmentation"""

  def __init__(self, scriptedEffect):
    scriptedEffect.name = 'Wrap Solidify'
    scriptedEffect.perSegment = True # this effect operates on all segments at once (not on a single selected segment)
    AbstractScriptedSegmentEditorEffect.__init__(self, scriptedEffect)

    self.logic = WrapSolidifyLogic()
    self.logic.logCallback = self.addLog

  def clone(self):
    # It should not be necessary to modify this method
    import qSlicerSegmentationsEditorEffectsPythonQt as effects
    clonedEffect = effects.qSlicerSegmentEditorScriptedEffect(None)
    clonedEffect.setPythonSource(__file__.replace('\\','/'))
    return clonedEffect

  def icon(self):
    # It should not be necessary to modify this method
    iconPath = os.path.join(os.path.dirname(__file__), 'SegmentEditorEffect.png')
    if os.path.exists(iconPath):
      return qt.QIcon(iconPath)
    return qt.QIcon()

  def helpText(self):
    return """<html>Create a solid segment from the outer surface or an internal surface of a segment. It is using a combination of shrinkwrapping, projection and solidification algorithms.<br>
    For further information, license, disclaimers and possible research partnerships visit <a href="https://github.com/sebastianandress/Slicer-SurfaceWrapSolidify">this</a> github repository.
    </html>"""

  def activate(self):
    pass

  def deactivate(self):
    self.cleanup()
  
  def cleanup(self):
    pass

  def setupOptionsFrame(self):

    # Load widget from .ui file. This .ui file can be edited using Qt Designer
    # (Edit / Application Settings / Developer / Qt Designer -> launch).
    uiWidget = slicer.util.loadUI(os.path.join(os.path.dirname(__file__), "SegmentEditorEffect.ui"))
    self.scriptedEffect.addOptionsWidget(uiWidget)
    self.ui = slicer.util.childWidgetVariables(uiWidget)

    # Set scene in MRML widgets. Make sure that in Qt designer
    # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
    # "setMRMLScene(vtkMRMLScene*)" slot.
    uiWidget.setMRMLScene(slicer.mrmlScene)

    self.ui.modeGroup = qt.QButtonGroup()
    self.ui.modeGroup.addButton(self.ui.modeOuterSurfaceRadioButton)
    self.ui.modeGroup.addButton(self.ui.modeInnerSurfaceRadioButton)

    self.ui.regionGroup = qt.QButtonGroup()
    self.ui.regionGroup.addButton(self.ui.regionAutoRadioButton)
    self.ui.regionGroup.addButton(self.ui.regionManualRadioButton)

    self.outputTypeGroup = qt.QButtonGroup()
    self.outputTypeGroup.addButton(self.ui.outputSegmentationRadioButton)
    self.outputTypeGroup.addButton(self.ui.outputModelRadioButton)

    # Widget to arguments mapping
    self.valueEditWidgets = {
      ARG_MODE: self.modeGroup,
      ARG_REGION: self.regionGroup,
      ARG_REGION_SEGMENT_ID: self.ui.regionSegmentSelector,
      ARG_REGION_EXTEND: self.ui.extendRegionCheckBox,
      ARG_REGION_EXTEND_DIAMETER, self.ui.regionExtendDiameterSlider,
      ARG_CREATE_SHELL: self.ui.createShellCheckBox,
      ARG_SHELL_THICKNESS: self.ui.shellThicknessSlider,
      ARG_OUTPUT_TYPE: self.outputTypeGroup,
      ARG_OUTPUT_MODEL_NODE: self.outputModelNodeSelector
      ARG_SMOOTHING_FACTOR: self.ui.smoothingFactorSlider,
      ARG_REMESH_OVERSAMPLING: self.ui.remeshOversamplingSlider,
      ARG_SHRINKWRAP_ITERATIONS: self.ui.iterationsSlider,
      ARG_SAVE_INTERMEDIATE_RESULTS: saveIntermediateResultsCheckBox
    }

    # Add connections

    for argName, widget in self.valueEditWidgets.items():
      widgetClassName = widget.metaObject().getClassName()
      if widgetClassName=="qMRMLSliderWidget":
        widget.connect("valueChanged(double)", self.updateMRMLFromGUI)
      elif widgetClassName=="QCheckBox":
        widget.connect("stateChanged(int)", self.updateMRMLFromGUI)
      elif widgetClassName=="qMRMLNodeComboBox":
        widget.connect("currentNodeChanged(vtkMRMLNode*)", self.updateMRMLFromGUI)
      elif widgetClassName=="QButtonGroup":
        widget.connect("buttonClicked(int)", self.updateMRMLFromGUI)
      elif widgetClassName=="qMRMLSegmentSelectorWidget":
        widget.connect("currentSegmentChanged(QString)", self.updateMRMLFromGUI)

    self.ui.applyButton.connect('clicked()', self.onApply)

  def createCursor(self, widget):
    return slicer.util.mainWindow().cursor

  def layoutChanged(self):
    pass

  def processInteractionEvents(self, callerInteractor, eventId, viewWidget):
    return False # For the sake of example

  def processViewNodeEvents(self, callerViewNode, eventId, viewWidget):
    pass # For the sake of example

  def setMRMLDefaults(self):
    for (argName, defaultValue) in ARG_DEFAULTS:
      self.scriptedEffect.setParameterDefault(argName, defaultValue)

  def updateGUIFromMRML(self):
    parameterNode = self.getParameterNode()
    if not parameterNode:
      return

    segmentationNode = self.scriptedEffect.parameterSetNode().GetSegmentationNode().GetSegmentation()

    # Update values in widgets
    for argName, widget in self.valueEditWidgets.items():
      widgetClassName = widget.metaObject().getClassName()
      oldBlockSignalsState = widget.blockSignals(True)
      if widgetClassName=="qMRMLSliderWidget":
        widget.value = self.scriptedEffect.doubleParameter(argName)
      elif widgetClassName=="QCheckBox":
        widget.setChecked(self.scriptedEffect.parameter(argName)=='True')
      elif widgetClassName=="qMRMLNodeComboBox":
        widget.setCurrentNodeID(parameterNode.GetNodeReferenceID(parameterName))
      elif widgetClassName=="QButtonGroup":
        optionIndex = ARG_OPTIONS[argName].index(self.scriptedEffect.parameter(argName))
        widget.button(-2-optionIndex).setChecked(True)
      elif widgetClassName=="qMRMLSegmentSelectorWidget":
        segmentId = self.scriptedEffect.parameter(argName)
        if widget.currentNode() != segmentationNode:
          widget.setCurrentNode(segmentationNode)
        widget.setCurrentSegmentID(segmentId)
      else:
        raise Exception("Unexpected widget class: {0}".format(widgetClassName))
      widget.blockSignals(oldBlockSignalsState)

    # Enable/disable dependent widgets
    regionExtend = self.scriptedEffect.parameter(ARG_REGION_EXTEND)
    self.valueEditWidgets[ARG_REGION_EXTEND_DIAMETER].enabled = regionExtend
    createShell = self.scriptedEffect.parameter(ARG_CREATE_SHELL)
    self.valueEditWidgets[ARG_SHELL_THICKNESS].enabled = createShell
    self.valueEditWidgets[ARG_SHELL_PRESERVE_CRACKS].enabled = createShell

    # Update output section
    region = self.scriptedEffect.parameter(ARG_REGION)
    if region == REGION_SEGMENT:
      outputSegmentId = self.scriptedEffect.parameterSetNode().GetSelectedSegmentID()
      self.ui.outputSegmentNameLabel.text = segmentationNode.GetSegment(outputSegmentId).GetName()
    else:  # inner surface
      outputSegmentId = self.scriptedEffect.parameter(ARG_REGION_SEGMENT_ID)
      if outputSegmentId == self.scriptedEffect.parameterSetNode().GetSelectedSegmentID():
        self.ui.outputSegmentNameLabel.text = "error - region segment cannot be the same as the current segment"
      else:
        self.ui.outputSegmentNameLabel.text = segmentationNode.GetSegment(outputSegmentId).GetName()

  def updateMRMLFromGUI(self):
    wasModified = self.scriptedEffect.parameterSetNode().StartModify()
    for argName, widget in self.valueEditWidgets.items():
      widgetClassName = widget.metaObject().getClassName()
      if widgetClassName=="qMRMLSliderWidget":
        self.scriptedEffect.setParameter(argName, widget.value)
      elif widgetClassName=="QCheckBox":
        self.scriptedEffect.setParameter(argName, "True" if widget.isChecked() else "False")
      elif widgetClassName=="qMRMLNodeComboBox":
        self.scriptedEffect.parameterSetNode().SetNodeReferenceID(argName, widget.currentNodeID)
      elif widgetClassName=="QButtonGroup":
        optionName = ARG_OPTIONS[argName][-2-widget.checkedId()]
        self.scriptedEffect.setParameter(argName, optionName)
      elif widgetClassName=="qMRMLSegmentSelectorWidget":
        segmentId = self.scriptedEffect.parameter(argName)
        if widget.currentNode() != segmentationNode:
          widget.setCurrentNode(segmentationNode)
        widget.setCurrentSegmentID(segmentId)
      else:
        raise Exception("Unexpected widget class: {0}".format(widgetClassName))
    self.scriptedEffect.parameterSetNode().EndModify(wasModified)

  def addLog(self, text):
    slicer.util.showStatusMessage(text)
    slicer.app.processEvents() # force update

  def onApply(self):

    if self.ui.applyButton.text == 'Cancel':
      self.logic.requestCancel()
      return

    self.scriptedEffect.saveStateForUndo()

    errorMessage = None
    self.ui.applyButton.text = 'Cancel'
    qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
    try:
      self.logic.segmentationNode = self.scriptedEffect.parameterSetNode().GetSegmentationNode()
      self.logic.segmentId = currentSegmentId = self.scriptedEffect.parameterSetNode().GetSelectedSegmentID()
      self.logic.mode = self.scriptedEffect.parameter(ARG_MODE)
      self.logic.region = self.scriptedEffect.parameter(ARG_REGION)
      self.logic.regionSegmentId = self.scriptedEffect.parameter(ARG_REGION_SEGMENT_ID)
      self.logic.regionExtend = (self.scriptedEffect.parameter(ARG_REGION_EXTEND] == "True")
      self.logic.regionExtendDiameter = self.scriptedEffect.doubleParameter(ARG_REGION_EXTEND_DIAMETER)
      self.logic.createShell = (self.scriptedEffect.parameter(ARG_CREATE_SHELL] == "True")
      self.logic.shellThickness = self.scriptedEffect.doubleParameter(ARG_SHELL_THICKNESS)
      self.logic.shellPreserveCracks = (self.scriptedEffect.parameter(ARG_SHELL_PRESERVE_CRACKS] == "True")
      self.logic.outputType = self.scriptedEffect.parameter(ARG_OUTPUT_TYPE)
      self.logic.outputModelNode = self.scriptedEffect.parameterSetNode().GetNodeReferenceID(ARG_OUTPUT_MODEL_NODE)
      self.logic.remeshOversampling = self.scriptedEffect.doubleParameter(ARG_REMESH_OVERSAMPLING)
      self.logic.smoothingFactor = self.scriptedEffect.doubleParameter(ARG_SMOOTHINGFACTOR)
      self.logic.shrinkwrapIterations = self.scriptedEffect.integerParameter(ARG_SHRINKWRAP_ITERATIONS)
      self.logic.saveIntermediateResults = (self.scriptedEffect.parameter(ARG_SAVE_INTERMEDIATE_RESULTS) == "True")
      self.logic.applyWrapSolidify():
    except Exception as e:
      import traceback
      traceback.print_exc()
      errorMessage = str(e)
    self.ui.applyButton.text = 'Apply'
    qt.QApplication.restoreOverrideCursor()
    if errorMessage:
      slicer.util.errorDisplay("Wrap solidfy failed: " + self.logic.resultMessage)


class WrapSolidifyLogic(object):

  def __init__(self):
    self.logCallback = None
    self.cancelRequested = False

    # Inputs
    self.segmentationNode = None
    self.segmentId = None
    self.mode = ARG_DEFAULTS[ARG_MODE]
    self.region = ARG_DEFAULTS[ARG_REGION]
    self.regionSegmentId = None
    self.regionExtend = ARG_DEFAULTS[ARG_REGION_EXTEND]
    self.regionExtendDiameter = ARG_DEFAULTS[ARG_REGION_EXTEND_DIAMETER]
    self.createShell = ARG_DEFAULTS[ARG_CREATE_SHELL]
    self.shellThickness = ARG_DEFAULTS[ARG_SHELL_THICKNESS]
    self.shellPreserveCracks = ARG_DEFAULTS[ARG_SHELL_PRESERVE_CRACKS]
    self.outputType = ARG_DEFAULTS[ARG_OUTPUT_TYPE]
    self.outputModelNode = None
    self.remeshOversampling = ARG_DEFAULTS[ARG_REMESH_OVERSAMPLING]
    self.smoothingFactor = ARG_DEFAULTS[ARG_SMOOTHINGFACTOR]
    self.shrinkwrapIterations = ARG_DEFAULTS[ARG_SHRINKWRAP_ITERATIONS]
    self.saveIntermediateResults = ARG_DEFAULTS[ARG_SAVE_INTERMEDIATE_RESULTS]

    # Temporary variables
    self._inputPd = None
    self._inputPdCellLocator = None
    self._inputSpacing = None

  def requestCancel(self):
    logging.info("User requested cancelling.")
    self.cancelRequested = True

  def _log(self, message):
    if self.logCallback:
      self.logCallback(message)

  def _checkCancelRequested(self):
    if self.cancelRequested:
      self.checkCancelRequested = False
      raise ValueError("Cancel requested")

  def applyWrapSolidify(self):
    """Applies the Shrinkwrap-Raycast-Shrinkwrap Filter, a surface filter, to the selected passed segment.
    """
    self.cancelRequested = False

    self.intermediateResultCounter = 0
    self.previousIntermediateResult = None

    try:
      self._log('Initializing Filtering Process...')

      self._updateInputPd()
      regionPd = self._getInitialRegionPd()

      # TODO: implement this
      # if self.regionExtendDiameter:
      #   regionPd = self._extendRegion(regionPd)

      shrunkenPd = self._shrinkWrap(regionPd)

      self._log('Smoothing...')
      shrinkModelPD.DeepCopy(self._smoothPolydata(shrinkModelPD, options[ARG_SMOOTHINGFACTOR]))
      self._saveIntermediateResult("Smoothed", shrinkModelPD)

      # Create shell
      if options[ARG_CREATESHELL]:
        if options[ARG_PRESERVECRACKS]:
          self._checkCancelRequested()
          self._log('Shell - preserving cracks...')
          self._shellPreserveCracks(inputPD, shrinkModelPD, spacing)
          self._saveIntermediateResult("ShellRemovedCaps", shrinkModelPD)
        if options[ARG_SHELLTHICKNESS] > 0:
          self._checkCancelRequested()
          self._log('Shell - solidifying...')
          self._shellSolidify(shrinkModelPD, options[ARG_SHELLTHICKNESS])
          self._saveIntermediateResult("ShellSolidified", shrinkModelPD)

      # Write output to target node
      self._log('Save result...')
      name = segment.GetName()
      color = segment.GetColor()
      if options[ARG_OUTPUTTYPE] == OUTPUT_SEGMENTATION:
        if options[ARG_MODE] == MODE_INNER_SURFACE:
          # update seed from segment (that was be grown)
          updatedSegmentID = seedSegmentID
        else:
          updatedSegmentID = segmentID
        self._polydataToSegment(shrinkModelPD, segmentationNode, updatedSegmentID)
      elif options[ARG_OUTPUTTYPE] == OUTPUT_MODEL:
        self._polydataToModel(shrinkModelPD, name, color)
      else:
        raise ValueError('Unknown Output Type')

    finally:
      self._cleanup()


  def _cleanup(self):
    if self.previousIntermediateResult:
      self.previousIntermediateResult.GetDisplayNode().SetVisibility(False)
    self._log('')

  def _updateInputPd(self):

    segment = self.segmentationNode.GetSegmentation().GetSegment(self.segmentId)
    self._inputPd = vtk.vtkPolyData()

    # Get input polydata (inputPD) and input spacing
    if self.segmentationNode.GetSegmentation().GetMasterRepresentationName() == slicer.vtkSegmentationConverter().GetSegmentationBinaryLabelmapRepresentationName():
      # Master representation is binary labelmap
      # Reconvert to closed surface using chosen chosen smoothing factor
      originalSurfaceSmoothing = float(self.segmentationNode.GetSegmentation().GetConversionParameter(
        slicer.vtkBinaryLabelmapToClosedSurfaceConversionRule().GetSmoothingFactorParameterName()))
      if abs(originalSurfaceSmoothing-options[ARG_SMOOTHINGFACTOR]) > 0.001:
        self.segmentationNode.GetSegmentation().SetConversionParameter(
          slicer.vtkBinaryLabelmapToClosedSurfaceConversionRule().GetSmoothingFactorParameterName(), str(options[ARG_SMOOTHINGFACTOR]))
        self.segmentationNode.RemoveClosedSurfaceRepresentation()
        self.segmentationNode.CreateClosedSurfaceRepresentation()
      self.segmentationNode.GetClosedSurfaceRepresentation(segmentID, inputPD)
      # Get input spacing
      inputLabelmap = slicer.vtkOrientedImageData()
      self.segmentationNode.GetBinaryLabelmapRepresentation(segmentID, inputLabelmap)
      self._inputSpacing = math.sqrt(np.sum(np.array(inputLabelmap.GetSpacing())**2))
    else:
      # Representation is already closed surface
      self.segmentationNode.CreateClosedSurfaceRepresentation()
      self.segmentationNode.GetClosedSurfaceRepresentation(segmentID, inputPD)
      # set spacing to have an approxmately 250^3 volume
      # this size is not too large for average computing hardware yet
      # it is sufficiently detailed for many applications
      preferredVolumeSizeInVoxels = 250 * 250 * 250
      bounds = np.zeros(6)
      inputPD.GetBounds(bounds)
      volumeSizeInMm3 = (bounds[1] - bounds[0]) * (bounds[3] - bounds[2]) * (bounds[5] - bounds[4])
      self._inputSpacing = pow(volumeSizeInMm3 / preferredVolumeSizeInVoxels, 1 / 3.)

    self._inputPdCellLocator = vtk.vtkCellLocator()
    self._inputPdCellLocator.SetDataSet(inputPD)
    self._inputPdCellLocator.BuildLocator()


  def _getInitialRegionPd(self):

    # Get seed polydata (seedPD)
    if self.region == REGION_LARGEST:
      # TODO: implement largest cavity find
      # create sphere as seed (that will be shrunken)
      bounds = np.zeros(6)
      inputPD.GetBounds(bounds)
      diameters = np.array([bounds[1]-bounds[0],bounds[3]-bounds[2],bounds[5]-bounds[4]])
      maxRadius = max(diameters)/2.0
      sphereSource = vtk.vtkSphereSource()
      # to make sure the volume is fully included in the sphere, radius must be sqrt(2) times larger
      sphereSource.SetRadius(maxRadius*1.5)
      sphereSource.SetPhiResolution(max(int(maxRadius/5), 10))
      sphereSource.SetThetaResolution(max(int(maxRadius/5), 10))
      sphereSource.SetCenter((bounds[0]+bounds[1])/2.0, (bounds[2]+bounds[3])/2.0, (bounds[4]+bounds[5])/2.0)
      sphereSource.Update()
      initialRegionPd = sphereSource.GetOutput()
    elif self.region == REGION_SEGMENT:
      # create seed from segment (that will be grown)
      if not self.regionSegmentId:
        raise ValueError("Seed segment is not set")
      if self.regionSegmentId == self.segmentId:
        raise ValueError("Seed segment cannot be the same as the solidified segment")
      initialRegionPd = segmentationNode.GetClosedSurfaceInternalRepresentation(seedSegmentID)
      if not initialRegionPd or initialRegionPd.GetNumberOfPoints() == 0:
        raise ValueError("Region segment is empty")
      #seedPD = self._remeshPolydata(seedPD, inputSpacing*5.0)  # simplify the mesh

    cleanPolyData = vtk.vtkCleanPolyData()
    cleanPolyData.SetInputData(initialRegionPd)
    cleanPolyData.Update()
    initialRegionPd = cleanPolyData.GetOutput()

    self._saveIntermediateResult("InitialRegion", initialRegionPd)
    return initialRegionPd


  def _shrinkWrap(self, regionPd)

    spacing = self.inputSpacing / self.oversampling

    # Main shrinkwrap loop
    for iterationIndex in range(self.shrinkwrapIterations):
      # shrink
      self._checkCancelRequested()
      self._log('Shrinking %s/%s...' %(iterationIndex+1, shrinkwrapIterations))
      shrinkModelPD.DeepCopy(self._shrinkPolydata(shrinkModelPD, inputPD, inputPDCellLocator))
      self._saveIntermediateResult("Shrunken", shrinkModelPD)
      # remesh
      self._checkCancelRequested()
      self._log('Remeshing %s/%s...' %(iterationIndex+1, shrinkwrapIterations))
      remeshedPD = self._remeshPolydata(shrinkModelPD, spacing)
      shrinkModelPD = vtk.vtkPolyData()
      shrinkModelPD.DeepCopy(remeshedPD)

      self._saveIntermediateResult("Remeshed", shrinkModelPD)
      shrunkenPd


  def _polydataToModel(self, polydata, name, color):
    if not self.outputModel:
      self.outputModel = slicer.modules.models.logic().AddModel(polydata)
      self.outputModel.SetName(name)
      self.outputModel.GetDisplayNode().SliceIntersectionVisibilityOn()
    else:
      self.outputModel.SetAndObservePolyData(polydata)
      
    self.outputModel.GetDisplayNode().SetColor(color)

    return True

  def _polydataToSegment(self, polydata, segmentationNode, segmentID):
    # Get existing representations
    segmentation = segmentationNode.GetSegmentation()
    masterRepresentationName = segmentationNode.GetSegmentation().GetMasterRepresentationName()
    representationNames = []
    segmentation.GetContainedRepresentationNames(representationNames)
    # Update
    slicer.vtkSlicerSegmentationsModuleLogic.ClearSegment(segmentationNode, segmentID)
    wasModified = segmentationNode.StartModify()
    segment = segmentation.GetSegment(segmentID)
    segment.RemoveAllRepresentations()
    segment.AddRepresentation(vtkSegmentationCorePython.vtkSegmentationConverter.GetSegmentationClosedSurfaceRepresentationName(), polydata)
    segmentation.CreateRepresentation(masterRepresentationName)
    for representationName in representationNames:
      if representationName:
        # already converted
        continue
      segmentation.CreateRepresentation(representationName)
    segmentationNode.EndModify(wasModified)

  def _saveIntermediateResult(self, name, polydata, color=None):
    if not self.saveIntermediateResults:
      return

    # Show the last intermediate result only, hide previous
    if self.previousIntermediateResult:
      self.previousIntermediateResult.GetDisplayNode().SetVisibility(False)

    polyDataCopy = vtk.vtkPolyData()
    polyDataCopy.DeepCopy(polydata)
    outputModel = slicer.modules.models.logic().AddModel(polyDataCopy)
    outputModel.SetName("WrapSolidify-{0}-{1}".format(self.intermediateResultCounter, name))
    self.intermediateResultCounter += 1
    outputModel.GetDisplayNode().SliceIntersectionVisibilityOn()
    outputModel.GetDisplayNode().SetEdgeVisibility(True)
    outputModel.GetDisplayNode().SetBackfaceCulling(False)
    if color:
      outputModel.GetDisplayNode().SetColor(color)
    else:
      outputModel.GetDisplayNode().SetColor(1.0,1.0,0)
    self.previousIntermediateResult = outputModel

  def _shrinkPolydata(self, polydata, inputPolydata, inputPDCellLocator, distance=0):
    if distance == 0:
      smoothFilter = vtk.vtkSmoothPolyDataFilter()
      smoothFilter.SetInputData(0, polydata)
      smoothFilter.SetInputData(1, inputPolydata)
      smoothFilter.Update()
      return smoothFilter.GetOutput()

    points = polydata.GetPoints()

    for i in range(points.GetNumberOfPoints()):
      originPoint = np.array(points.GetPoint(i))
      closestPoint = np.array([0.0,0.0,0.0])
      cell = vtk.vtkGenericCell()
      cellId = vtk.mutable(0)
      subId = vtk.mutable(0)
      closestPointDist2 = vtk.mutable(0)

      inputPDCellLocator.FindClosestPoint(originPoint, closestPoint, cell, cellId, subId, closestPointDist2)

      vector = closestPoint - originPoint
      vectorLength = np.linalg.norm(vector)

      if distance > 0 and vectorLength > 0.01:
        if distance < vectorLength:  # do not go back inside the surface if the point was too close
          points.SetPoint(i, closestPoint - ((vector/vectorLength) * distance))
      else:
        points.SetPoint(i, closestPoint)
      
    points.Modified() #polydata.SetPoints(points)
    return polydata
  
  def _remeshPolydata(self, polydata, inputSpacing):
    bounds = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    polydata.GetBounds(bounds)

    whiteImage = vtk.vtkImageData()

    spacing = np.ones(3) * inputSpacing
    dim = [0, 0, 0]
    for i in range(3):
      # Add 3 to the dimensions to have at least 1 voxel thickness and 1 voxel boundary on both sides
      dim[i] = int(math.ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i])) + 3

    origin = [0.0, 0.0, 0.0]
    origin[0] = bounds[0] - spacing[0]
    origin[1] = bounds[2] - spacing[1]
    origin[2] = bounds[4] - spacing[2]

    whiteImage.SetOrigin(origin)
    whiteImage.SetSpacing(spacing)
    whiteImage.SetExtent(0, dim[0], 0, dim[1], 0, dim[2])
    whiteImage.AllocateScalars(vtk.VTK_UNSIGNED_CHAR,1)

    pol2stenc = vtk.vtkPolyDataToImageStencil()
    pol2stenc.SetInputData(polydata)
    pol2stenc.SetOutputOrigin(origin)
    pol2stenc.SetOutputSpacing(spacing)
    pol2stenc.SetOutputWholeExtent(whiteImage.GetExtent())
    pol2stenc.Update()

    imgstenc = vtk.vtkImageStencil()
    imgstenc.SetInputData(whiteImage)
    imgstenc.SetStencilConnection(pol2stenc.GetOutputPort())
    imgstenc.ReverseStencilOn()
    imgstenc.SetBackgroundValue(1)
    imgstenc.Update()

    revimgstenc = vtk.vtkImageStencil()
    revimgstenc.SetInputData(imgstenc.GetOutput())
    revimgstenc.SetStencilConnection(pol2stenc.GetOutputPort())
    revimgstenc.ReverseStencilOff()
    revimgstenc.SetBackgroundValue(0)
    revimgstenc.Update()

    discreteCubes = vtk.vtkDiscreteMarchingCubes()
    discreteCubes.SetInputConnection(revimgstenc.GetOutputPort())
    discreteCubes.GenerateValues(1,0,0)
    discreteCubes.Update()

    reverse = vtk.vtkReverseSense()
    reverse.SetInputConnection(discreteCubes.GetOutputPort())
    reverse.ReverseCellsOn()
    reverse.ReverseNormalsOn()
    reverse.Update()

    return reverse.GetOutput()

  def _smoothPolydata(self, polydata, smoothingFactor):
    passBand = pow(10.0, -4.0 * smoothingFactor)
    smootherSinc = vtk.vtkWindowedSincPolyDataFilter()
    smootherSinc.SetInputData(polydata)
    smootherSinc.SetNumberOfIterations(20)
    smootherSinc.FeatureEdgeSmoothingOff()
    smootherSinc.BoundarySmoothingOff()
    smootherSinc.NonManifoldSmoothingOn()
    smootherSinc.NormalizeCoordinatesOn()
    smootherSinc.Update()
    return smootherSinc.GetOutput()

  def _carveCavities(self, inputPD, inputPDCellLocator, seedPD, inputSpacing, shrink, cavitiesDiameter, cavitiesDepth):
    # This method expands the mesh in surface normal direction at regions where cavities are larger than the specified
    # diameter. To determine regions where large cavities are, we compute a zero-offset and a cavitiesDiameter-offset
    # surface. CavitiesDiameter-offset surface is "sunken" where the surface has holes larger than CavitiesDiameter,
    # therefore the places where large holes are those where these surfaces are closer to each other than CavitiesDiameter.

    # Few iterations of low-resolution shrinkwrap
    shrinkModelPD = vtk.vtkPolyData()
    shrinkModelPD.DeepCopy(seedPD)

    for mainIteration in range(15):

      numberOfIterations = 3
      # We need to be able to distinguish cavities that large because they are snapped to the edge of the cavity
      # therefore we need to make the resampling spacing smaller than the minimum cavity diameter.
      spacing = max(cavitiesDiameter * 0.5, inputSpacing)
      # need to get close to the surface so that the points snap to the edge of the cavity
      # and so the cavity will appear as a large cell
      offset = cavitiesDiameter
      for iterationIndex in range(numberOfIterations):
        # shrink
        self._checkCancelRequested()
        self._log('Carve cavities shrinking %s/%s...' %(iterationIndex+1, numberOfIterations))
        if iterationIndex < 2:
          shrunkdPD = self._shrinkPolydata(shrinkModelPD, inputPD, inputPDCellLocator, offset)
        else:
          shrunkdPD = self._shrinkPolydata(shrinkModelPD, inputPD, inputPDCellLocator, 0)
        if shrunkdPD.GetNumberOfPoints() == 0:
          # the polydata dissolved (probably spacing was too large), there are no small holes we could carve into
          shrinkModelPD.DeepCopy(seedPD)
          return shrinkModelPD
        shrinkModelPD = vtk.vtkPolyData()
        shrinkModelPD.DeepCopy(shrunkdPD)
        self._saveIntermediateResult("CarveCavitiesShrunken", shrinkModelPD)
        if iterationIndex == 1:
          shrinkModelPDWithDistance = vtk.vtkPolyData()
          shrinkModelPDWithDistance.DeepCopy(shrinkModelPD)
        elif iterationIndex == 2:
          # no need to remesh, we use the carved output
          break
        # remesh
        self._checkCancelRequested()
        self._log('Carve cavities remeshing %s/%s...' %(iterationIndex+1, numberOfIterations))
        remeshedPD = self._remeshPolydata(shrinkModelPD, spacing)
        shrinkModelPD = vtk.vtkPolyData()
        if remeshedPD.GetNumberOfPoints() == 0:
          # the polydata dissolved (probably spacing was too large), there are no small holes we could carve into
          shrinkModelPD.DeepCopy(seedPD)
          return shrinkModelPD
        shrinkModelPD.DeepCopy(remeshedPD)
        self._saveIntermediateResult("CarveCavitiesRemeshed", shrinkModelPD)

      # Raycast

      self._checkCancelRequested()
      self._log('Carve cavities raycasting...')

      shrinkModelPDCellLocator = vtk.vtkCellLocator()
      shrinkModelPDCellLocator.SetDataSet(shrinkModelPD)
      shrinkModelPDCellLocator.BuildLocator()

      points = shrinkModelPDWithDistance.GetPoints()
      numberOfPoints = points.GetNumberOfPoints()
      # Distance between convex hull and offset surface is smaller than
      # cavitiesDiameter at holes that have larger diameter than cavitiesDiameter,
      # (elsewhere it is approximately same as cavitiesDiameter).
      # We set the threshold slightly (10%) lower than cavitiesDiameter to prevent
      # false detection due to inaccuracies in distance estimation.
      minDistance2 = pow(offset, 2) * 0.9
      cavityPointIds = []
      distances=[]
      for pointId in range(numberOfPoints):
        convexHullPoint = points.GetPoint(pointId)
        closestPoint = [0.0, 0.0, 0.0]
        cell = vtk.vtkGenericCell()
        cellId = vtk.mutable(0)
        subId = vtk.mutable(0)
        closestPointDist2 = vtk.mutable(0)
        shrinkModelPDCellLocator.FindClosestPoint(convexHullPoint, closestPoint, cell, cellId, subId, closestPointDist2)
        distances.append(closestPointDist2)
        #if closestPointDist2 < minDistance2:
        cavityPointIds.append(pointId)

      # Smooth to make the normals less noisy
      smoothFilter = vtk.vtkSmoothPolyDataFilter()
      smoothFilter.SetInputData(shrinkModelPDWithDistance)
      smoothFilter.SetRelaxationFactor(0.1)
      smoothFilter.Update()

      # Generate normals
      normals = vtk.vtkPolyDataNormals()
      normals.ComputePointNormalsOn()
      normals.ComputeCellNormalsOff()
      normals.SplittingOff()
      normals.SetInputConnection(smoothFilter.GetOutputPort())
      normals.AutoOrientNormalsOn()
      normals.SetFlipNormals(shrink)
      normals.Update()

      if shrink:
        # smoothing reduces surface mesh size - expand it by 20% to compensate for that
        # (to prevent clipping of extruding sharp edges)
        expand = vtk.vtkWarpVector()
        expand.SetInputConnection(normals.GetOutputPort())
        expand.SetInputArrayToProcess(0, 0, 0,
          vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, vtk.vtkDataSetAttributes.NORMALS)
        bounds = np.zeros(6)
        shrinkModelPDWithDistance.GetBounds(bounds)
        maxDiameter = max([bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]])
        maxDiameter = 0
        expand.SetScaleFactor(-maxDiameter * 0.1)  # 20% extrusion
        expand.Update()
        shrinkModelPD.DeepCopy(expand.GetOutput())
        shrinkModelPD.GetPointData().AddArray(normals.GetOutput().GetPointData().GetArray("Normals"))
      else:
        shrinkModelPD.DeepCopy(normals.GetOutput())
      self._saveIntermediateResult("CarveCavitiesSmoothed", shrinkModelPD)

      # projection
      #numberOfCavityPoints = cavityPointIds.GetNumberOfIds()
      normalArray = shrinkModelPD.GetPointData().GetArray('Normals')
      points = shrinkModelPD.GetPoints()
      #for cavityPointIndex in range(numberOfCavityPoints):
      #  pointId = cavityPointIds.GetId(cavityPointIndex)
      minParametricDistance = offset / (cavitiesDepth + offset)
      for pointId in cavityPointIds:
        normal = np.array(normalArray.GetTuple(pointId))
        vector = normal * (cavitiesDepth + offset) # max Length greater 0, checked above
        # find intersection (point + vector) and label model
        originalPointPos = points.GetPoint(pointId)
        a1 = originalPointPos + vector
        tol = 1.0

        t = vtk.mutable(0)
        intersectionPos = np.array([0.0,0.0,0.0]) #global
        par = np.array([0.0,0.0,0.0]) #parametric
        cell = vtk.vtkGenericCell()
        cellId = vtk.mutable(0)
        subId = vtk.mutable(0)
        if inputPDCellLocator.IntersectWithLine(originalPointPos, a1, tol, t, intersectionPos, par, subId, cellId, cell):
          if t>minParametricDistance:
            if t > minParametricDistance * 2:
              t = minParametricDistance * 2
            updatedPos = originalPointPos + (t-minParametricDistance) * vector
            points.SetPoint(pointId, updatedPos)

      points.Modified()
      self._saveIntermediateResult("CarveCavitiesResult", shrinkModelPD)


    return shrinkModelPD

  def _shellPreserveCracks(self, inputPD, shrinkModelPD, spacing):
    """Remove cells of the mesh that are far from the original surface"""

    # implicit distance, add point ids with larger distance to ids
    implicitDistance = vtk.vtkImplicitPolyDataDistance()
    implicitDistance.SetInput(inputPD)
    
    # delete cells in great distance
    nonsolidPolyData = vtk.vtkPolyData()
    nonsolidPolyData.DeepCopy(shrinkModelPD)
    nonsolidPolyData.BuildLinks()

    shellDistance = 0.9  # a bit more than half diameter of a cube (0.5*sqrt(3))
    maxDistance = shellDistance * spacing

    for c in range(nonsolidPolyData.GetNumberOfCells()):
      cell = nonsolidPolyData.GetCell(c)
      points = cell.GetPoints()
      for p in range(points.GetNumberOfPoints()):
        point = points.GetPoint(p)
        distance = implicitDistance.EvaluateFunction(point)

        if abs(distance) > maxDistance:
          nonsolidPolyData.DeleteCell(c)
          break

    nonsolidPolyData.RemoveDeletedCells()
    shrinkModelPD.DeepCopy(nonsolidPolyData)


  def _shellSolidify(self, shrinkModelPD, shellThickness):
    """Create a thick shell from a surface by extruding in surface normal direction"""

    # remove double vertices
    cleanPolyData = vtk.vtkCleanPolyData()
    cleanPolyData.SetInputData(shrinkModelPD)
    cleanPolyData.Update()

    # create normals
    normals = vtk.vtkPolyDataNormals()
    normals.SetComputeCellNormals(1)
    normals.SetInputData(cleanPolyData.GetOutput())
    normals.SplittingOff()
    normals.Update()

    shrinkModelPD.DeepCopy(normals.GetOutput())
    numberOfPoints = shrinkModelPD.GetNumberOfPoints()

    # get boundary edges, used later
    featureEdges = vtk.vtkFeatureEdges()
    featureEdges.BoundaryEdgesOn()
    featureEdges.ColoringOff()
    featureEdges.FeatureEdgesOff()
    featureEdges.NonManifoldEdgesOff()
    featureEdges.ManifoldEdgesOff()
    featureEdges.SetInputData(normals.GetOutput())
    featureEdges.Update()

    addingPoints = []
    addingPolys = []

    for pointID in range(numberOfPoints):
      cellIDs = vtk.vtkIdList()
      shrinkModelPD.GetPointCells(pointID, cellIDs)
      normalsArray = []

      # ilterate through all cells/faces which contain point
      for i in range(cellIDs.GetNumberOfIds()):
        n = []
        n.append(shrinkModelPD.GetCellData().GetArray('Normals').GetValue(cellIDs.GetId(i)*3))
        n.append(shrinkModelPD.GetCellData().GetArray('Normals').GetValue(cellIDs.GetId(i)*3 + 1))
        n.append(shrinkModelPD.GetCellData().GetArray('Normals').GetValue(cellIDs.GetId(i)*3 + 2))

        normalsArray.append(np.array(n) * (-1))

      # calculate position of new vert
      dir_vec = np.zeros(3)

      for n in normalsArray:
        dir_vec = dir_vec + np.array(n)

      dir_vec_norm = dir_vec / np.linalg.norm(dir_vec)
      proj_length = np.dot(dir_vec_norm, np.array(normalsArray[0]))
      dir_vec_finallenght = dir_vec_norm * proj_length
      vertex_neu = np.array(shrinkModelPD.GetPoint(pointID)) + (dir_vec_finallenght * shellThickness)

      # append point
      addingPoints.append(vertex_neu)

    for cellID in range(shrinkModelPD.GetNumberOfCells()):
      pointIDs = vtk.vtkIdList()
      shrinkModelPD.GetCellPoints(cellID, pointIDs)

      newPointIDs = vtk.vtkIdList()
      for i in reversed(range(pointIDs.GetNumberOfIds())):
        newPointIDs.InsertNextId(int(pointIDs.GetId(i) + numberOfPoints))

      addingPolys.append(newPointIDs)

    doubleSurfacePoints = vtk.vtkPoints()
    doubleSurfacePolys = vtk.vtkCellArray()

    doubleSurfacePoints.DeepCopy(shrinkModelPD.GetPoints())
    doubleSurfacePolys.DeepCopy(shrinkModelPD.GetPolys())

    for p in addingPoints:
      doubleSurfacePoints.InsertNextPoint(p)
    for p in addingPolys:
      doubleSurfacePolys.InsertNextCell(p)

    doubleSurfacePD = vtk.vtkPolyData()
    doubleSurfacePD.SetPoints(doubleSurfacePoints)
    doubleSurfacePD.SetPolys(doubleSurfacePolys)

    # add faces to boundary edges
    mergePoints = vtk.vtkMergePoints()
    mergePoints.InitPointInsertion(doubleSurfacePD.GetPoints(), doubleSurfacePD.GetBounds())
    mergePoints.SetDataSet(doubleSurfacePD)
    mergePoints.BuildLocator()

    manifoldPolys = vtk.vtkCellArray()
    manifoldPolys.DeepCopy(doubleSurfacePD.GetPolys())
    manifoldPoints = vtk.vtkPoints()
    manifoldPoints.DeepCopy(doubleSurfacePD.GetPoints())

    for e in range(featureEdges.GetOutput().GetNumberOfCells()):
      pointIDs = vtk.vtkIdList()
      featureEdges.GetOutput().GetCellPoints(e, pointIDs)
      if pointIDs.GetNumberOfIds() == 2: # -> Edge
        matchingPointIDs = []
        newPointIDs = vtk.vtkIdList()
        for p in range(2):
          matchingPointIDs.append(mergePoints.IsInsertedPoint(featureEdges.GetOutput().GetPoint(pointIDs.GetId(p))))
        if not (-1) in matchingPointIDs: # edge vertex not found in original pd, should not happen
          newPointIDs.InsertNextId(matchingPointIDs[1])
          newPointIDs.InsertNextId(matchingPointIDs[0])
          newPointIDs.InsertNextId(matchingPointIDs[0]+numberOfPoints)
          newPointIDs.InsertNextId(matchingPointIDs[1]+numberOfPoints)
          manifoldPolys.InsertNextCell(newPointIDs)

    manifoldPD = vtk.vtkPolyData()
    manifoldPD.SetPoints(manifoldPoints)
    manifoldPD.SetPolys(manifoldPolys)

    triangleFilter = vtk.vtkTriangleFilter()
    triangleFilter.SetInputData(manifoldPD)
    triangleFilter.Update()

    shrinkModelPD.DeepCopy(triangleFilter.GetOutput())

ARG_DEFAULTS = {}

ARG_MODE = 'mode'
MODE_OUTER_SURFACE = 'outerSurface'
MODE_INNER_SURFACE = 'innerSurface'
ARG_OPTIONS[ARG_MODE] = [MODE_OUTER_SURFACE, MODE_INNER_SURFACE]
ARG_DEFAULTS[ARG_MODE] = MODE_OUTER_SURFACE

ARG_REGION = 'region'
REGION_LARGEST = 'largest'
REGION_SEGMENT = 'segment'
ARG_OPTIONS[ARG_REGION] = [REGION_LARGEST, REGION_SEGMENT]
ARG_DEFAULTS[ARG_REGION] = REGION_LARGEST

ARG_REGION_SEGMENT_ID = 'regionSegmentID'
ARG_DEFAULTS[ARG_REGION_SEGMENT_ID] = ''

ARG_REGION_EXTEND = 'extendRegion'
ARG_DEFAULTS[ARG_REGION_EXTEND] = False

ARG_REGION_EXTEND_DIAMETER = 'regionExtendDiameter'
ARG_DEFAULTS[ARG_REGION_EXTEND_DIAMETER] = 10.0

ARG_CREATE_SHELL = 'createShell'
ARG_DEFAULTS[ARG_CREATE_SHELL] = False

ARG_SHELL_THICKNESS = 'shellThickness'
ARG_DEFAULTS[ARG_SHELL_THICKNESS] = 1.5

ARG_SHELL_PRESERVE_CRACKS = 'preserveCracks'
ARG_DEFAULTS[ARG_SHELL_PRESERVE_CRACKS] = True

ARG_OUTPUT_TYPE = 'outputType'
OUTPUT_MODEL = 'model'
OUTPUT_SEGMENTATION = 'segmentation'
ARG_OPTIONS[ARG_OUTPUT_TYPE] = [OUTPUT_MODEL, OUTPUT_SEGMENTATION]
ARG_DEFAULTS[ARG_OUTPUT_TYPE] = OUTPUT_SEGMENTATION

ARG_OUTPUT_MODEL_NODE = 'WrapSolidify.OutputModelNodeID'

ARG_REMESH_OVERSAMPLING = 'remeshOversampling'
ARG_DEFAULTS[ARG_REMESH_OVERSAMPLING] = 1.5  # 1.5x oversampling by default

ARG_SMOOTHINGFACTOR = 'smoothingFactor'
ARG_DEFAULTS[ARG_SMOOTHINGFACTOR] = 0.2

ARG_ITERATIONS = 'iterations'
ARG_DEFAULTS[ARG_ITERATIONS] = 6

ARG_SAVE_INTERMEDIATE_RESULTS = 'saveIntermediateResults'
ARG_DEFAULTS[ARG_SAVE_INTERMEDIATE_RESULTS] = False
