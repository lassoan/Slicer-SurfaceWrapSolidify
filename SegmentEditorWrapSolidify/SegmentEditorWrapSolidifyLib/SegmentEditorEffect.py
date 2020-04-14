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

    self.logic = WrapSolidifyLogic(scriptedEffect)
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

    # Set defaults
    if DEFAULT_MODE == MODE_OUTER_SURFACE:
      self.ui.modeOuterSurfaceRadioButton.setChecked(True)
    elif DEFAULT_MODE == MODE_INNER_SURFACE:
      self.ui.modeInnerSurfaceRadioButton.setChecked(True)
    self.ui.smoothingFactorSlider.value = DEFAULT_SMOOTHINGFACTOR
    self.ui.cavitiesCheckBox.setChecked(DEFAULT_CARVECAVITIES)
    self.ui.cavitiesDiameterSlider.value = DEFAULT_CAVITIESDIAMETER
    self.ui.cavitiesDepthSlider.value = DEFAULT_CAVITIESDEPTH
    self.ui.createShellCheckBox.setChecked(DEFAULT_CREATESHELL)
    self.ui.shellPreserveCracksCheckBox.setChecked(DEFAULT_PRESERVECRACKS)
    self.ui.saveIntermediateResultsCheckBox.setChecked(DEFAULT_SAVEINTERMEDIATERESULTS)

    self.ui.shellThicknessSlider.value = DEFAULT_SHELLTHICKNESS
    if DEFAULT_OUTPUTTYPE == OUTPUT_SEGMENTATION:
      self.ui.outputSegmentationRadioButton.setChecked(True)
    elif DEFAULT_OUTPUTTYPE == OUTPUT_MODEL:
      self.ui.outputModelRadioButton.setChecked(True)
    self.ui.iterationsSlider.value = DEFAULT_ITERATIONS
    self.ui.remeshOversamplingSlider.value = DEFAULT_SPACING

    # Add connections

    self.ui.modeGroup = qt.QButtonGroup()
    self.ui.modeGroup.addButton(self.ui.modeOuterSurfaceRadioButton)
    self.ui.modeGroup.addButton(self.ui.modeInnerSurfaceRadioButton)

    self.outputTypeGroup = qt.QButtonGroup()
    self.outputTypeGroup.addButton(self.ui.outputSegmentationRadioButton)
    self.outputTypeGroup.addButton(self.ui.outputModelRadioButton)

    self.ui.modeGroup.connect('buttonClicked(int)', self.updateMRMLFromGUI)
    self.ui.innerSurfaceSeedSegmentSelector.connect("currentSegmentChanged(QString)", self.updateMRMLFromGUI)
    self.ui.smoothingFactorSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.ui.cavitiesCheckBox.connect('stateChanged(int)', self.updateMRMLFromGUI)
    self.ui.cavitiesDiameterSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.ui.cavitiesDepthSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.ui.createShellCheckBox.connect('stateChanged(int)', self.updateMRMLFromGUI)
    self.ui.shellPreserveCracksCheckBox.connect('stateChanged(int)', self.updateMRMLFromGUI)
    self.ui.shellThicknessSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.outputTypeGroup.connect('buttonClicked(int)', self.updateMRMLFromGUI)
    self.ui.outputModelSelector.connect('currentNodeChanged(vtkMRMLNode*)', self.updateMRMLFromGUI)
    self.ui.iterationsSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)    
    self.ui.remeshOversamplingSlider.connect('valueChanged(double)', self.updateMRMLFromGUI)
    self.ui.applyButton.connect('clicked()', self.onApply)
    self.ui.saveIntermediateResultsCheckBox.connect('stateChanged(int)', self.updateMRMLFromGUI)

  def createCursor(self, widget):
    return slicer.util.mainWindow().cursor

  def layoutChanged(self):
    pass

  def processInteractionEvents(self, callerInteractor, eventId, viewWidget):
    return False # For the sake of example

  def processViewNodeEvents(self, callerViewNode, eventId, viewWidget):
    pass # For the sake of example

  def setMRMLDefaults(self):
    self.scriptedEffect.setParameterDefault(ARG_MODE, DEFAULT_MODE)
    self.scriptedEffect.setParameterDefault(ARG_INNER_SURFACE_SEED_SEGMENT_ID, None)

    self.scriptedEffect.setParameterDefault(ARG_CARVECAVITIES, DEFAULT_CARVECAVITIES)
    self.scriptedEffect.setParameterDefault(ARG_CAVITIESDIAMETER, DEFAULT_CAVITIESDIAMETER)
    self.scriptedEffect.setParameterDefault(ARG_CAVITIESDEPTH, DEFAULT_CAVITIESDEPTH)

    self.scriptedEffect.setParameterDefault(ARG_CREATESHELL, DEFAULT_CREATESHELL)
    self.scriptedEffect.setParameterDefault(ARG_SHELLTHICKNESS, DEFAULT_SHELLTHICKNESS)
    self.scriptedEffect.setParameterDefault(ARG_PRESERVECRACKS, DEFAULT_PRESERVECRACKS)

    self.scriptedEffect.setParameterDefault(ARG_OUTPUTTYPE, DEFAULT_OUTPUTTYPE)
    self.scriptedEffect.setParameterDefault(ARG_SMOOTHINGFACTOR, DEFAULT_SMOOTHINGFACTOR)
    self.scriptedEffect.parameterSetNode().SetNodeReferenceID(ARG_OUTPUTMODELNODE, None)
    
    self.scriptedEffect.setParameterDefault(ARG_ITERATIONS, DEFAULT_ITERATIONS)
    self.scriptedEffect.setParameterDefault(ARG_SPACING, DEFAULT_SPACING)
    self.scriptedEffect.setParameterDefault(ARG_SAVEINTERMEDIATERESULTS, DEFAULT_SAVEINTERMEDIATERESULTS)

  def updateGUIFromMRML(self):
    for pId, pElement in [
        (ARG_SMOOTHINGFACTOR, self.ui.smoothingFactorSlider),
        (ARG_CAVITIESDIAMETER, self.ui.cavitiesDiameterSlider),
        (ARG_CAVITIESDEPTH, self.ui.cavitiesDepthSlider),
        (ARG_SHELLTHICKNESS, self.ui.shellThicknessSlider),
        (ARG_ITERATIONS, self.ui.iterationsSlider),
        (ARG_SPACING, self.ui.remeshOversamplingSlider)
      ]:
      value = self.scriptedEffect.doubleParameter(pId)
      wasBlocked = pElement.blockSignals(True)
      pElement.value = value
      pElement.blockSignals(wasBlocked)

    wasBlocked1 = self.ui.modeOuterSurfaceRadioButton.blockSignals(True)
    wasBlocked2 = self.ui.modeInnerSurfaceRadioButton.blockSignals(True)
    modeID = self.scriptedEffect.parameter(ARG_MODE)
    self.ui.modeOuterSurfaceRadioButton.setChecked(modeID == MODE_OUTER_SURFACE)
    self.ui.modeInnerSurfaceRadioButton.setChecked(modeID == MODE_INNER_SURFACE)
    self.ui.modeOuterSurfaceRadioButton.blockSignals(wasBlocked1)
    self.ui.modeInnerSurfaceRadioButton.blockSignals(wasBlocked2)
    self.ui.innerSurfaceSeedSegmentSelector.enabled = (modeID == MODE_INNER_SURFACE)

    wasBlocked = self.ui.innerSurfaceSeedSegmentSelector.blockSignals(True)
    self.ui.innerSurfaceSeedSegmentSelector.setCurrentNode(self.scriptedEffect.parameterSetNode().GetSegmentationNode())
    innerSurfaceSeedSegmentId = self.scriptedEffect.parameter(ARG_INNER_SURFACE_SEED_SEGMENT_ID)
    if innerSurfaceSeedSegmentId:
      self.ui.innerSurfaceSeedSegmentSelector.setCurrentSegmentID(innerSurfaceSeedSegmentId)
    self.ui.innerSurfaceSeedSegmentSelector.blockSignals(wasBlocked)

    wasBlocked = self.ui.cavitiesCheckBox.blockSignals(True)
    carveCavities = self.scriptedEffect.parameter(ARG_CARVECAVITIES)=='True'
    self.ui.cavitiesCheckBox.setChecked(carveCavities)
    self.ui.cavitiesCheckBox.blockSignals(wasBlocked)
    self.ui.cavitiesDiameterSlider.enabled = carveCavities
    self.ui.cavitiesDepthSlider.enabled = carveCavities

    wasBlocked = self.ui.createShellCheckBox.blockSignals(True)
    createShell = self.scriptedEffect.parameter(ARG_CREATESHELL) == 'True'
    self.ui.createShellCheckBox.setChecked(createShell)
    self.ui.createShellCheckBox.blockSignals(wasBlocked)
    self.ui.shellThicknessSlider.enabled = createShell
    self.ui.shellPreserveCracksCheckBox.enabled = createShell

    wasBlocked = self.ui.shellPreserveCracksCheckBox.blockSignals(True)
    preserveCracks = self.scriptedEffect.parameter(ARG_PRESERVECRACKS) == 'True'
    self.ui.shellPreserveCracksCheckBox.setChecked(preserveCracks)
    self.ui.shellPreserveCracksCheckBox.blockSignals(wasBlocked)

    wasBlocked = self.ui.saveIntermediateResultsCheckBox.blockSignals(True)
    self.ui.saveIntermediateResultsCheckBox.setChecked(self.scriptedEffect.parameter(ARG_SAVEINTERMEDIATERESULTS)=='True')
    self.ui.saveIntermediateResultsCheckBox.blockSignals(wasBlocked)

    wasBlocked1 = self.ui.outputSegmentationRadioButton.blockSignals(True)
    wasBlocked2 = self.ui.outputModelRadioButton.blockSignals(True)
    outputTypeID = self.scriptedEffect.parameter(ARG_OUTPUTTYPE)
    self.ui.outputSegmentationRadioButton.setChecked(outputTypeID == OUTPUT_SEGMENTATION)
    self.ui.outputModelRadioButton.setChecked(outputTypeID != OUTPUT_SEGMENTATION)
    self.ui.outputSegmentationRadioButton.blockSignals(wasBlocked1)
    self.ui.outputModelRadioButton.blockSignals(wasBlocked2)

    wasBlocked = self.ui.outputModelSelector.blockSignals(True)
    self.ui.outputModelSelector.setCurrentNode(slicer.mrmlScene.GetNodeByID(self.scriptedEffect.parameter(ARG_OUTPUTMODELNODE)))
    self.ui.outputModelSelector.blockSignals(wasBlocked)

    segmentation = self.scriptedEffect.parameterSetNode().GetSegmentationNode().GetSegmentation()
    if modeID == MODE_OUTER_SURFACE:
      outputSegmentId = self.scriptedEffect.parameterSetNode().GetSelectedSegmentID()
      self.ui.outputSegmentNameLabel.text = segmentation.GetSegment(outputSegmentId).GetName()
    else:  # inner surface
      outputSegmentId = innerSurfaceSeedSegmentId
      if outputSegmentId == self.scriptedEffect.parameterSetNode().GetSelectedSegmentID():
        self.ui.outputSegmentNameLabel.text = "error - seed segment cannot be the same as the current segment"
      else:
        self.ui.outputSegmentNameLabel.text = segmentation.GetSegment(outputSegmentId).GetName()

    outputMod = slicer.mrmlScene.GetNodeByID(self.scriptedEffect.parameter(ARG_OUTPUTMODELNODE))


  def updateMRMLFromGUI(self):
    wasModified = self.scriptedEffect.parameterSetNode().StartModify()
    for pId, pElement in [
        (ARG_SMOOTHINGFACTOR, self.ui.smoothingFactorSlider),
        (ARG_CAVITIESDIAMETER, self.ui.cavitiesDiameterSlider),
        (ARG_CAVITIESDEPTH, self.ui.cavitiesDepthSlider),
        (ARG_SHELLTHICKNESS, self.ui.shellThicknessSlider),
        (ARG_ITERATIONS, self.ui.iterationsSlider),
        (ARG_SPACING, self.ui.remeshOversamplingSlider)
        ]:
      self.scriptedEffect.setParameter(pId, pElement.value)

    self.scriptedEffect.setParameter(ARG_MODE, MODE_OUTER_SURFACE if self.ui.modeOuterSurfaceRadioButton.isChecked() else MODE_INNER_SURFACE)
    innerSurfaceSeedSegmentId = self.ui.innerSurfaceSeedSegmentSelector.currentSegmentID()
    if innerSurfaceSeedSegmentId:
      self.scriptedEffect.setParameter(ARG_INNER_SURFACE_SEED_SEGMENT_ID, innerSurfaceSeedSegmentId)
    self.scriptedEffect.setParameter(ARG_CARVECAVITIES, self.ui.cavitiesCheckBox.isChecked())
    self.scriptedEffect.setParameter(ARG_CREATESHELL, self.ui.createShellCheckBox.isChecked())
    self.scriptedEffect.setParameter(ARG_PRESERVECRACKS, self.ui.shellPreserveCracksCheckBox.isChecked())
    self.scriptedEffect.setParameter(ARG_OUTPUTTYPE, OUTPUT_SEGMENTATION if self.ui.outputSegmentationRadioButton.isChecked() else OUTPUT_MODEL)
    self.scriptedEffect.parameterSetNode().SetNodeReferenceID(ARG_OUTPUTMODELNODE, self.ui.outputModelSelector.currentNodeID)
    self.scriptedEffect.setParameter(ARG_SAVEINTERMEDIATERESULTS, self.ui.saveIntermediateResultsCheckBox.isChecked())

    self.scriptedEffect.parameterSetNode().EndModify(wasModified)

  #
  # Effect specific methods (the above ones are the API methods to override)
  #

  def onApply(self):

    if self.ui.applyButton.text == 'Cancel':
      self.logic.requestCancel()
      return

    self.scriptedEffect.saveStateForUndo()

    self.ui.applyButton.text = 'Cancel'
    qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)

    seg = self.scriptedEffect.parameterSetNode().GetSegmentationNode()
    segID = self.scriptedEffect.parameterSetNode().GetSelectedSegmentID()
    modID = self.scriptedEffect.parameter(ARG_OUTPUTMODELNODE)
    outputMod = slicer.mrmlScene.GetNodeByID(self.scriptedEffect.parameter(ARG_OUTPUTMODELNODE))

    kwargs = {
      ARG_MODE : self.scriptedEffect.parameter(ARG_MODE),
      ARG_INNER_SURFACE_SEED_SEGMENT_ID : self.scriptedEffect.parameter(ARG_INNER_SURFACE_SEED_SEGMENT_ID),
      ARG_CARVECAVITIES : self.scriptedEffect.parameter(ARG_CARVECAVITIES)=='True',
      ARG_CAVITIESDIAMETER : self.scriptedEffect.doubleParameter(ARG_CAVITIESDIAMETER),
      ARG_CAVITIESDEPTH : self.scriptedEffect.doubleParameter(ARG_CAVITIESDEPTH),

      ARG_CREATESHELL : self.scriptedEffect.parameter(ARG_CREATESHELL)=='True',
      ARG_SHELLTHICKNESS : self.scriptedEffect.doubleParameter(ARG_SHELLTHICKNESS),
      ARG_PRESERVECRACKS: self.scriptedEffect.parameter(ARG_PRESERVECRACKS) == 'True',

      ARG_OUTPUTTYPE : self.scriptedEffect.parameter(ARG_OUTPUTTYPE),
      ARG_SMOOTHINGFACTOR : self.scriptedEffect.doubleParameter(ARG_SMOOTHINGFACTOR),

      ARG_ITERATIONS : self.scriptedEffect.doubleParameter(ARG_ITERATIONS),
      ARG_SPACING : self.scriptedEffect.doubleParameter(ARG_SPACING),
      ARG_SAVEINTERMEDIATERESULTS : self.scriptedEffect.parameter(ARG_SAVEINTERMEDIATERESULTS)=='True'
    }
    
    if not self.logic.applyWrapSolidify(seg, segID, outputMod, **kwargs):
      slicer.util.errorDisplay("Wrap solidfy failed: " + self.logic.resultMessage)

    qt.QApplication.restoreOverrideCursor()
    self.ui.applyButton.text = 'Apply'


  def addLog(self, text):
    slicer.util.showStatusMessage(text)
    slicer.app.processEvents() # force update



class WrapSolidifyLogic(object):

  def __init__(self, scriptedEffect):
    self.scriptedEffect = scriptedEffect
    self.logCallback = None
    self.cancelRequested = False
  
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

  def applyWrapSolidify(self, segmentationNode, segmentID, modelNode=None, **kwargs):
    """Applies the Shrinkwrap-Raycast-Shrinkwrap Filter, a surface filter, to the selected passed segment.
    
    Arguments:
        segmentationNode (vtkMRMLSegmentationNode): Segmentation Node the filter gets applied to.
        segmentID (string): ID of the Segment the filter gets applied to. WARNING: The filter replaces the former segmentation.
        modelNode (vtkMRMLModelNode): If 'MODEL' is selected as outputType, the polydata of this node gets updated.
        **mode (string): outerSurface, innerSurface
        **outputType (string): Possible options: 'MODEL', 'SEGMENTATION'
        **smoothingFactor (double): 0-1, Smoothing of the surface representation of the input used by this filter. This factor is also used on the output surface model and segmentation.
        **carveCavities (bool): Solidification process also carves out deeper caves.
        **cavitiesDiameter (double): 0.1-100, Entrance diameter of caves. Only used if carveCavities == True.
        **cavitiesDepth (double): 0.1-1000, Depth of caves. Only used if carveCavities == True.
        **createShell (bool): Whether or not the solidification will be done only for a outer shell.
        **preserveCracks (bool): Add gaps where the input surface has cracks.
        **shellThickness (-0.1-10): Thickness of the shell. Only used if createShell == True. WARNING: If <0, a nonmanifold mesh gets created, Segmentation will fail.
        **iterationsNr (1-20): Nr. of iterations of the shrinkwrap process.
        **spacing (>0.1): Spacing of remesh process. 1.0 means remesh at current resolution
    Returns:
        bool: return True on success.
    """

    self.cancelRequested = False
    self.segLogic = slicer.vtkSlicerSegmentationsModuleLogic
    self.modelsLogic = slicer.modules.models.logic()

    options = {
      ARG_MODE : DEFAULT_MODE,

      ARG_OUTPUTTYPE : DEFAULT_OUTPUTTYPE,
      ARG_SMOOTHINGFACTOR : DEFAULT_SMOOTHINGFACTOR,

      ARG_CARVECAVITIES : DEFAULT_CARVECAVITIES,
      ARG_CAVITIESDIAMETER : DEFAULT_CAVITIESDIAMETER,
      ARG_CAVITIESDEPTH : DEFAULT_CAVITIESDEPTH,

      ARG_CREATESHELL : DEFAULT_CREATESHELL,
      ARG_SHELLTHICKNESS : DEFAULT_SHELLTHICKNESS,
      ARG_PRESERVECRACKS: DEFAULT_PRESERVECRACKS,

      ARG_ITERATIONS : DEFAULT_ITERATIONS,
      ARG_SPACING : DEFAULT_SPACING,
      ARG_SAVEINTERMEDIATERESULTS : DEFAULT_SAVEINTERMEDIATERESULTS
      }

    options.update(kwargs)

    self.saveIntermediateResults = options[ARG_SAVEINTERMEDIATERESULTS]
    self.intermediateResultCounter = 0
    self.previousIntermediateResult = None
    success = True
    self.resultMessage = ""

    try:
      self._log('Initializing Filtering Process...')

      segment = segmentationNode.GetSegmentation().GetSegment(segmentID)
      inputPD = vtk.vtkPolyData()

      # Get input polydata (inputPD) and input spacing
      if segmentationNode.GetSegmentation().GetMasterRepresentationName() == slicer.vtkSegmentationConverter().GetSegmentationBinaryLabelmapRepresentationName():
        # Master representation is binary labelmap
        # Reconvert to closed surface using chosen chosen smoothing factor
        originalSurfaceSmoothing = float(segmentationNode.GetSegmentation().GetConversionParameter(slicer.vtkBinaryLabelmapToClosedSurfaceConversionRule().GetSmoothingFactorParameterName()))
        if abs(originalSurfaceSmoothing-options[ARG_SMOOTHINGFACTOR]) > 0.001:
          segmentationNode.GetSegmentation().SetConversionParameter(slicer.vtkBinaryLabelmapToClosedSurfaceConversionRule().GetSmoothingFactorParameterName(), str(options[ARG_SMOOTHINGFACTOR]))
          segmentationNode.RemoveClosedSurfaceRepresentation()
          segmentationNode.CreateClosedSurfaceRepresentation()
        segmentationNode.GetClosedSurfaceRepresentation(segmentID, inputPD)
        # Get input spacing
        inputLabelmap = slicer.vtkOrientedImageData()
        segmentationNode.GetBinaryLabelmapRepresentation(segmentID, inputLabelmap)
        inputSpacing = math.sqrt(np.sum(np.array(inputLabelmap.GetSpacing())**2))
      else:
        # Representation is already closed surface
        segmentationNode.CreateClosedSurfaceRepresentation()
        segmentationNode.GetClosedSurfaceRepresentation(segmentID, inputPD)
        # set spacing to have an approxmately 250^3 volume
        # this size is not too large for average computing hardware yet
        # it is sufficiently detailed for many applications
        preferredVolumeSizeInVoxels = 250 * 250 * 250
        bounds = np.zeros(6)
        inputPD.GetBounds(bounds)
        volumeSizeInMm3 = (bounds[1] - bounds[0]) * (bounds[3] - bounds[2]) * (bounds[5] - bounds[4])
        inputSpacing = pow(volumeSizeInMm3 / preferredVolumeSizeInVoxels, 1 / 3.)

      inputPDCellLocator = vtk.vtkCellLocator()
      inputPDCellLocator.SetDataSet(inputPD)
      inputPDCellLocator.BuildLocator()
      self.outputModel = modelNode

      # Get seed polydata (seedPD)
      if options[ARG_MODE] == MODE_OUTER_SURFACE:
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
        seedPD = sphereSource.GetOutput()
      elif options[ARG_MODE] == MODE_INNER_SURFACE:
        # create seed from segment (that will be grown)
        seedSegmentID = options[ARG_INNER_SURFACE_SEED_SEGMENT_ID]
        if not seedSegmentID:
          raise ValueError("Seed segment is not set")
        if seedSegmentID == segmentID:
          raise ValueError("Seed segment cannot be the same as the solidified segment")
        seedPD = segmentationNode.GetClosedSurfaceInternalRepresentation(seedSegmentID)
        if not seedPD or seedPD.GetNumberOfPoints() == 0:
          raise ValueError("Seed segment is empty")
        #seedPD = self._remeshPolydata(seedPD, inputSpacing*5.0)  # simplify the mesh

      self._saveIntermediateResult("Seed", seedPD)
      cleanPolyData = vtk.vtkCleanPolyData()
      cleanPolyData.SetInputData(seedPD)
      cleanPolyData.Update()

      shrinkModelPD = vtk.vtkPolyData()
      shrinkModelPD.DeepCopy(seedPD)

      # Carve cavities
      spacing = inputSpacing / options[ARG_SPACING]
      if options[ARG_CARVECAVITIES]:
        shrinkModelPD = self._carveCavities(inputPD, inputPDCellLocator, shrinkModelPD, spacing, options[ARG_MODE] == MODE_OUTER_SURFACE,
                            options[ARG_CAVITIESDIAMETER], options[ARG_CAVITIESDEPTH])

      # Main shrinkwrap loop
      numberOfIterations = int(options[ARG_ITERATIONS])
      for iterationIndex in range(numberOfIterations):
        # shrink
        self._checkCancelRequested()
        self._log('Shrinking %s/%s...' %(iterationIndex+1, numberOfIterations))
        shrinkModelPD.DeepCopy(self._shrinkPolydata(shrinkModelPD, inputPD, inputPDCellLocator))
        self._saveIntermediateResult("Shrunken", shrinkModelPD)
        # remesh
        self._checkCancelRequested()
        self._log('Remeshing %s/%s...' %(iterationIndex+1, numberOfIterations))
        remeshedPD = self._remeshPolydata(shrinkModelPD, spacing)
        shrinkModelPD = vtk.vtkPolyData()
        shrinkModelPD.DeepCopy(remeshedPD)

        self._saveIntermediateResult("Remeshed", shrinkModelPD)

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

    except Exception as e:
      import traceback
      traceback.print_exc()
      self.resultMessage = str(e)
      success = False
    finally:
      self._cleanup()

    return success

  def _cleanup(self):
    if self.previousIntermediateResult:
      self.previousIntermediateResult.GetDisplayNode().SetVisibility(False)
    self._log('')

  def _polydataToModel(self, polydata, name, color):
    if not self.outputModel:
      self.outputModel = self.modelsLogic.AddModel(polydata)
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
    outputModel = self.modelsLogic.AddModel(polyDataCopy)
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
      if closestPointDist2 < minDistance2:
        cavityPointIds.append(pointId)

    # Smooth to make the normals less noisy
    smoothFilter = vtk.vtkSmoothPolyDataFilter()
    smoothFilter.SetInputData(shrinkModelPDWithDistance)
    smoothFilter.SetRelaxationFactor(0.3)
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
        points.SetPoint(pointId, intersectionPos[0], intersectionPos[1], intersectionPos[2])
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



ARG_MODE = 'mode'
MODE_OUTER_SURFACE = 'outerSurface'
MODE_INNER_SURFACE = 'innerSurface'
DEFAULT_MODE = MODE_OUTER_SURFACE

ARG_INNER_SURFACE_SEED_SEGMENT_ID = 'innerSurfaceSegmentID'

ARG_SMOOTHINGFACTOR = 'smoothingFactor'
DEFAULT_SMOOTHINGFACTOR = 0.2

ARG_CARVECAVITIES = 'carveCavities'
DEFAULT_CARVECAVITIES = False
ARG_CAVITIESDIAMETER = 'cavitiesDiameter'
DEFAULT_CAVITIESDIAMETER = 20.0
ARG_CAVITIESDEPTH = 'cavitiesDepth'
DEFAULT_CAVITIESDEPTH = 100.0

ARG_CREATESHELL = 'createShell'
DEFAULT_CREATESHELL = False
ARG_SHELLTHICKNESS = 'shellThickness'
DEFAULT_SHELLTHICKNESS = 1.5
ARG_PRESERVECRACKS = 'preserveCracks'
DEFAULT_PRESERVECRACKS = True

ARG_OUTPUTTYPE = 'outputType'
OUTPUT_MODEL = 'MODEL'
OUTPUT_SEGMENTATION = 'SEGMENTATION'
DEFAULT_OUTPUTTYPE = OUTPUT_SEGMENTATION

ARG_OUTPUTMODELNODE = 'WrapSolidify.OutputModelNodeID'
DEFAULT_OUTPUTMODELNODE = ''

ARG_ITERATIONS = 'iterationsNr'
DEFAULT_ITERATIONS = 6
ARG_SPACING = 'spacing'
DEFAULT_SPACING = 1.5  # 1.5x oversampling by default
ARG_SAVEINTERMEDIATERESULTS = 'saveIntermediateResults'
DEFAULT_SAVEINTERMEDIATERESULTS = False
