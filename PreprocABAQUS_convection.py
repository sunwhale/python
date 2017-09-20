# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 15:09:26 2014

@author: Sun,Jinyu
"""

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
from Constants import *
from PreprocABAQUSParameters import *

executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(referenceRepresentation=ON)
Mdb()
#: A new model database has been created.
#: The model "Model-1" has been created.
#==============================================================================
# create part
#==============================================================================
session.viewports['Viewport: 1'].setValues(displayedObject=None)
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.sketchOptions.setValues(viewStyle=AXISYM)
s.setPrimaryObject(option=STANDALONE)
s.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
s.FixedConstraint(entity=g[2])
s.rectangle(point1=(5.00, 0.0), point2=(6.00, 1.0))
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=AXISYMMETRIC, type=DEFORMABLE_BODY, twist=ON)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']









#==============================================================================
# create material and section
#==============================================================================
mdb.models['Model-1'].Material(name='Material-1')
mdb.models['Model-1'].materials['Material-1'].Depvar(n=200)
mdb.models['Model-1'].materials['Material-1'].UserMaterial(mechanicalConstants=(0.0, ))
mdb.models['Model-1'].materials['Material-1'].Conductivity(table=((10.0, ), ))
mdb.models['Model-1'].materials['Material-1'].Expansion(table=((0.0e-06, ), ))
#mdb.models['Model-1'].materials['Material-1'].Expansion(table=((7.5e-06, ), ))
#mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1', material='Material-1', thickness=None)

materialname='GH4169'
mdb.models['Model-1'].Material(name=materialname)
mdb.models['Model-1'].materials[materialname].Density(temperatureDependency=ON, table=((8.24E-9, 20), ))
mdb.models['Model-1'].materials[materialname].Conductivity(temperatureDependency=ON, table=((13.4, 11), (14.7, 100), (15.9, 200), (17.8, 300), (18.3, 400), (19.6, 500), (21.2, 600), (22.8, 700), (23.6, 800), (27.6, 900), (30.4, 1000), ))
mdb.models['Model-1'].materials[materialname].SpecificHeat(temperatureDependency=ON, table=((481.4E6, 300), (493.9E6, 400), (514.8E6, 500), (539E6, 600), (573.4E6, 700), (615.3E6, 800), (657.2E6, 900), (707.4E6, 1000), ))
mdb.models['Model-1'].materials[materialname].Expansion(temperatureDependency=ON, table=((9.1E-6, -253), (11E-6, -183), (11.8E-6, 20), (11.8E-6, 100), (13E-6, 200), (13.5E-6, 300), (14.1E-6, 400), (14.4E-6, 500), (14.8E-6, 600), (15.4E-6, 700), (17E-6, 800), (18.4E-6, 900), (18.7E-6, 1000), ))
mdb.models['Model-1'].materials[materialname].Elastic(temperatureDependency=ON, table=((204E3, 0.3, 20), (181E3, 0.3, 300), (176E3, 0.31, 400), (160E3, 0.32, 500), (160E3, 0.32, 550), (150E3, 0.32, 600), (146E3, 0.325, 650), (141E3, 0.33, 700), ))
mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1', material=materialname, thickness=None)






#==============================================================================
# assign section
#==============================================================================
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
region = p.Set(faces=faces, name='Set-1')
p = mdb.models['Model-1'].parts['Part-1']
p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)









#==============================================================================
# assembly
#==============================================================================
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(coordSysType=CYLINDRICAL, origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), point2=(0.0, 0.0, -1.0))
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1-1', part=p, dependent=ON)









#==============================================================================
# mesh
#==============================================================================
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, engineeringFeatures=OFF, mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(meshTechnique=ON)
p1 = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
# seed edges 
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#1 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=5, constraint=FINER)

p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#8 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=5, constraint=FINER)
# element type
elemType1 = mesh.ElemType(elemCode=CGAX8T, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=CGAX6MT, elemLibrary=STANDARD, 
    hourglassControl=DEFAULT)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
pickedRegions =(faces, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
# generate mesh
p = mdb.models['Model-1'].parts['Part-1']
p.generateMesh()









#==============================================================================
# node set
#==============================================================================
p = mdb.models['Model-1'].parts['Part-1']
n = p.nodes
nodes = n.getSequenceFromMask(mask=('[#ffffffff:3 ]', ), )
p.Set(nodes=nodes, name='NALL')
p = mdb.models['Model-1'].parts['Part-1']
p.SetFromNodeLabels('N36',(36,))
p.SetFromNodeLabels('N33',(33,))
p.SetFromNodeLabels('N92',(92,))
#: The set 'N36' has been created (1 node).








#==============================================================================
# step
#==============================================================================
mdb.models['Model-1'].CoupledTempDisplacementStep(name='Step-1', 
    previous='Initial', response=STEADY_STATE, timePeriod=time_period,
    maxNumInc=1000000, initialInc=initial_inc, minInc=min_inc, maxInc=max_inc, deltmx=None, 
    cetol=None, creepIntegration=None, amplitude=RAMP
    , nlgeom=nonlinear
    )
#session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
    








#==============================================================================
# amplitude
#==============================================================================
# COS
mdb.models['Model-1'].PeriodicAmplitude(name='COS', timeSpan=STEP, frequency=1.0, start=0.0, a_0=0.0, data=((1.0, 0.0), ))
# SIN
mdb.models['Model-1'].PeriodicAmplitude(name='SIN', timeSpan=STEP, frequency=1.0, start=0.0, a_0=0.0, data=((0.0, 1.0), ))
# SIN Phase 0
mdb.models['Model-1'].PeriodicAmplitude(name='SIN0', timeSpan=STEP, frequency=1.0, start=0.0, a_0=0.0, data=((0.0, 1.0), ))
# SIN Phase 90
mdb.models['Model-1'].PeriodicAmplitude(name='SIN90', timeSpan=STEP, frequency=1.0, start=0.25, a_0=0.0, data=((0.0, 1.0), ))
# Displacement Triangular Wave
mdb.models['Model-1'].TabularAmplitude(name='DispTriangularWave', timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1.0, 1.0)))
# Pressure Triangular Wave
mdb.models['Model-1'].TabularAmplitude(name='PresTriangularWave', timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1.0, 1.0)))
# Rotation Triangular Wave
mdb.models['Model-1'].TabularAmplitude(name='RotaTriangularWave', timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1.0, 1.0)))
# Torque Triangular Wave
mdb.models['Model-1'].TabularAmplitude(name='TorqTriangularWave', timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1.0, 1.0)))
# Temperature Triangular Wave
mdb.models['Model-1'].TabularAmplitude(name='TempTriangularWave', timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1.0, 1.0)))









#==============================================================================
# predifined field
#==============================================================================
a = mdb.models['Model-1'].rootAssembly
region = a.instances['Part-1-1'].sets['Set-1']
mdb.models['Model-1'].Temperature(name='Predefined Field-1', createStepName='Initial', region=region, distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(predefined_temperature, ))









#==============================================================================
# boundary conditions
#==============================================================================
# bottom fix y direction
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-1-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#1 ]', ), )
region = regionToolset.Region(edges=edges1)
mdb.models['Model-1'].DisplacementBC(name='BC-1', createStepName='Step-1', region=region, u1=UNSET, u2=0.0, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
# bottom fix rotation y direction
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-1-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#1 ]', ), )
region = regionToolset.Region(edges=edges1)
mdb.models['Model-1'].DisplacementBC(name='BC-2', createStepName='Step-1', region=region, u1=UNSET, u2=UNSET, ur2=0.0, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
# top y direction tension
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-1-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#4 ]', ), )
region = regionToolset.Region(edges=edges1)
mdb.models['Model-1'].DisplacementBC(name='BC-3', createStepName='Step-1', region=region, u1=UNSET, u2=1.0, ur2=UNSET, ur3=UNSET, amplitude='DispTriangularWave', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
# top y direction rotation
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-1-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#4 ]', ), )
region = regionToolset.Region(edges=edges1)
mdb.models['Model-1'].DisplacementBC(name='BC-4', createStepName='Step-1', region=region, u1=UNSET, u2=UNSET, ur2=1.0, ur3=UNSET, amplitude='RotaTriangularWave', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
# temprature
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-1-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
region = regionToolset.Region(edges=edges1)
mdb.models['Model-1'].TemperatureBC(name='BC-5', createStepName='Step-1', region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', magnitude=1.0, amplitude='TempTriangularWave')
# convection
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['Part-1-1'].edges
side1Edges1 = s1.getSequenceFromMask(mask=('[#8 ]', ), )
region=regionToolset.Region(side1Edges=side1Edges1)
mdb.models['Model-1'].FilmCondition(name='inner_convection', 
    createStepName='Step-1', surface=region, definition=EMBEDDED_COEFF, 
    filmCoeff=0.3, filmCoeffAmplitude='', sinkTemperature=20.0, 
    sinkAmplitude='', sinkDistributionType=UNIFORM, sinkFieldName='')
    






#==============================================================================
# load
#==============================================================================
## pressure
#a = mdb.models['Model-1'].rootAssembly
#s1 = a.instances['Part-1-1'].edges
#side1Edges1 = s1.getSequenceFromMask(mask=('[#4 ]', ), )
#region = regionToolset.Region(side1Edges=side1Edges1)
#mdb.models['Model-1'].Pressure(name='Load-1', createStepName='Step-1', \
#    region=region, distributionType=UNIFORM, field='', magnitude=1.0, \
#    amplitude='PresTriangularWave')
## torque
## first create reference point
#p = mdb.models['Model-1'].parts['Part-1']
#p.ReferencePoint(point=(0.0, 1.0, 0.0))
## create constraint to couple reference point and top surface
#a = mdb.models['Model-1'].rootAssembly
#r1 = a.instances['Part-1-1'].referencePoints
##refPoints1=(r1[7], )
#refPoints1=(r1[r1.keys()[0]], )
#region1=regionToolset.Region(referencePoints=refPoints1)
#a = mdb.models['Model-1'].rootAssembly
#s1 = a.instances['Part-1-1'].edges
#side1Edges1 = s1.getSequenceFromMask(mask=('[#4 ]', ), )
#region2=regionToolset.Region(side1Edges=side1Edges1)
#mdb.models['Model-1'].Coupling(name='Constraint-1', controlPoint=region1, \
#    surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, \
#    localCsys=None, u1=OFF, u2=OFF, ur2=ON, ur3=OFF)
## add torque on reference point
#a = mdb.models['Model-1'].rootAssembly
#r1 = a.instances['Part-1-1'].referencePoints
##refPoints1=(r1[7], )
#refPoints1=(r1[r1.keys()[0]], )
#region = regionToolset.Region(referencePoints=refPoints1)
#mdb.models['Model-1'].Moment(name='Load-2', createStepName='Step-1', \
#    region=region, cm2=1.0, amplitude='TorqTriangularWave', \
#    distributionType=UNIFORM, field='', localCsys=None)









#==============================================================================
# field output
#==============================================================================
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'E', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 'NT', 'TEMP', 'HFL', 'RFL'))
#==============================================================================
# history output
#==============================================================================
mdb.models['Model-1'].historyOutputRequests['H-Output-1'].setValues(variables=PRESELECT)
#==============================================================================
# job
#==============================================================================
mdb.Job(name=JobName, model='Model-1', \
    description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, \
    queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, \
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, \
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', \
    scratch='', multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
#==============================================================================
# write input file
#==============================================================================
mdb.jobs[JobName].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "AxisymmetricTensionTorsionThermal.inp".
#==============================================================================
# save cae file
#==============================================================================
pathName = AbaqusWorkDirectory+CAEName
mdb.saveAs(pathName=pathName)