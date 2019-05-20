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
import json
import numpy as np

def read_json_file(file_name):
    """
    Writing JSON data to file.
    """
    with open(file_name,'r') as data_file:
        return json.loads(data_file.read())

def write_json_file(file_name, data):
    """
    Reading JSON data to file.
    """
    with open(file_name,'w') as data_file:
        return json.dump(data, data_file)

executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(referenceRepresentation=ON)
Mdb()
#: A new model database has been created.
#: The model "Model-1" has been created.
#==============================================================================
# create part
#==============================================================================
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.Line(point1=(4.25, 15.0), point2=(4.25, 0.0))
s.VerticalConstraint(entity=g[2], addUndoState=False)
s.Line(point1=(4.25, 0.0), point2=(3.25, 0.0))
s.HorizontalConstraint(entity=g[3], addUndoState=False)
s.PerpendicularConstraint(entity1=g[2], entity2=g[3], addUndoState=False)
s.Line(point1=(3.25, 0.0), point2=(3.25, 60.0))
s.VerticalConstraint(entity=g[4], addUndoState=False)
s.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
s.Line(point1=(3.25, 60.0), point2=(6.0, 60.0))
s.HorizontalConstraint(entity=g[5], addUndoState=False)
s.PerpendicularConstraint(entity1=g[4], entity2=g[5], addUndoState=False)
s.Line(point1=(6.0, 60.0), point2=(6.0, 23.95))
s.VerticalConstraint(entity=g[6], addUndoState=False)
s.PerpendicularConstraint(entity1=g[5], entity2=g[6], addUndoState=False)
s.ArcByStartEndTangent(point1=(4.25, 15.0), point2=(6.0, 23.95), entity=g[2])
mdb.models['Model-1'].sketches.changeKey(fromName='__profile__', 
    toName='Specimen_half')
s.unsetPrimaryObject()

s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.sketchOptions.setValues(viewStyle=AXISYM)
s.setPrimaryObject(option=STANDALONE)
s.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
s.FixedConstraint(entity=g[2])
s.retrieveSketch(sketch=mdb.models['Model-1'].sketches['Specimen_half'])
session.viewports['Viewport: 1'].view.fitView()
#: Info: 8 entities copied from Clamp.
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=AXISYMMETRIC, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']

p = mdb.models['Model-1'].parts['Part-1']
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=15.0)
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=23.95)
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=40.0)
f, e, d = p.faces, p.edges, p.datums

for i in [2,3,4]:
    p.PartitionFaceByDatumPlane(datumPlane=d[i], faces=f)

s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.Line(point1=(20.0, 40.0), point2=(20.0, 0.0))
s.VerticalConstraint(entity=g[2], addUndoState=False)
s.Line(point1=(20.0, 0.0), point2=(6.0, 0.0))
s.HorizontalConstraint(entity=g[3], addUndoState=False)
s.PerpendicularConstraint(entity1=g[2], entity2=g[3], addUndoState=False)
s.Line(point1=(6.0, 0.0), point2=(6.0, 20.0))
s.VerticalConstraint(entity=g[4], addUndoState=False)
s.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
s.Line(point1=(6.0, 20.0), point2=(3.25, 20.0))
s.HorizontalConstraint(entity=g[5], addUndoState=False)
s.PerpendicularConstraint(entity1=g[4], entity2=g[5], addUndoState=False)
s.Line(point1=(3.25, 20.0), point2=(3.25, 150.0))
s.VerticalConstraint(entity=g[6], addUndoState=False)
s.PerpendicularConstraint(entity1=g[5], entity2=g[6], addUndoState=False)
s.Line(point1=(3.25, 150.0), point2=(12.5, 150.0))
s.HorizontalConstraint(entity=g[7], addUndoState=False)
s.PerpendicularConstraint(entity1=g[6], entity2=g[7], addUndoState=False)
s.Line(point1=(12.5, 150.0), point2=(12.5, 60.0))
s.VerticalConstraint(entity=g[8], addUndoState=False)
s.PerpendicularConstraint(entity1=g[7], entity2=g[8], addUndoState=False)
s.ArcByStartEndTangent(point1=(12.5, 60.0), point2=(20.0, 40.0), entity=g[8])
mdb.models['Model-1'].sketches.changeKey(fromName='__profile__', 
    toName='Clamp')
s.unsetPrimaryObject()

s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.sketchOptions.setValues(viewStyle=AXISYM)
s.setPrimaryObject(option=STANDALONE)
s.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
s.FixedConstraint(entity=g[2])
s.retrieveSketch(sketch=mdb.models['Model-1'].sketches['Clamp'])
session.viewports['Viewport: 1'].view.fitView()
#: Info: 8 entities copied from Clamp.
p = mdb.models['Model-1'].Part(name='Part-2', dimensionality=AXISYMMETRIC, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-2']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-2']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']


p = mdb.models['Model-1'].parts['Part-2']
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=20.0)
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=40.0)
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=60.0)
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=100.0)
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=6.0)
f, e, d = p.faces, p.edges, p.datums

for i in [2,3,4,5,6]:
    p.PartitionFaceByDatumPlane(datumPlane=d[i], faces=f)

mdb.models['Model-1'].setValues(stefanBoltzmann=5.67e-8*1e-3)
mdb.models['Model-1'].setValues(absoluteZero=-273.15)






#==============================================================================
# create material and section
#==============================================================================
materialname='User'
mdb.models['Model-1'].Material(name=materialname)
mdb.models['Model-1'].materials[materialname].Depvar(n=200)
mdb.models['Model-1'].materials[materialname].UserMaterial(mechanicalConstants=(0.0, ))
#mdb.models['Model-1'].materials[materialname].Density(temperatureDependency=ON, table=((8.24E-9, 20), ))
mdb.models['Model-1'].materials[materialname].Conductivity(temperatureDependency=ON, table=((13.4, 11), (14.7, 100), (15.9, 200), (17.8, 300), (18.3, 400), (19.6, 500), (21.2, 600), (22.8, 700), (23.6, 800), (27.6, 900), (30.4, 1000), ))
#mdb.models['Model-1'].materials[materialname].Conductivity(temperatureDependency=ON, table=((10.0, 11), (14.7, 100), (15.9, 200), (17.8, 300), (18.3, 400), (19.6, 500), (21.2, 600), (22.8, 700), (23.6, 800), (27.6, 900), (30.4, 1000), ))
#mdb.models['Model-1'].materials[materialname].SpecificHeat(temperatureDependency=ON, table=((481.4E6, 300), (493.9E6, 400), (514.8E6, 500), (539E6, 600), (573.4E6, 700), (615.3E6, 800), (657.2E6, 900), (707.4E6, 1000), ))
#mdb.models['Model-1'].materials[materialname].Expansion(temperatureDependency=ON, table=((9.1E-6, -253), (11E-6, -183), (11.8E-6, 20), (11.8E-6, 100), (13E-6, 200), (13.5E-6, 300), (14.1E-6, 400), (14.4E-6, 500), (14.8E-6, 600), (15.4E-6, 700), (17E-6, 800), (18.4E-6, 900), (18.7E-6, 1000), ))
#mdb.models['Model-1'].materials[materialname].Expansion(table=((0.0e-06, ), ))
mdb.models['Model-1'].materials[materialname].Expansion(table=((7.5e-06, ), ))
mdb.models['Model-1'].HomogeneousSolidSection(name='Section-User', material=materialname, thickness=None)

materialname='GH4169'
mdb.models['Model-1'].Material(name=materialname)
mdb.models['Model-1'].materials[materialname].Density(temperatureDependency=ON, table=((8.24E-9, 20), ))
mdb.models['Model-1'].materials[materialname].Conductivity(temperatureDependency=ON, table=((13.4, 11), (14.7, 100), (15.9, 200), (17.8, 300), (18.3, 400), (19.6, 500), (21.2, 600), (22.8, 700), (23.6, 800), (27.6, 900), (30.4, 1000), ))
mdb.models['Model-1'].materials[materialname].SpecificHeat(temperatureDependency=ON, table=((481.4E6, 300), (493.9E6, 400), (514.8E6, 500), (539E6, 600), (573.4E6, 700), (615.3E6, 800), (657.2E6, 900), (707.4E6, 1000), ))
mdb.models['Model-1'].materials[materialname].Expansion(temperatureDependency=ON, table=((9.1E-6, -253), (11E-6, -183), (11.8E-6, 20), (11.8E-6, 100), (13E-6, 200), (13.5E-6, 300), (14.1E-6, 400), (14.4E-6, 500), (14.8E-6, 600), (15.4E-6, 700), (17E-6, 800), (18.4E-6, 900), (18.7E-6, 1000), ))
mdb.models['Model-1'].materials[materialname].Elastic(temperatureDependency=ON, table=((204E3, 0.3, 20), (181E3, 0.3, 300), (176E3, 0.31, 400), (160E3, 0.32, 500), (160E3, 0.32, 550), (150E3, 0.32, 600), (146E3, 0.325, 650), (141E3, 0.33, 700), ))
mdb.models['Model-1'].HomogeneousSolidSection(name='Section-GH4169', material=materialname, thickness=None)






#==============================================================================
# assign section
#==============================================================================
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#f ]', ), )
region = p.Set(faces=faces, name='Set-1')
p = mdb.models['Model-1'].parts['Part-1']
p.SectionAssignment(region=region, sectionName='Section-GH4169', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

p = mdb.models['Model-1'].parts['Part-2']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#1ff ]', ), )
region = regionToolset.Region(faces=faces)
p = mdb.models['Model-1'].parts['Part-2']
p.SectionAssignment(region=region, sectionName='Section-GH4169', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)









#==============================================================================
# assembly
#==============================================================================
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(coordSysType=CYLINDRICAL, origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), point2=(0.0, 0.0, -1.0))
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1-1', part=p, dependent=ON)

a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Part-2']
a.Instance(name='Part-2-1', part=p, dependent=ON)
a.translate(instanceList=('Part-2-1', ), vector=(0.0, 40.0, 0.0))
#: The instance Part-2-1 was translated by 0., 40., 0. with respect to the assembly coordinate system










#==============================================================================
# mesh
#==============================================================================
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, engineeringFeatures=OFF, mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(meshTechnique=ON)
p1 = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
# seed edges 
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#f ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#1fff ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=5, constraint=FINER)

# element type
elemType1 = mesh.ElemType(elemCode=DCAX8, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=DCAX6, elemLibrary=STANDARD)

p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
p.setMeshControls(regions=f, technique=STRUCTURED)
faces = f.getSequenceFromMask(mask=('[#f ]', ), )
pickedRegions =(faces, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
# generate mesh
p.generateMesh()

p = mdb.models['Model-1'].parts['Part-2']
f, e, d = p.faces, p.edges, p.datums
p.setMeshControls(regions=f, technique=STRUCTURED)
p.seedPart(size=5.0, deviationFactor=0.1, minSizeFactor=0.1)
pickedEdges = e.getSequenceFromMask(mask=('[#1ffffff ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=2, constraint=FINER)
pickedEdges = e.getSequenceFromMask(mask=('[#5a0055 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=4, constraint=FINER)
elemType1 = mesh.ElemType(elemCode=DCAX8, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=DCAX6, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['Part-2']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#1ff ]', ), )
pickedRegions =(faces, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
p.generateMesh()








#==============================================================================
# node set
#==============================================================================
p = mdb.models['Model-1'].parts['Part-1']
n = p.nodes
nodes = n.getSequenceFromMask(mask=('[#ffffffff:4 #3ffff ]', ), )
p.Set(nodes=nodes, name='NALL')
p = mdb.models['Model-1'].parts['Part-1']
p.SetFromNodeLabels('N36',(7,))
p.SetFromNodeLabels('N33',(8,))
p.SetFromNodeLabels('N92',(252,))
#: The set 'N36' has been created (1 node).








#==============================================================================
# step
#==============================================================================
#mdb.models['Model-1'].CoupledTempDisplacementStep(name='Step-1', 
#    previous='Initial', response=STEADY_STATE, timePeriod=time_period,
#    maxNumInc=1000000, initialInc=initial_inc, minInc=min_inc, maxInc=max_inc, deltmx=None, 
#    cetol=None, creepIntegration=None, amplitude=RAMP
#    , nlgeom=nonlinear
#    )

mdb.models['Model-1'].HeatTransferStep(name='Step-1', previous='Initial', 
    response=STEADY_STATE, maxNumInc=1000, initialInc=0.1, amplitude=RAMP)
mdb.models['Model-1'].steps['Step-1'].setValues(response=TRANSIENT, 
    timePeriod=600.0, minInc=0.0006, maxInc=60.0, deltmx=10.0, amplitude=STEP)







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

# Displacement Triangular Wave 1
mdb.models['Model-1'].TabularAmplitude(name='Displacement_1', timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (60.0, 1.0*0.0099)))
# Displacement Triangular Wave 2
mdb.models['Model-1'].TabularAmplitude(name='Displacement_2', timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (60.0, 0.0*0.0099)))
# Displacement Triangular Wave 3
mdb.models['Model-1'].TabularAmplitude(name='Displacement_3', timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (60.0, -1.0*0.0099)))
# Displacement Triangular Wave 4
mdb.models['Model-1'].TabularAmplitude(name='Displacement_4', timeSpan=STEP, smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (60.0, 0.0*0.0099)))

mdb.models['Model-1'].TabularAmplitude(name='Amp-Heatflux', timeSpan=STEP, 
    smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (5.0, 0.0), (6, 1.0), (31.0, 1.0), (32, 0.0), (
    10000.0, 0.0)))

ambient_temperature = read_json_file('ambient_temperature.txt')

mdb.models['Model-1'].TabularAmplitude(name='Amp-Ambient-Temperature', timeSpan=STEP, 
    smooth=SOLVER_DEFAULT, data=ambient_temperature)





#==============================================================================
# predifined field
#==============================================================================
#a = mdb.models['Model-1'].rootAssembly
#region = a.instances['Part-1-1'].sets['Set-1']
#mdb.models['Model-1'].Temperature(name='Predefined Field-1', createStepName='Initial', region=region, distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(predefined_temperature, ))

a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-1-1'].faces
faces1 = f1.getSequenceFromMask(mask=('[#f ]', ), )
e1 = a.instances['Part-1-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#1fff ]', ), )
v1 = a.instances['Part-1-1'].vertices
verts1 = v1.getSequenceFromMask(mask=('[#3ff ]', ), )
f2 = a.instances['Part-2-1'].faces
faces2 = f2.getSequenceFromMask(mask=('[#1ff ]', ), )
e2 = a.instances['Part-2-1'].edges
edges2 = e2.getSequenceFromMask(mask=('[#1ffffff ]', ), )
v2 = a.instances['Part-2-1'].vertices
verts2 = v2.getSequenceFromMask(mask=('[#1ffff ]', ), )
region = regionToolset.Region(vertices=verts1+verts2, edges=edges1+edges2, 
    faces=faces1+faces2)
mdb.models['Model-1'].Temperature(name='Predefined Field-1', createStepName='Initial', region=region, distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(predefined_temperature, ))









#==============================================================================
# boundary conditions
#==============================================================================
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['Part-1-1'].edges
side1Edges1 = s1.getSequenceFromMask(mask=('[#248 ]', ), )
region = a.Surface(side1Edges=side1Edges1, name='Surf-Outer')
mdb.models['Model-1'].SurfaceHeatFlux(name='Outer-Heatflux', 
    createStepName='Step-1', region=region, magnitude=heat_flux, amplitude='Amp-Heatflux')
mdb.models['Model-1'].ExpressionField(name='AnalyticalField-1', localCsys=None, 
    description='', expression='3.25 /  X ')
mdb.models['Model-1'].loads['Outer-Heatflux'].setValues(distributionType=FIELD, 
    field='AnalyticalField-1')

#a = mdb.models['Model-1'].rootAssembly
#e1 = a.instances['Part-1-1'].edges
#edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
#region = regionToolset.Region(edges=edges1)
#mdb.models['Model-1'].TemperatureBC(name='BC-Outer-Temperature', createStepName='Step-1', 
#    region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
#    magnitude=outer_temperature, amplitude=UNSET)

a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['Part-1-1'].edges
side1Edges1 = s1.getSequenceFromMask(mask=('[#1092 ]', ), )
region=a.Surface(side1Edges=side1Edges1, name='Surf-Inner')
mdb.models['Model-1'].FilmCondition(name='Int-Inner-Convection', 
    createStepName='Step-1', surface=region, definition=EMBEDDED_COEFF, 
    filmCoeff=film_coefficient_inner, filmCoeffAmplitude='', sinkTemperature=sink_temperature_inner, 
    sinkAmplitude='', sinkDistributionType=UNIFORM, sinkFieldName='')

a = mdb.models['Model-1'].rootAssembly
region=a.surfaces['Surf-Outer']
mdb.models['Model-1'].FilmCondition(name='Int-Outer-Convection', 
    createStepName='Step-1', surface=region, definition=EMBEDDED_COEFF, 
    filmCoeff=film_coefficient_outer, filmCoeffAmplitude='', sinkTemperature=1.0, 
    sinkAmplitude='Amp-Ambient-Temperature', sinkDistributionType=UNIFORM, sinkFieldName='')

a = mdb.models['Model-1'].rootAssembly
region=a.surfaces['Surf-Outer']
mdb.models['Model-1'].RadiationToAmbient(name='Int-Outer-Radiation', createStepName='Step-1', 
    surface=region, radiationType=AMBIENT, distributionType=UNIFORM, field='', 
    emissivity=emissivity, ambientTemperature=1.0, ambientTemperatureAmp='Amp-Ambient-Temperature')

mdb.models['Model-1'].ContactProperty('IntProp-Conduction')
mdb.models['Model-1'].interactionProperties['IntProp-Conduction'].ThermalConductance(
    definition=TABULAR, clearanceDependency=ON, pressureDependency=OFF, 
    temperatureDependencyC=OFF, massFlowRateDependencyC=OFF, dependenciesC=0, 
    clearanceDepTable=((0.5, 0.0), (0.0, 1.0)))
#: The interaction property "IntProp-1" has been created.
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['Part-1-1'].edges
side1Edges1 = s1.getSequenceFromMask(mask=('[#400 ]', ), )
region1=regionToolset.Region(side1Edges=side1Edges1)
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['Part-2-1'].edges
side1Edges1 = s1.getSequenceFromMask(mask=('[#100000 ]', ), )
region2=regionToolset.Region(side1Edges=side1Edges1)
mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Int-Conduction', 
    createStepName='Initial', master=region1, slave=region2, sliding=FINITE, 
    thickness=ON, interactionProperty='IntProp-Conduction', adjustMethod=NONE, 
    initialClearance=OMIT, datumAxis=None, clearanceRegion=None)

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-2-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#200 ]', ), )
region = regionToolset.Region(edges=edges1)
mdb.models['Model-1'].TemperatureBC(name='BC-1', createStepName='Step-1', 
    region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    magnitude=20.0, amplitude=UNSET)

#: The interaction "Int-4" has been created.

#a = mdb.models['Model-1'].rootAssembly
#region=a.surfaces['Surf-Inner']
#mdb.models['Model-1'].RadiationToAmbient(name='Int-Inner-Radiation', createStepName='Step-1', 
#    surface=region, radiationType=AMBIENT, distributionType=UNIFORM, field='', 
#    emissivity=emissivity, ambientTemperature=sink_temperature_inner, ambientTemperatureAmp='')

#a = mdb.models['Model-1'].rootAssembly
#s1 = a.instances['Part-1-1'].edges
#side1Edges1 = s1.getSequenceFromMask(mask=('[#4 ]', ), )
#region = a.Surface(side1Edges=side1Edges1, name='Surf-Top')
#mdb.models['Model-1'].SurfaceHeatFlux(name='Top-Heatflux', createStepName='Step-1', 
#    region=region, magnitude=-1000.0)




#==============================================================================
# field output
#==============================================================================
#mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'E', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 'NT', 'TEMP', 'HFL', 'RFL'))
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=('NT', 'TEMP', 'HFL', 'RFL'))
#==============================================================================
# history output
#==============================================================================
#mdb.models['Model-1'].historyOutputRequests['H-Output-1'].setValues(variables=PRESELECT)
mdb.models['Model-1'].HistoryOutputRequest(name='H-Output-1', createStepName='Step-1', variables=('FTEMP', 'HFLA', 'HTL'))
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