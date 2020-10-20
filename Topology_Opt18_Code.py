# Serena Ferraro
# ------------------------------------------------------------------------
# Do not delete the following import lines
import sys, os, glob
import numpy as np
from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
from textRepr import prettyPrint
from numpy.oldnumeric import array, Int32, Float64
from numpy import linalg
import time
from os.path import exists
	


def runFullAnalysis(xtemp,ytemp,ModelName,JobName):

	# CPU
	NumDomains = 8

	#-----------------------------------------------------------
	# Faiulre strength parameters 2-ply laminates (Units: N-mm)
	#-----------------------------------------------------------
	sF1t = 76.16
	sF1c = 34.50
	sF3 = 14.55
	sF4 = 3.26
	sF6 = 1.10
	#-----------------------------------------------------------
	# Faiulre strength parameters 4-ply laminates (Units: N-mm)
	#-----------------------------------------------------------
	# sF1t = 148.85
	# sF1c = 54.87
	# sF3 = 14.27
	# sF4 = 6.84
	# sF6 = 3.49
	#
	# calculating failure coefficients
	#---------------------------------
	sJ1 = 1/sF1t - 1/sF1c
	sK11 = 1/(sF1t*sF1c)
	sK33 = 1/pow(sF3,2)
	sK44 = 1/pow(sF4,2)
	sK66 = 1/pow(sF6,2)
	sK12 = -sK11/2


	#-------------------------------------------------------------------------
	#----------------------Create Parts---------------------------------------
	#-------------------------------------------------------------------------
	#def createParts():
	mdb.Model(ModelName)

	# central reference point: RP_center
	RP_center = mdb.models[ModelName].Part(name='RP_center', dimensionality=THREE_D, 
		type=DEFORMABLE_BODY)
	RP_center.ReferencePoint(point=(0.0, 130.1960361, 0))
	RP_center = mdb.models[ModelName].parts['RP_center']
	# Left reference point: RP_left
	RP_left = mdb.models[ModelName].Part(name='RP_left', dimensionality=THREE_D, 
		type=DEFORMABLE_BODY)
	RP_left.ReferencePoint(point=(0.0, 0.0, 130.1960361))
	RP_left = mdb.models[ModelName].parts['RP_left']
	# Right reference point: RP_right
	RP_right = mdb.models[ModelName].Part(name='RP_right', dimensionality=THREE_D, 
		type=DEFORMABLE_BODY)
	RP_right.ReferencePoint(point=(0.0, 0.0, -130.1960361))
	RP_right = mdb.models[ModelName].parts['RP_right']
	# Flexural joint geometry imported from CAD model
	acis = mdb.openAcis(
		'D:/Serena/Repeat_Topology_Studies/Opt18_cutouts_001-09/smooth_joint_U5_L3.SAT', 
		scaleFromFile=OFF)
	mdb.models[ModelName].PartFromGeometryFile(name='smooth_joint_U5_L3', 
		geometryFile=acis, combine=True, dimensionality=THREE_D, 
		type=DEFORMABLE_BODY)
	joint = mdb.models[ModelName].parts['smooth_joint_U5_L3']

			
	#--------------------------------------------------------------------------
	#------------------Create 2,3,4,6 Plies Sets-------------------------------
	#--------------------------------------------------------------------------
	#def plySets():
	f = joint.faces

	# 2 plies
	faces = f.getSequenceFromMask(mask=('[#0 #738019bc #3efc1c ]', ), )
	joint.Set(faces=faces, name='two_plies')
	# 3 plies
	faces = f.getSequenceFromMask(mask=('[#d7bb272f #39c002 ]', ), )
	joint.Set(faces=faces, name='three_plies')
	# 4 plies
	faces = f.getSequenceFromMask(mask=('[#0 #8c000640 #c103e3 ]', ), )
	joint.Set(faces=faces, name='four_plies')
	# 6 plies
	faces = f.getSequenceFromMask(mask=('[#2844d8d0 #462001 ]', ), )
	joint.Set(faces=faces, name='six_plies')

			
	#--------------------------------------------------------------------------
	#-------------------Create and Assign Stiffness Matrices-------------------
	#--------------------------------------------------------------------------
	#def stiffnessMatrices():
	### ABD Matrices
	mdb.models[ModelName].GeneralStiffnessSection(name='Astroquartz- [45/0/45]', 
		referenceTemperature=None, stiffnessMatrix=(5640.9, 2121.9, 5640.9, 
		0.0, 0.0, 2402.1, 0.0, 0.0, 0.0, 10.0013, 0.0, 0.0, 0.0, 5.31263, 
		10.0013, 0.0, 0.0, 0.0, 0.0, 0.0, 5.86543), applyThermalStress=0, 
		poissonDefinition=DEFAULT, useDensity=ON, density=3.2e-10)
	mdb.models[ModelName].GeneralStiffnessSection(
		name='Astroquartz- [45/45/0/0/45/45]', referenceTemperature=None, 
		stiffnessMatrix=(9565.1, 3598.0, 9565.1, 0.0, 0.0, 4073.1, 0.0, 0.0, 
		0.0, 69.7918, 0.0, 0.0, 0.0, 37.073, 69.7918, 0.0, 0.0, 0.0, 0.0, 0.0, 
		40.9306), applyThermalStress=0, poissonDefinition=DEFAULT, 
		useDensity=ON, density=6.1e-10)
	mdb.models[ModelName].GeneralStiffnessSection(
		name='Astroquartz- [45/45/45/45]', referenceTemperature=None, 
		stiffnessMatrix=(6085.0, 3365.4, 6085.0, 0.0, 0.0, 3706.5, 0.0, 0.0, 
		0.0, 23.0122, 0.0, 0.0, 0.0, 12.7273, 23.0122, 0.0, 0.0, 0.0, 0.0, 0.0, 
		14.0174), applyThermalStress=0, poissonDefinition=DEFAULT, 
		useDensity=ON, density=4.15e-10)
	mdb.models[ModelName].GeneralStiffnessSection(name='Astroquartz- [45/45]', 
		referenceTemperature=None, stiffnessMatrix=(3477.1, 1923.1, 3477.1, 
		0.0, 0.0, 2118.0, 0.0, 0.0, 0.0, 3.09537, 0.0, 0.0, 0.0, 1.71087, 
		3.09537, 0.0, 0.0, 0.0, 0.0, 0.0, 1.87778), applyThermalStress=0, 
		poissonDefinition=DEFAULT, useDensity=ON, density=2.1e-10)
	### Section Assignment	
	region = joint.sets['four_plies']
	joint.SectionAssignment(region=region, sectionName='Astroquartz- [45/45/45/45]', 
		offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', 
		thicknessAssignment=FROM_SECTION)
	region = joint.sets['six_plies']
	joint.SectionAssignment(region=region, 
		sectionName='Astroquartz- [45/45/0/0/45/45]', offset=0.0, 
		offsetType=MIDDLE_SURFACE, offsetField='', 
		thicknessAssignment=FROM_SECTION)
	region = joint.sets['three_plies']
	joint.SectionAssignment(region=region, sectionName='Astroquartz- [45/0/45]', 
		offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', 
		thicknessAssignment=FROM_SECTION)
	region = joint.sets['two_plies']
	joint.SectionAssignment(region=region, sectionName='Astroquartz- [45/45]', 
		offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', 
		thicknessAssignment=FROM_SECTION)

	#---------------------------------------------------------------------------
	#---------Create Datum Coordinate Sistems and Local Orientations------------
	#--------------------------------------------------------------------------- 
	#def localOrientation():
	# XYZ_left
	joint.DatumCsysByThreePoints(name='XYZ_left', coordSysType=CARTESIAN, 
			origin=(0.0, 0.0, 130.1960361), point1=(1.0, 0.0, 130.1960361), point2=(0.0, 11.22532015, 141.4213563))	

	faces = f.getSequenceFromMask(mask=('[#1fc03fc3 #8303ea8c #5e057f ]', ), )
	region = regionToolset.Region(faces=faces)
	d = joint.datums
	orientation = d[d.keys()[-1]]
	joint.MaterialOrientation(region=region, orientationType=SYSTEM, axis=AXIS_2, 
			localCsys=orientation, fieldName='', additionalRotationType=ROTATION_ANGLE, additionalRotationField='', 
			angle=180.0)		
			
	# XYZ_right
	joint.DatumCsysByThreePoints(name='XYZ_right', coordSysType=CARTESIAN, 
			origin=(0.0, 0.0, -130.1960361), point1=(1.0, 0.0, -130.1960361), point2=(0.0, 11.22532015, -141.4213563))	
			
	faces = f.getSequenceFromMask(mask=('[#e03fc03c #7cfc1573 #a1fa80 ]', ), )
	region = regionToolset.Region(faces=faces)
	d = joint.datums
	orientation = d[d.keys()[-1]]
	joint.MaterialOrientation(region=region, orientationType=SYSTEM, axis=AXIS_2, localCsys=orientation, fieldName='', 
			additionalRotationType=ROTATION_NONE, angle=0.0, additionalRotationField='')		
		
	#----------------------------------------------------------------------------
	#--------------------------------Create Cutout-------------------------------
	#----------------------------------------------------------------------------
	#def createCutout():
	# Sketch Plane - cutout_plane
	joint.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=15.875)
	joint.features.changeKey(
		fromName='Datum plane-1', toName='cutout_plane')
		
	# Cutout
	# spot are points that I choose on sketch.
	e, d = joint.edges, joint.datums
	cutoutPlane = d[d.keys()[-1]]
	t = joint.MakeSketchTransform(sketchPlane=d[d.keys()[-1]], sketchUpEdge=e[178],
			sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(15.875, 111.0, 0.0))
	sketch1 = mdb.models[ModelName].ConstrainedSketch(name='__profile__',
			sheetSize=720.01, gridSpacing=18.0, transform=t)
	g, v, dim, c = sketch1.geometry, sketch1.vertices, sketch1.dimensions, sketch1.constraints
	sketch1.setPrimaryObject(option=SUPERIMPOSE)
	joint.projectReferencesOntoSketch(sketch=sketch1, filter=COPLANAR_EDGES)

	pts = []
	for II in range(0,len(xtemp)):
		for JJ in range(0,len(xtemp[II])):
			
			sketch1.Spot(point=[xtemp[II][JJ][0],ytemp[II][JJ][0]])
			pts.append([xtemp[II][JJ][0],ytemp[II][JJ][0]])
		
		pts.append([xtemp[II][0][0],ytemp[II][0][0]])
		Pts = tuple(pts)
		# spline
		sketch1.Spline(points=Pts)
		pts=[]

	# project cutout on joint
	joint.CutExtrude(sketchPlane=d[d.keys()[-1]], sketchUpEdge=e[178], sketchPlaneSide=SIDE1, 
			sketchOrientation=RIGHT, sketch=sketch1, flipExtrudeDirection=OFF)
	sketch1.unsetPrimaryObject()

	# #-------------------------------------------------------------------------------
	# #----------------------------Create Mesh----------------------------------------
	# #-------------------------------------------------------------------------------	
	#def meshPart():

	# part seed
	joint.seedPart(size=2.0, deviationFactor=0.1, minSizeFactor=0.1)
	e = joint.edges

	# S4R 0.125
	one = e.getByBoundingBox(-20.0,90.0,-24.0,20.0,200.0,24.0)
	joint.seedEdgeBySize(edges=one, size=0.125, deviationFactor=0.1, minSizeFactor=0.1, constraint=FINER)

	two = e.getClosest(coordinates=((12.5,114.031003,29.973683),))
	joint.seedEdgeBySize(edges=(two[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	three=e.getClosest(coordinates=((7.500017,127.761975,22.192651),))
	joint.seedEdgeBySize(edges=(three[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	four=e.getClosest(coordinates=((0.,129.665706,22.98097),))
	joint.seedEdgeBySize(edges=(four[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	five=e.getClosest(coordinates=((-7.5,127.761981,22.192657),))
	joint.seedEdgeBySize(edges=(five[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	six=e.getClosest(coordinates=((-12.5,114.031003,29.973683),))
	joint.seedEdgeBySize(edges=(six[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	seven=e.getClosest(coordinates=((-12.5,100.211879,16.154558),))
	joint.seedEdgeBySize(edges=(seven[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	eight=e.getClosest(coordinates=((-7.5,98.596174,11.832451),))
	joint.seedEdgeBySize(edges=(eight[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	nine=e.getClosest(coordinates=((0.,96.696852,11.048543),))
	joint.seedEdgeBySize(edges=(nine[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	ten=e.getClosest(coordinates=((7.500011,98.596179,11.832456),))
	joint.seedEdgeBySize(edges=(ten[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	eleven=e.getClosest(coordinates=((12.5,100.211879,16.154558),))
	joint.seedEdgeBySize(edges=(eleven[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)

	two = e.getClosest(coordinates=((12.5,114.031003,-29.973683),))
	joint.seedEdgeBySize(edges=(two[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	three=e.getClosest(coordinates=((7.500017,127.761975,-22.192651),))
	joint.seedEdgeBySize(edges=(three[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	four=e.getClosest(coordinates=((0.,129.665706,-22.98097),))
	joint.seedEdgeBySize(edges=(four[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	five=e.getClosest(coordinates=((-7.5,127.761981,-22.192657),))
	joint.seedEdgeBySize(edges=(five[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	six=e.getClosest(coordinates=((-12.5,114.031003,-29.973683),))
	joint.seedEdgeBySize(edges=(six[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	seven=e.getClosest(coordinates=((-12.5,100.211879,-16.154558),))
	joint.seedEdgeBySize(edges=(seven[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	eight=e.getClosest(coordinates=((-7.5,98.596174,-11.832451),))
	joint.seedEdgeBySize(edges=(eight[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	nine=e.getClosest(coordinates=((0.,96.696852,-11.048543),))
	joint.seedEdgeBySize(edges=(nine[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	ten=e.getClosest(coordinates=((7.500011,98.596179,-11.832456),))
	joint.seedEdgeBySize(edges=(ten[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)
	eleven=e.getClosest(coordinates=((12.5,100.211879,-16.154558),))
	joint.seedEdgeBySize(edges=(eleven[0][0], ), size=0.125, deviationFactor=0.1, minSizeFactor=0.1)

	# S4R 0.7
	two = e.getClosest(coordinates=((12.5,92.088687,24.294825),))
	joint.seedEdgeBySize(edges=(two[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	three=e.getClosest(coordinates=((15.874999,98.991558,31.197696),))
	joint.seedEdgeBySize(edges=(three[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	four=e.getClosest(coordinates=((12.5,100.587144,43.399883),))
	joint.seedEdgeBySize(edges=(four[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	five=e.getClosest(coordinates=((7.500021,103.560735,46.373474),))
	joint.seedEdgeBySize(edges=(five[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	six=e.getClosest(coordinates=((10.220416,107.564949,39.771086),))
	joint.seedEdgeBySize(edges=(six[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	seven=e.getClosest(coordinates=((0.,104.916969,47.729708),))
	joint.seedEdgeBySize(edges=(seven[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	eight=e.getClosest(coordinates=((3.866397,109.866944,42.073081),))
	joint.seedEdgeBySize(edges=(eight[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	nine=e.getClosest(coordinates=((-7.5,103.560743,46.373482),))
	joint.seedEdgeBySize(edges=(nine[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	ten=e.getClosest(coordinates=((-3.866393,109.866944,42.073082),))
	joint.seedEdgeBySize(edges=(ten[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	eleven=e.getClosest(coordinates=((-12.5,100.587144,43.399883),))
	joint.seedEdgeBySize(edges=(eleven[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	twelve=e.getClosest(coordinates=((-10.220396,107.564961,39.771098),))
	joint.seedEdgeBySize(edges=(twelve[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	thirteen=e.getClosest(coordinates=((-15.874999,98.991548,31.197686),))
	joint.seedEdgeBySize(edges=(thirteen[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	fourteen=e.getClosest(coordinates=((-12.5,86.790022,29.60276),))
	joint.seedEdgeBySize(edges=(fourteen[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	fifteen=e.getClosest(coordinates=((-10.22178,90.416504,22.622641),))
	joint.seedEdgeBySize(edges=(fifteen[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	sixteen=e.getClosest(coordinates=((-7.5,83.817008,26.629747),))
	joint.seedEdgeBySize(edges=(sixteen[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	seventeen=e.getClosest(coordinates=((-3.865661,88.120436,20.326573),))
	joint.seedEdgeBySize(edges=(seventeen[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	eighteen=e.getClosest(coordinates=((0.,82.466328,25.279067),))
	joint.seedEdgeBySize(edges=(eighteen[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	nineteen=e.getClosest(coordinates=((3.865576,88.12042,20.326558),))
	joint.seedEdgeBySize(edges=(nineteen[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	twenty=e.getClosest(coordinates=((7.500006,83.81701,26.629749),))
	joint.seedEdgeBySize(edges=(twenty[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	twentyone=e.getClosest(coordinates=((10.221751,90.416487,22.622624),))
	joint.seedEdgeBySize(edges=(twentyone[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)

	two = e.getClosest(coordinates=((12.5,92.088687,-24.294825),))
	joint.seedEdgeBySize(edges=(two[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	three=e.getClosest(coordinates=((15.874999,98.991558,-31.197696),))
	joint.seedEdgeBySize(edges=(three[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	four=e.getClosest(coordinates=((12.5,100.587144,-43.399883),))
	joint.seedEdgeBySize(edges=(four[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	five=e.getClosest(coordinates=((7.500021,103.560735,-46.373474),))
	joint.seedEdgeBySize(edges=(five[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	six=e.getClosest(coordinates=((10.220416,107.564949,-39.771086),))
	joint.seedEdgeBySize(edges=(six[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	seven=e.getClosest(coordinates=((0.,104.916969,-47.729708),))
	joint.seedEdgeBySize(edges=(seven[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	eight=e.getClosest(coordinates=((3.866397,109.866944,-42.073081),))
	joint.seedEdgeBySize(edges=(eight[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	nine=e.getClosest(coordinates=((-7.5,103.560743,-46.373482),))
	joint.seedEdgeBySize(edges=(nine[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	ten=e.getClosest(coordinates=((-3.866393,109.866944,-42.073082),))
	joint.seedEdgeBySize(edges=(ten[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	eleven=e.getClosest(coordinates=((-12.5,100.587144,-43.399883),))
	joint.seedEdgeBySize(edges=(eleven[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	twelve=e.getClosest(coordinates=((-10.220396,107.564961,-39.771098),))
	joint.seedEdgeBySize(edges=(twelve[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	thirteen=e.getClosest(coordinates=((-15.874999,98.991548,-31.197686),))
	joint.seedEdgeBySize(edges=(thirteen[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	fourteen=e.getClosest(coordinates=((-12.5,86.790022,-29.60276),))
	joint.seedEdgeBySize(edges=(fourteen[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	fifteen=e.getClosest(coordinates=((-10.22178,90.416504,-22.622641),))
	joint.seedEdgeBySize(edges=(fifteen[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	sixteen=e.getClosest(coordinates=((-7.5,83.817008,-26.629747),))
	joint.seedEdgeBySize(edges=(sixteen[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	seventeen=e.getClosest(coordinates=((-3.865661,88.120436,-20.326573),))
	joint.seedEdgeBySize(edges=(seventeen[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	eighteen=e.getClosest(coordinates=((0.,82.466328,-25.279067),))
	joint.seedEdgeBySize(edges=(eighteen[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	nineteen=e.getClosest(coordinates=((3.865576,88.12042,-20.326558),))
	joint.seedEdgeBySize(edges=(nineteen[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	twenty=e.getClosest(coordinates=((7.500006,83.81701,-26.629749),))
	joint.seedEdgeBySize(edges=(twenty[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)
	twentyone=e.getClosest(coordinates=((10.221751,90.416487,-22.622624),))
	joint.seedEdgeBySize(edges=(twentyone[0][0], ), size=0.7, deviationFactor=0.1, minSizeFactor=0.1)

	# mesh control
	f = joint.faces
	pickedFaces = f.getByBoundingBox(-20.0,-20.0,-150.0,20.0,200.0,150.0)
	joint.setMeshControls(regions=pickedFaces, elemShape=QUAD)

	# element type
	elemType1 = mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD)
	elemType2 = mesh.ElemType(elemCode=S3, elemLibrary=STANDARD, secondOrderAccuracy=OFF)
	pickedFaces =(pickedFaces, )
	joint.setElementType(regions=pickedFaces, elemTypes=(elemType1, elemType2))
	joint.generateMesh()
		
	#-------------------------------------------------------------------------------
	#----------------------------Create Assembly------------------------------------
	#-------------------------------------------------------------------------------	
	#def createAssembly():
	a = mdb.models[ModelName].rootAssembly
	a.DatumCsysByDefault(CARTESIAN)

	# generate features
	a.Instance(name='RP_center-1', part=RP_center, dependent=ON)
	a.Instance(name='RP_left-1', part=RP_left, dependent=ON)
	a.Instance(name='RP_right-1', part=RP_right, dependent=ON)
	a.Instance(name='smooth_joint_U5_L3-1', part=joint, dependent=ON)

	# copy sets from part module for reference poits
	r1 = a.instances['RP_center-1'].referencePoints
	refPoints1=(r1[1], )
	a.Set(referencePoints=refPoints1, name='RP_center')
	r1 = a.instances['RP_left-1'].referencePoints
	refPoints1=(r1[1], )
	a.Set(referencePoints=refPoints1, name='RP_left')
	r1 = a.instances['RP_right-1'].referencePoints
	refPoints1=(r1[1], )
	a.Set(referencePoints=refPoints1, name='RP_right')

	# generate sets for constraints and BCs
	v1 = a.instances['smooth_joint_U5_L3-1'].vertices
	verts1 = v1.findAt(((0.,150.574628,0.),))
	a.Set(vertices=verts1, name='fixed_top_point')
	verts1 = v1.findAt(((0.,-11.22532,118.970716),))
	a.Set(vertices=verts1, name='fixed_left_point')
	verts1 = v1.findAt(((0.,-11.22532,-118.970716),))
	a.Set(vertices=verts1, name='fixed_right_point')
	f = a.instances['smooth_joint_U5_L3-1'].faces
	face = f.getByBoundingSphere((15,0,-129),50)
	a.Set(faces=face, name='patch_right')
	face = f.getByBoundingSphere((15,0,129),50)
	a.Set(faces=face, name='patch_left')

	#-------------------------------------------------------------------------------
	#----------------------------Create Coupling Constraints------------------------
	#-------------------------------------------------------------------------------	
	#def createConstraints():
	region1=a.sets['RP_left']
	region2=a.sets['patch_left']
	mdb.models[ModelName].Coupling(name='coupling_leftpatch_to_RP', 
		controlPoint=region1, surface=region2, influenceRadius=WHOLE_SURFACE, 
		couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=OFF, u3=OFF, ur1=ON, 
		ur2=OFF, ur3=OFF)
	region1=a.sets['RP_right']
	region2=a.sets['patch_right']
	mdb.models[ModelName].Coupling(name='coupling_rightpatch_to_RP', 
		controlPoint=region1, surface=region2, influenceRadius=WHOLE_SURFACE, 
		couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=OFF, u3=OFF, ur1=ON, 
		ur2=OFF, ur3=OFF)
	mdb.models[ModelName].Equation(name='connect_rotations', terms=((-1.0, 'RP_left', 
		4), (1.0, 'RP_right', 4), (1.0, 'RP_center', 4)))
		
	#-------------------------------------------------------------------------------
	#----------------------------Create Amplitudes----------------------------------
	#-------------------------------------------------------------------------------	
	mdb.models[ModelName].SmoothStepAmplitude(name='amplitude_rotation', 
		timeSpan=STEP, data=((0.0, 0.0), (1.9, 1.0), (2.0, 1.0)))
	mdb.models[ModelName].SmoothStepAmplitude(name='amplitude_viscous_pressure', 
		timeSpan=STEP, data=((0.0, 1.0), (1.9, 1.0), (2.0, 0.0)))

	#-------------------------------------------------------------------------------
	#----------------------------Create Steps---------------------------------------
	#-------------------------------------------------------------------------------	
	#def defineStep():
	### Step 1 = rotation
	mdb.models[ModelName].StaticStep(name='rotation', previous='Initial', 
			timePeriod=2.0, stabilizationMethod=DISSIPATED_ENERGY_FRACTION, 
			continueDampingFactors=True, adaptiveDampingRatio=0.05, initialInc=2.0, 
			minInc=2e-05, maxInc=2.0, nlgeom=ON)

	# Field Output
	mdb.models[ModelName].fieldOutputRequests['F-Output-1'].setValues(variables=('S', 
			'SE', 'U', 'SF'))
	mdb.models[ModelName].fieldOutputRequests.changeKey(fromName='F-Output-1', 
		toName='Whole_Model')

	# History Output
	mdb.models[ModelName].historyOutputRequests.changeKey(fromName='H-Output-1', 
		toName='Energy')
	mdb.models[ModelName].historyOutputRequests['Energy'].setValues(
			variables=PRESELECT)
	regionDef=mdb.models[ModelName].rootAssembly.sets['RP_center']
	mdb.models[ModelName].HistoryOutputRequest(name='RP_center', 
			createStepName='rotation', variables=('UR', 'RM'), region=regionDef, 
			sectionPoints=DEFAULT, rebar=EXCLUDE)

	# Boundary Conditions
	region = a.sets['RP_center']
	mdb.models[ModelName].DisplacementBC(name='RP_center', createStepName='rotation', 
		region=region, u1=0.0, u2=0.0, u3=0.0, ur1=0.8, ur2=0.0, ur3=0.0, 
		amplitude='amplitude_rotation', fixed=OFF, distributionType=UNIFORM, 
		fieldName='', localCsys=None)
	region = a.sets['RP_left']
	mdb.models[ModelName].DisplacementBC(name='RP_left', createStepName='rotation', 
		region=region, u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=0.0, ur3=0.0, 
		amplitude='amplitude_rotation', fixed=OFF, distributionType=UNIFORM, 
		fieldName='', localCsys=None)
	region = a.sets['RP_right']
	mdb.models[ModelName].DisplacementBC(name='RP_right', createStepName='rotation', 
		region=region, u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=0.0, ur3=0.0, 
		amplitude='amplitude_rotation', fixed=OFF, distributionType=UNIFORM, 
		fieldName='', localCsys=None)
	region = a.sets['fixed_left_point']
	mdb.models[ModelName].DisplacementBC(name='fix_Ux_left', 
		createStepName='rotation', region=region, u1=0.0, u2=UNSET, u3=UNSET, 
		ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='amplitude_rotation', 
		fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
	region = a.sets['fixed_right_point']
	mdb.models[ModelName].DisplacementBC(name='fix_Ux_right', 
		createStepName='rotation', region=region, u1=0.0, u2=UNSET, u3=UNSET, 
		ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='amplitude_rotation', 
		fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
	region = a.sets['fixed_top_point']
	mdb.models[ModelName].DisplacementBC(name='fix_all_translations_top', 
		createStepName='rotation', region=region, u1=0.0, u2=0.0, u3=0.0, 
		ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude='amplitude_rotation', 
		fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
		
	#------------------------------------------------------------------------------
	#------------------------Create and Submit Job---------------------------------
	#------------------------------------------------------------------------------	

	myJob = mdb.Job(name=JobName, model=ModelName, description='', type=ANALYSIS, 
			atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
			memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE, 
			nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF, 
			contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
			resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=NumDomains, 
			activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=NumDomains)
	myJob.writeInput(consistencyChecking=OFF)
	myJob.submit(consistencyChecking=OFF)



	myJob.waitForCompletion()



	#------------------------------------------------------------------------------
	#---------------------Initialize Global Variables------------------------------
	#------------------------------------------------------------------------------	
	odb = None
	instance = None
	odbName = None

	allSteps = False# consider all steps or single step
	allFrames = False# consider all frames or only final frame

	debug = False
	if os.environ.get('ABQ_COMPOSITEBEAM_DBG'):
		debug = True

	#----------------------------------------------------------
	# User input: odb data
	#----------------------------------------------------------
	odbPath = 'D:/Serena/Repeat_Topology_Studies/Opt18_cutouts_001-09/'+JobName+'.odb'
	stepName = 'rotation'
	stepNames = ['rotation','deployment']
	instanceName = 'smooth_joint_U5_L3-1'

	instanceName = instanceName.upper()

	if debug:
		print 'Instance Name'+instanceName


	#------------------------------------------------------------------------------
	#-----------------------------Failure Calculation------------------------------
	#------------------------------------------------------------------------------	


	def addSectionFrames(step, allFrames):
		"""
		Add FieldOutput to each frame.
		"""


		if allFrames:
			fCount = 0
			for frame in step.frames:
		
				addSectionOdbField(frame)
				if debug:
					fCount+=1
					print 'Frame No.%s is added'%fCount
		else:

			frame = step.frames[-1]
			addSectionOdbField(frame)
			
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	def addSectionOdbField(frame):
		"""
		se (N, M or FI-)is the name of the new fieldOutput
		"""     
		global dataFI
		
		# containters for element labels and stress-resultants
		labels = []
		dataF = []
		dataM = []
		dataFI = []

		# create reference to append fuctions
		labelsAppend = labels.append
		dataFAppend = dataF.append
		dataMAppend = dataM.append
		dataFIAppend = dataFI.append

		lastLabel = -1
		fieldType = VECTOR

		#take field output SF and SM only from 2-ply sets
		force = frame.fieldOutputs['SF']
		moment = frame.fieldOutputs['SM']

		twoply = odb.rootAssembly.instances[instanceName].elementSets['TWO_PLIES']
		force_subset_twoply = force.getSubset(region=twoply)
		moment_subset_twoply = moment.getSubset(region=twoply)

		fValues = force_subset_twoply.values
		mValues = moment_subset_twoply.values
		
		for i in range (len(fValues)):

			fValue = fValues[i]
			mValue = mValues[i]
			
			label = fValue.elementLabel
			if lastLabel != label:
				labelsAppend(label)
			lastLabel = label

			M1, M2, M12 = mValue.dataDouble
			N1, N2, N3, N12, N13, N23 = fValue.dataDouble
		

		
			Nx = 0.5*(N1 + N2) + N12
			Ny = 0.5*(N1 + N2) - N12
			Nxy = 0.5*(N2 - N1)
			
			Mx = 0.5*(M1 + M2) + M12
			My = 0.5*(M1 + M2) - M12
			Mxy = 0.5*(M2 - M1)

			dataValueM = (Mx, My, Mxy)
			dataValueF = (Nx, Ny, Nxy)

	#-----------------------        
	#       In-plane loads
	#-----------------------
			fIndexIP = sJ1*(Nx+Ny) + sK11*( pow(Nx,2)+ pow(Ny,2) ) + sK12*Nx*Ny + sK33*pow(Nxy,2)
		
	#-----------------------
	#       Moments
	#-----------------------
			if (abs(Mx) >= abs(My)):
				M = Mx
			else:
				M = My
			
			fIndexM = sK44*pow(M,2) + sK66*pow(Mxy,2)

	#-----------------------
	#       Interaction
	#-----------------------
			fIndexC = 0
			fIndexCx = 0
			fIndexCy = 0
		
			if (fIndexIP >= 1):
				fIndexC = 1
			
			else:
				# Calculating axial strengths 
				
				# Calculating sFx
				sFx = 0.0001            
				delta = pow((sJ1 + sK12*Ny),2) - 4*sK11*(sJ1*Ny + sK11*pow(Ny,2)+ sK33*pow(Nxy,2)-1)
				
				if delta < 0:
					if debug:
							print 'Value error, deltaX, b^2-4ac = %s'%delta
				elif Ny >= 0:
					sFx = (-(sJ1 + sK12*Ny) + sqrt(delta))/(2*sK11)
				else:
					sFx = (-(sJ1 + sK12*Ny) - sqrt(delta))/(2*sK11)        
				
				fIndexCx = Nx/sFx + abs(M)/sF4

				# Calculating sFy
				sFy = 0.0001
				delta = pow((sJ1 + sK12*Nx),2) - 4*sK11*(sJ1*Nx + sK11*pow(Nx,2)+ sK33*pow(Nxy,2)-1)
				
				if delta < 0:
					if debug:
							print 'Value error deltaY, b^2-4ac = %s'%delta
				elif Nx >= 0:
					sFy = (-(sJ1 + sK12*Nx) + sqrt(delta))/(2*sK11)
				else:
					sFy = (-(sJ1 + sK12*Nx) - sqrt(delta))/(2*sK11)        
				
				fIndexCy = Ny/sFy + abs(M)/sF4
				
			   # selecting maximum index
			
				if (fIndexCx >= fIndexCy):
					fIndexC = fIndexCx 
				else:
					fIndexC = fIndexCy
							
			dataValueFI = (fIndexIP, fIndexM, fIndexC)
		
			dataFAppend(dataValueF)
			dataMAppend(dataValueM)
			dataFIAppend(dataValueFI)

		labels = array(labels, Int32)
		dataF = array(dataF, Float64)
		dataM = array(dataM, Float64)
		dataFI = array(dataFI, Float64)

		if debug:
			print 'len (labels):%s, len(dataF):%s'%(len(labels), len(dataF))

		
		# #----------------------
		# # Add force-resultants
		# #----------------------
		# se = 'Nn'
		# descript = 'Transformed force resultants'
		# componentLabels = ['Nx', 'Ny', 'Nxy']
		# #    componentLabels = ['N1', 'N2', 'N12']    
		# #    componentLabels = [(se + cl) for cl in componentLabels]     
		# fo = frame.FieldOutput(name=se, 
			# description=descript, 
			# type=fieldType, 
			# componentLabels=componentLabels)

		# fo.addData(position=INTEGRATION_POINT, 
			# instance=instance, 
			# labels=labels, 
			# data=dataF)    

		# #----------------------
		# # Add moment-resultants
		# #----------------------
		# se = 'Mn' 
		# #    componentLabels = ['M1', 'M2', 'M12']
		# componentLabels = ['Mx', 'My', 'Mxy']
		# descript = 'Transformed moment resultants'    

		# fo = frame.FieldOutput(name=se, 
			# description=descript, 
			# type=fieldType, 
			# componentLabels=componentLabels)

		# fo.addData(position=INTEGRATION_POINT, 
			# instance=instance, 
			# labels=labels, 
			# data=dataM)  
			
		#----------------------
		# Add failure indices
		#----------------------
		se = 'FI_M'
		componentLabels = ['FI_M-1', 'FI_M-2', 'FI_M-3']
		descript = 'Failure Indices'

		fo = frame.FieldOutput(name=se, 
			description=descript, 
			type=fieldType, 
			componentLabels=componentLabels)

		fo.addData(position=INTEGRATION_POINT, 
			instance=instance, 
			labels=labels, 
			data=dataFI) 
			

	#==========================================================================
	# S T A R T
	#========================================================================== 
	
	if not exists(path+'J_Iter-'+str(iter)+'.lck'):		
		o1 = session.openOdb(name= odbPath, readOnly=0)
		odb = session.odbs[odbPath]
		odbName = odb.name

		assembly = odb.rootAssembly
		instance = assembly.instances[instanceName]
		
		
		# #------------------------------------------------------------------------------
		# #------------------------------Bending Stiffness-------------------------------
		# #------------------------------------------------------------------------------		
		global bendingStiffness
		
		UR1object = odb.steps[stepName].historyRegions['Node RP_CENTER-1.1'].historyOutputs['UR1'].data
		UR1array = array(UR1object, Float64)
		UR1 = UR1array[:,1]

		RM1object = odb.steps[stepName].historyRegions['Node RP_CENTER-1.1'].historyOutputs['RM1'].data
		RM1array = array(RM1object, Float64)
		RM1 = RM1array[:,1]

		bendingStiffness = RM1[1]/UR1[1]

		
		#================================
		#  Add new results 
		#================================    
		if allSteps:
			for sName in stepNames:
				step = odb.steps[sName]
				addSectionFrames(step, allFrames)    
		else:
			step = odb.steps[stepName]    
			addSectionFrames(step, allFrames)    

		odb.save()
		session.odbs[odbName].close()
		print('Saved to ODB')
		
		os.remove('D:/Serena/Repeat_Topology_Studies/Opt18_cutouts_001-09/'+JobName+'.dat')
		os.remove('D:/Serena/Repeat_Topology_Studies/Opt18_cutouts_001-09/'+JobName+'.ipm')
		os.remove('D:/Serena/Repeat_Topology_Studies/Opt18_cutouts_001-09/'+JobName+'.com')
		os.remove('D:/Serena/Repeat_Topology_Studies/Opt18_cutouts_001-09/'+JobName+'.msg')
		os.remove('D:/Serena/Repeat_Topology_Studies/Opt18_cutouts_001-09/'+JobName+'.prt')
		os.remove('D:/Serena/Repeat_Topology_Studies/Opt18_cutouts_001-09/'+JobName+'.sim')
		os.remove('D:/Serena/Repeat_Topology_Studies/Opt18_cutouts_001-09/'+JobName+'.sta')
		#os.remove('D:/Serena/Repeat_Topology_Studies/Opt18_cutouts_001-09/'+JobName+'.log')
	
def func4d(Design_Variables):
	
	global f
	global iter
	iter = iter + 1
	ModelName = 'Iter-'+str(iter)
	JobName = 'J_Iter-'+str(iter)
	
	#------------------------------------------------------------------------------------------------------
	# Load cutouts from Matlab
	global path
	path = 'D:/Serena/Repeat_Topology_Studies/Opt18_cutouts_001-09/'

	cut = 1
	filename = 'cutout-'+str(cut)+'-iter'+str(iter)
	xtemp = []
	ytemp = []

	# Wait for Matlab to generate cutouts files
	while not exists(path+filename):
		time.sleep(5)
	
	# Read cutouts files and generate points (xtemp, ytemp)
	while exists(path+filename):
	 
		with open(filename) as data:
			substrings = data.read().split()
			values = [map(float, substring.split(',')) for substring in substrings]
			length = len(values)/2
			xtemp.append(values[0:length])
			ytemp.append(values[length:])
		
		#os.remove(path+filename)
		cut = cut + 1
		filename = 'cutout-'+str(cut)+'-iter'+str(iter)
		
	#-----------------------------------------------------------------------------------------------------=
	
	# Evaluate objective function FI1/FI1_initial - K/K_initial
	runFullAnalysis(xtemp,ytemp,ModelName,JobName)	
	
	if not exists(path+JobName+'.lck'):
		FI1 = dataFI[:,0]
		FI1newnorm = linalg.norm(FI1)
		FI1max = np.amax(FI1)
		
		FI1max_rep.append(FI1max)
		FI1_initial = FI1max_rep[0]
		B_S_rep.append(bendingStiffness)
		K_initial = B_S_rep[0]
		
		elem_num_1 = 0
		elem_num_1p1 = 0
		elem_num_1p3 = 0
		elem_num_1p5 = 0
		elem_num_2 = 0
		elem_num_3 = 0
		elem_num_5 = 0
		elem_num_10 = 0
		elem_num_15 = 0
		elem_num_20 = 0
		elem_num_25 = 0
		for ii in FI1:
			if 1.0 <= ii:
				elem_num_1 = elem_num_1 + 1
			if 1.1 <= ii:
				elem_num_1p1 = elem_num_1p1 + 1
			if 1.3 <= ii:
				elem_num_1p3 = elem_num_1p3 + 1
			if 1.5 <= ii:
				elem_num_1p5 = elem_num_1p5 + 1
			if 2.0 <= ii:
				elem_num_2 = elem_num_2 + 1
			if 3.0 <= ii:
				elem_num_3 = elem_num_3 + 1
			if 5.0 <= ii:
				elem_num_5 = elem_num_5 + 1
			if 10.0 <= ii:
				elem_num_10 = elem_num_10 + 1
			if 15.0 <= ii:
				elem_num_15 = elem_num_15 + 1
			if 20.0 <= ii:
				elem_num_20 = elem_num_20 + 1
			if 25.0 <= ii:
				elem_num_25 = elem_num_25 + 1
		
		area_elem = pow(mesh_size,2)
		area = elem_num_1*area_elem
		
		area_rep.append(area)
		areaInitial = area_rep[0]
				
		f = (FI1max/FI1_initial) - (bendingStiffness/K_initial)
		f_1 = (FI1max/FI1_initial)
		f_2 = - (bendingStiffness/K_initial)
		f_3 = area
		
		# write data to file
		FUN = open("dataElements_1.txt","a")
		FUN.write("%r %s %r %s\n" % (elem_num_1, 'elements', area, 'mm^2'))
		FUN.close()
		
		FUN = open("dataElements_1p1.txt","a")
		FUN.write("%r %s %r %s\n" % (elem_num_1p1, 'elements', (elem_num_1p1*area_elem), 'mm^2'))
		FUN.close()
		
		FUN = open("dataElements_1p3.txt","a")
		FUN.write("%r %s %r %s\n" % (elem_num_1p3, 'elements', (elem_num_1p3*area_elem), 'mm^2'))
		FUN.close()
		
		FUN = open("dataElements_1p5.txt","a")
		FUN.write("%r %s %r %s\n" % (elem_num_1p5, 'elements', (elem_num_1p5*area_elem), 'mm^2'))
		FUN.close()
		
		FUN = open("dataElements_2.txt","a")
		FUN.write("%r %s %r %s\n" % (elem_num_2, 'elements', (elem_num_2*area_elem), 'mm^2'))
		FUN.close()
		
		FUN = open("dataElements_3.txt","a")
		FUN.write("%r %s %r %s\n" % (elem_num_3, 'elements', (elem_num_3*area_elem), 'mm^2'))
		FUN.close()
				
		FUN = open("dataElements_5.txt","a")
		FUN.write("%r %s %r %s\n" % (elem_num_5, 'elements', (elem_num_5*area_elem), 'mm^2'))
		FUN.close()	

		FUN = open("dataElements_10.txt","a")
		FUN.write("%r %s %r %s\n" % (elem_num_10, 'elements', (elem_num_10*area_elem), 'mm^2'))
		FUN.close()	
				
		FUN = open("dataElements_15.txt","a")
		FUN.write("%r %s %r %s\n" % (elem_num_15, 'elements', (elem_num_15*area_elem), 'mm^2'))
		FUN.close()

		FUN = open("dataElements_20.txt","a")
		FUN.write("%r %s %r %s\n" % (elem_num_20, 'elements', (elem_num_20*area_elem), 'mm^2'))
		FUN.close()

		FUN = open("dataElements_25.txt","a")
		FUN.write("%r %s %r %s\n" % (elem_num_25, 'elements', (elem_num_25*area_elem), 'mm^2'))
		FUN.close()
	
		FUN = open("dataCostFunction.txt","a")
		FUN.write("%r\n" % (f))
		FUN.close()
		
		FUN = open("dataCostFunction_1.txt","a")
		FUN.write("%r\n" % (f_1))
		FUN.close()
		
		FUN = open("dataCostFunction_2.txt","a")
		FUN.write("%r\n" % (f_2))
		FUN.close()
		
		FUN = open("dataCostFunction_3.txt","a")
		FUN.write("%r\n" % (f_3))
		FUN.close()
		
		FUN = open("dataStiffness.txt","a")
		FUN.write("%r\n" % (bendingStiffness))
		FUN.close()
		
		FUN = open("dataFINorm.txt","a")
		FUN.write("%r\n" % (FI1newnorm))
		FUN.close()
		
		FUN = open("designVariables.txt","a")
		FUN.write("%r\n" % (Design_Variables[0]+0.1))
		FUN.close()
		
		FUN = open("cutoutPlane.txt","a")
		FUN.write("%r\n" % (Design_Variables))
		FUN.close()
		
		FUN = open("dataFIMax.txt","a")
		FUN.write("%r\n" % (FI1max))
		FUN.close()
		
	else:
		f = 0.0
		
		FUN = open("failedConvergence.txt","a")
		FUN.write("%r %r %r\n" % (iter, f, Design_Variables))
		FUN.close()
		
		FUN = open("designVariables.txt","a")
		FUN.write("%r\n" % (Design_Variables[0]+0.1))
		FUN.close()
	
	return f

# def print_fun(x, f, accepted):
	
	# # write data to file
	# FUN = open("dataCostFunction_BH.txt","a")
	# FUN.write("%r\n" % (f))
	# FUN.close()

	# FUN = open("dataDesignVariables_BH.txt","a")
	# FUN.write("%r\n" % (x))
	# FUN.close()
	
	# FUN = open("dataAcceptedBH.txt","a")
	# FUN.write("%r\n" % (accepted))
	# FUN.close()

# cons = ({'type':'ineq','fun': lambda Design_Variables: Design_Variables[0] - 0.1},			# cutPlane > 0.1 => cutPlane - 0.1 > 0
		# {'type':'ineq','fun': lambda Design_Variables: 1 - Design_Variables[0]})		# cutPlane < 1 => 1 - cutPlane > 0
			
	

global iter,mesh_size,FI1max_rep,area_rep,B_S_rep
iter = 0
mesh_size = 0.125

FI1max_rep = []
area_rep = []
B_S_rep = []

# minimizer_kwargs = {"method": "COBYLA", "constraints": cons}

for cutPlane in np.arange(0.01,0.91,0.1):

	Design_Variables = np.array([cutPlane])

	# ret = basinhopping(func4d, Design_Variables, minimizer_kwargs=minimizer_kwargs, niter=10, callback=print_fun)
	
	ret = func4d(Design_Variables)


