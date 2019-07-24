"""
**************************************************************
			MESH REFINEMENT NURBS OPTIMIZATION TEST (2D)
**************************************************************

This python script automates running the NURBS optimization test case with different 
parameters. 
The tests corresponding to mesh dimensions in MESH_DIM_LIST are run.
The order of the polynomials in the reference domain is set by P_LIST. P_LIST[i] has 
polynomial order corresponding to MESH_DIM_LIST[i]. 
NURBS_ORDER_P is the order of the NURBS basis function in the direction tangent to the wing surface. 
NURBS_NUM_CTRL_PTS_XI is the number of control points parametrizing the airfoil spline. The same
number of control points is used for the farfield circle. See NURBS_Parametric_Domain_Generator for details
on the NURBS parameters. 

The results are saved in directory specified by ARCHIVED_OUTPUT_PATH. 
These include the optimization results, the input files used to generate them (in the Parameters subdirectory)
as well as the console output saved in log.txt. 
"""

import os
import sys
import subprocess
import shlex
import re
import shutil 
import fnmatch

PROJECT_SRC_DIR="/home/vassili/Desktop/DPGSolver/DPGSolver/"
NURBS_GENERATION_REL_DIR="input/meshes/NURBS_Patch_Generator/"
sys.path.append(PROJECT_SRC_DIR+NURBS_GENERATION_REL_DIR)

#editor complains here because it cant process the sys path append
from NURBS_Parametric_Domain_Generator import output_file as nurbs_output_file
from Airfoil_Patch import get_patch_information

#TWO MAIN THINGS TO CONTROL. 
# MESH DIMENSIONS FOR WHICH YOU WANT TO RUN THE TESTS
# AND LIST OF ORDERS AT WHICH SHOULD BE COMPUTED
MESH_DIM_LIST=[[38,27], #1000 dof
				[19,14],
				[14,10], 
				[53,38], #2000 dof
				[28,19],
				[18,13],
				[65,46], #3000 dof
				[34,23],
				[23,16], 
#				[75,54], #4000
#				[38,27],
#				[26,18]
]
#MESH_DIM_LIST is looped over. Length of this fixes the amount of tests
#that are run
#MESH_DIM_LIST=MESH_DIM_LIST#+MESH_DIM_LIST
MESH_DIM_LIST=[[53,38], #2000 dof
				[28,19],
				[18,13]]*3
P_LIST=[1,2,3]*4
#P_LIST=P_LIST+P_LIST
#MESH_DIM_LIST=[[10,10]]

#NURBS GEOMETRY PARAMETERS
NURBS_ORDER_P=3
#NURBS_NUM_CTRL_PTS_XI_LIST=[11]*9+[15]*9
NURBS_NUM_CTRL_PTS_XI_LIST=[13]*3+[17]*3+[19]*3

TEST_NAME_FORMAT="TEST_Python_Automated_Euler_NURBSAirfoil_TargetCLTestsReference_ParametricQUAD2D_mesh_%dby%d__ml%d__p%d"
BUILD_DIR=PROJECT_SRC_DIR+"build_debug_2D/"
LOG_OUTPUT_PATH=BUILD_DIR+"output/optimization/euler/steady/NURBS_Airfoil/Constrained_TargetCL/TestsReference/"
ARCHIVED_OUTPUT_PATH=PROJECT_SRC_DIR+"Archived_Output/Vassili/"

CTRL_SUBFOLDER="euler/NURBS_Airfoil/"
CONTROL_FILE_FOLDER=PROJECT_SRC_DIR+"input/testing/control_files/"+CTRL_SUBFOLDER
CONTROL_FILE_TEMPLATE="TEST_Euler_NURBSAirfoil_Python_Template.ctrl"
GEO_FILE_FOLDER=PROJECT_SRC_DIR+"input/meshes/n-cube/"
GEO_FILE_TEMPLATE="2d_stretched_python_template.geo"
GEO_FORMAT="2d_stretched_temp.geo"

NURBS_GEO_FILE_FOLDER=PROJECT_SRC_DIR+"input/input_files/euler/steady/NURBS_Airfoil/"
NURBS_EXTENSION_FORMAT="airfoil_P%d_NumPtsXi%d_Q1_NumPtsEta2"
NURBS_GEO_FILE_FORMAT="geometry_parameters_" + NURBS_EXTENSION_FORMAT +".geo"

ML=1 
#

def find(pattern, path):
	"""Find files that match pattern in the subdirectories of path
	https://stackoverflow.com/questions/1724693/find-a-file-in-python"""
	result = []
	for root, dirs, files in os.walk(path):
		for name in files:
			if fnmatch.fnmatch(name, pattern):
				result.append(os.path.join(root, name))
	return result

def copy_files(src_dir, dest_dir, fname_list=[], f_fulldir=0):
	"""Copy files from source directory to destination directory. 
	By default, copy all of them unless a file name list is specified.
	By default, files specified as file names. If full directory returned, 
		set f_fulldir to 1. """
	files_to_cp=[]
	if fname_list: #if list not empty
		files_to_cp=fname_list
	else:
		assert(os.path.exists(src_dir))
		files_to_cp = os.listdir(src_dir)
	for file_name in files_to_cp:
		if f_fulldir==0:
			full_file_name = os.path.join(src_dir, file_name)
		else:
			full_file_name=file_name
		assert(os.path.isfile(full_file_name))
		os.makedirs(os.path.dirname(dest_dir), exist_ok=True)
		shutil.copy(full_file_name, dest_dir)

def run_test(test_name, test_identifier):
	#TODO: GET LOG SAVING WORKING... 
	"""Runs command for test
	Command format obtained by ctest -R testname -V"""
	log_fpath=ARCHIVED_OUTPUT_PATH+test_identifier+"/log.txt"
	test_cmd_str=BUILD_DIR+"bin/test_integration_optimization " \
		+ CTRL_SUBFOLDER+ test_name+" petsc_options_gmres_default"#+" > "+\
		#	
	test_command=shlex.split(test_cmd_str)
	f=open(log_fpath, "w+")
	subprocess.call(test_command, cwd=BUILD_DIR+"bin/", stdout=f)
	f.close
	return "Test Command: "+test_cmd_str

def copy_files_to_archive(input_fname_list, test_identifier):
	"""Copy files from build output directory to the archived file directory
	Copies both the test output files and the control/geo files used for the test"""
	
	copy_files(src_dir=LOG_OUTPUT_PATH, dest_dir=ARCHIVED_OUTPUT_PATH+test_identifier+"/")
	copy_files("", dest_dir=ARCHIVED_OUTPUT_PATH+test_identifier+"/Parameters/", 
		fname_list=input_fname_list, f_fulldir=1)


def create_input_files(mesh_dims, ml,P, nurbs_p, nurbs_num_ctrl_pts_xi):
	"""Creates the geometry file and the control file that are used by the C code."""
	geo_fname=GEO_FILE_FOLDER+GEO_FORMAT
	ctrl_fname=CONTROL_FILE_FOLDER+TEST_NAME_FORMAT % (mesh_dims[0], mesh_dims[1], ML, P)+'.ctrl'
	nurbs_geo_fname=NURBS_GEO_FILE_FOLDER+NURBS_GEO_FILE_FORMAT % (nurbs_p, nurbs_num_ctrl_pts_xi)
	nurbs_geo_extension=NURBS_EXTENSION_FORMAT % (nurbs_p, nurbs_num_ctrl_pts_xi)

	geo_fname_nopath=GEO_FORMAT
	#CREATE NURBS GEOMETRY FILE
	patch_information=get_patch_information(NURBS_ORDER_P, nurbs_num_ctrl_pts_xi)
	nurbs_output_file(patch_information,nurbs_geo_fname)

	# CREATE GEOMETRY FILE
	f = open(GEO_FILE_FOLDER+GEO_FILE_TEMPLATE,'r')
	filedata = f.read()
	f.close()
	#need to delete after
	newdata = filedata.replace("PYTH_TEMP_NUM_X_ELEMENTS",str(mesh_dims[0]))
	newdata = newdata.replace("PYTH_TEMP_NUM_Y_ELEMENTS",str(mesh_dims[1]))
	f = open(geo_fname,'w')
	f.write(newdata)
	f.close()

	#CONTROL FILE
	f = open(CONTROL_FILE_FOLDER+CONTROL_FILE_TEMPLATE,'r')
	filedata = f.read()
	f.close()

	newdata = filedata.replace("PYTH_TEMP_ML",str(ml))
	newdata = newdata.replace("PYTH_TEMP_P",str(P))
	newdata = newdata.replace("PYTH_TEMP_GEO_FILE",geo_fname_nopath)
	newdata = newdata.replace("PYTH_TEMP_GEO_PARAMETERS_EXTENSION",nurbs_geo_extension)
	f = open(ctrl_fname,'w')
	f.write(newdata)
	f.close()


	return geo_fname, ctrl_fname, nurbs_geo_fname

def append_to_log(contents, fpath):
	f=open(fpath, "a+")
	f.write(contents)
	f.close

def delete_files_from_dir(d):
	"""Remove all files from directory specified by d
	https://stackoverflow.com/questions/1039711/remove-all-files-in-a-directory?lq=1"""
	if os.path.exists(d):
		filesToRemove = [os.path.join(d,f) for f in os.listdir(d)]
		for f in filesToRemove:
				os.remove(f) 

if __name__ == "__main__":
	""" Run ctest with TEST_NAME. Then copy the output files to appropriate directory in Archived_Output """
	#run_tests(TEST_NAME_LIST)
	for idx, mesh_dims in enumerate(MESH_DIM_LIST):
		delete_files_from_dir(LOG_OUTPUT_PATH)
		
		P=P_LIST[idx]
		nurbs_p = NURBS_ORDER_P
		num_ctrl_pts_xi=NURBS_NUM_CTRL_PTS_XI_LIST[idx]
		
		test_name=TEST_NAME_FORMAT % (mesh_dims[0], mesh_dims[1], ML, P)
		#test_identifier="M%dx%d_P%d_NURBS_P%d_NCTRL%d" % (mesh_dims[0], mesh_dims[1], P, nurbs_p, num_ctrl_pts_xi)
		test_identifier="NCTRL%d/M%dx%d_P%d_NURBS_P%d" % (num_ctrl_pts_xi, mesh_dims[0], mesh_dims[1], P, nurbs_p)

		geo_fname, ctrl_fname, nurbs_geo_fname=create_input_files(mesh_dims, ML, P, nurbs_p, num_ctrl_pts_xi)
		#create output directory
		os.makedirs(os.path.dirname(ARCHIVED_OUTPUT_PATH+test_identifier+"/"), exist_ok=True)
		subprocess.call(shlex.split("make all"), cwd=BUILD_DIR)
		subprocess.call(shlex.split("make meshes"), cwd=BUILD_DIR)
		test_command_to_log=run_test(test_name, test_identifier)
		copy_files_to_archive([geo_fname, ctrl_fname, nurbs_geo_fname], test_identifier)
		print(test_command_to_log + "\n Done. ")
		#append_to_log(test_identifier+" : "+
		#		test_command_to_log, ARCHIVED_OUTPUT_PATH+test_identifier+"/log.txt")
		os.remove(geo_fname)
		os.remove(ctrl_fname)
		os.remove(nurbs_geo_fname)





	
	
	
	
	


	

	

	#print(control_files)
	#copy_files(src_dir=LOG_OUTPUT_PATH, dest_dir=ARCHIVED_OUTPUT_PATH+TEST_NAME+"/")
	

			


	#python_call = Python_call(project_src_dir,mesh_name_full)

	#python_call.set_input(project_src_dir,mesh_name)
	#python_call.set_output(mesh_name_full)
	#python_call.call_function(mesh_name_full)