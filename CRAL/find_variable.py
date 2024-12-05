import os
import re

rootdir="/Users/sraghu/RAMCOMP/back/back_cral/CRAL"
varname="nvector"
varsearchdirs=["amr", "pm", "hydro"]
treesearchdirs=["amr", "pm", "hydro"]
verb=False


def leaf_functions(directory, variable_name, verbose) :
    sdir_path = []
    sub_array = []
    filepath_array = []
    inside_subroutine = False
    

    for sdir in varsearchdirs : 
        sdir_path.append ( directory + "/" + sdir)

    for root, dirs, files in os.walk(directory):
        for dirname in dirs:
            dir_path    = os.path.join(root, dirname)
            if dir_path in sdir_path : 
                #print(dir_path)            
                if os.path.exists(dir_path) and os.path.isdir(dir_path):
                    for filename in os.listdir(dir_path):
                        file_path = os.path.join(dir_path, filename)
                        if os.path.isfile(file_path) and filename.endswith('.f90'):
                            #print("---- : ", file_path)   
                            with open(file_path, 'r') as f:
                                lines = f.readlines()
                                inside_subroutine = False

                                for line_number,line in enumerate(lines, start=1) : 
                                    if (inside_subroutine == False) :
                                        if line.strip().lower().startswith('subroutine'):
                                            #print("----------", line.strip().lower())
                                            line = line.replace("subroutine", "").strip()
                                            subname_begin = line.split('(')[0].strip()
                                        elif line.strip().lower().startswith('recursive subroutine') :
                                            line = line.replace("recursive subroutine", "").strip()
                                            subname_begin = line.split('(')[0].strip()
                                             
                                        if varname in line : 
                                            inside_subroutine = True
                                            #if (verbose): 
                                                #print(">>>>>>>>>>>>>>> :", line_number, line)   
                                    if inside_subroutine == True :
                                        if line.strip().lower().startswith('end subroutine'):     
                                            line = line.replace("end subroutine", "").strip()
                                            subname_end = line #.split('(')[0].strip()
                                            if subname_end == subname_begin :
                                                if (verbose): 
                                                    print(line_number, "-------- inside_subroutine :", subname_begin, "in file:", file_path)  
                                                    print("##################################################")
                                                    print("##################################################")
                                                    print("")
                                                sub_array.append(subname_begin)
                                                filepath_array.append(file_path)
                                                inside_subroutine = False
                                              
    return [sub_array, filepath_array]


def parent_callers(funcname,  root_dir, tree_dirs, verbose) :

    sdir_path = []
    caller_array = []
    for sdir in tree_dirs : 
        sdir_path.append ( root_dir + "/" + sdir)

    #print("call Tree for function : ", funcname)

    for root, dirs, files in os.walk(root_dir):
        for dirname in dirs:
            dir_path    = os.path.join(root, dirname)
            if dir_path in sdir_path:
                if os.path.exists(dir_path) and os.path.isdir(dir_path):
                    func_call = "call " + funcname
                    for filename in os.listdir(dir_path):
                        file_path = os.path.join(dir_path, filename)
                        if os.path.isfile(file_path) and filename.endswith('.f90') and filename != "units.f90":
                            #print("---- : ", file_path)   
                            with open(file_path, 'r') as fp:
                                lines = fp.readlines()
                                inside_subroutine = False
                                for line_number,line in enumerate(lines, start=1) : 
                                    if inside_subroutine == False :
                                        if line.strip().lower().startswith('subroutine') : 
                                            line = line.replace("subroutine", "").strip()
                                            subname_begin = line.split('(')[0].strip()
                                        elif line.strip().lower().startswith('recursive subroutine') :
                                            line = line.replace("recursive subroutine", "").strip()
                                            subname_begin = line.split('(')[0].strip()
                            
                                        if func_call in line.strip().lower():
                                            inside_subroutine = True
                                            if(verbose) :
                                                print(func_call, "  ",line_number, "    ", file_path, "      ", subname_begin )

                                    if inside_subroutine == True :
                                        if line.strip().lower().startswith('end subroutine') :
                                            subname_end = line.replace("end subroutine", "").strip()
                                            if (subname_begin == subname_end) :
                                                inside_subroutine = False
                                                caller_array.append(subname_end)
                                        elif line.strip().lower().startswith('end program ramses') :
                                            subname_end = line.replace("end program", "").strip()
                                            inside_subroutine = False
                                            caller_array.append(subname_end)
        
    #print("callers for : ", funcname, " :: ", caller_array)
    return caller_array

    


#leaf_subroutines, filepath_arr = leaf_functions(rootdir, varname, verb)

#nleafs = len(leaf_subroutines)
#nfiles = len(filepath_arr)

#if (nleafs != nfiles):
#    print("something wrong ")

#print("")
#print("leaf functions for the variable :", varname)
#for ifunc in range (0, nleafs) :
#    print(ifunc, "      ", leaf_subroutines[ifunc], "            " ,filepath_arr[ifunc])        
#print("############################")
#print("")


#func_idx = 1
verb = False
#func_name = leaf_subroutines[func_idx]


func_name='flag_formation_sites'
for itr in range(0, 10): 
	caller_funcs = parent_callers(func_name,  rootdir, treesearchdirs, verb)
	print("caller functions for : ", func_name, " :: ", caller_funcs)
	func_idx=0
	func_name = caller_funcs[func_idx]
	if( func_name=='ramses' or len(caller_funcs) ==0  ) :
		print("Reached end of tree")
		break
