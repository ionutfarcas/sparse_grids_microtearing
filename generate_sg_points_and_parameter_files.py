###!/usr/bin/env python3
### -*- coding: utf-8 -*-
from sg_lib.grid.grid import *
from sg_lib.algebraic.multiindex import *
import numpy as np
import os
import shutil

# shat
left_shat  	= 3.8
right_shat  = 21.7

# coll
left_coll     = 6.5e-3
right_coll    = 8.7e-2

# amhd
left_amhd     = 0.8
right_amhd    = 20.0

# q0
left_q0     = 3.5
right_q0    = 9.0

# beta
left_beta     = 0.0002
right_beta    = 0.006

# temp i
left_temp_i    = 0.9
right_temp_i   = 3.5

# omn
left_omn    = 5.0
right_omn   = 240.0

# omt i
left_omt_i    = 14.0
right_omt_i   = 72.0

# omt e
left_omt_e    = 14.0
right_omt_e   = 360.0


if __name__ == '__main__':

    dim = 9 # number of parameters in the scan
    
    ############### only these three parameters should be modified ############
    level 	 = 4 # grid level; this tells us the maximum monomial degree in each direction (which is level - 1)
    out_dir  = lambda n: '/pscratch/sd/d/drhatch/miller_tests/' # in case you want to change the output directory to save each ky scan
    dir_name = lambda n: 'ky_scan_sg_point_' + str(n + 1) + '/' # folder names that contain the parameter files
    #########################################################################

    level_to_nodes  = 1
    weights 		= [lambda x: 1. for d in range(dim)]
    left_bounds    	= np.array([0 for d in range(dim)])
    right_bounds   	= np.array([1 for d in range(dim)])


    Grid_obj = Grid(dim, level, level_to_nodes, left_bounds, right_bounds, weights)	
    Multiindex_obj = Multiindex(dim)
    
    multiindex_set = Multiindex_obj.get_std_total_degree_mindex(level)

    sg_points 	= Grid_obj.get_std_sg_surplus_points(multiindex_set)
    n_sg 		= sg_points.shape[0] 

    print('total number of grid points for dim =', dim, 'and level =', level, ' is n =', n_sg)
   
    # map sg points from [0, 1] to the respective bounds for all d parameters
    map_rv = lambda left, right, x: left + (right - left)*x
   
    sg_shat     = map_rv(left_shat, right_shat, sg_points[:, 0])
    sg_coll     = map_rv(left_coll, right_coll, sg_points[:, 1])
    sg_amhd     = map_rv(left_amhd, right_amhd, sg_points[:, 2])
    sg_q0       = map_rv(left_q0, right_q0, sg_points[:, 3])
    sg_beta     = map_rv(left_beta, right_beta, sg_points[:, 4])
    sg_temp_i   = map_rv(left_temp_i, right_temp_i, sg_points[:, 5])
    sg_omn      = map_rv(left_omn, right_omn, sg_points[:, 6])
    sg_omt_i    = map_rv(left_omt_i, right_omt_i, sg_points[:, 7])
    sg_omt_e    = map_rv(left_omt_e, right_omt_e, sg_points[:, 8])

    for f in range(n_sg):
        if not os.path.isdir(dir_name(f)):
            os.mkdir(dir_name(f))
            
    for n in range(n_sg):
        source      = 'parameters'
        destination = dir_name(n) + source
        
        shutil.copy2(source, destination)
        
        command = './modify_params_9D ' + destination + ' ' + \
					out_dir(n) + ' ' + \
                    str(sg_shat[n]) + ' ' + \
					str(sg_coll[n]) + ' ' + \
					str(sg_amhd[n]) + ' ' + \
					str(sg_q0[n]) + ' ' + \
					str(sg_beta[n]) + ' ' + \
					str(sg_temp_i[n]) + ' ' + \
					str(sg_omn[n]) + ' ' + \
                    str(sg_omt_i[n]) + ' ' + \
                    str(sg_omt_e[n])

        print(command)
        os.system(command)