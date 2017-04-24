#!/Users/dora/Library/Enthought/Canopy_32bit/User/bin/python

import subprocess as subproc
import sys

# possible choices:   'torus_hd', 'torus_mhd', 'hkdisk_mhd', 'hkdisk_hd'

case_to_do =  self_grav_jns2d

mpi = False               

# mpi = True
#./configure --with-coord=cylindrical --with-problem=cylwindrotb


PATH_BASE='/Users/dora/WORK/ECLIPSE_SPACE/'
PATH_BASE = '/Users/dora/WORK/ECLIPSE_SPACE/AthenaWind'

import socket
name=socket.gethostname()

#print(len(sys.argv)); exit();

if name == 'atorus':
    PATH_BASE = '/local/data/atorus1/dora/PROJECTS'

if name == 'Antons-MacBook-Pro.local':
    PATH_BASE = "/Users/dora/WORK/ECLIPSE_SPACE"

PATH = PATH_BASE + '/AthenaWind/'

print("current Path: ",  PATH)

if  case_to_do==self_grav_jns2d:
    problemToConfig = '--with-problem=jeans'
    sgrav = '--with-gravity=fft'
    fft = '--enable-fft'
    
if case_to_do == 'hkdisk_mhd':    
    problemToConfig = '--with-problem=hkdisk'
    methodGasOrMHD =  '--with-gas=mhd'
    inputFile = '../tst/cylindrical/athinput.hkdisk-3D'
    
if case_to_do == 'hkdisk_hd':
    problemToConfig = '--with-problem=hkdisk'
    methodGasOrMHD =  '--with-gas=hydro'
    inputFile = '../tst/cylindrical/athinput.hkdisk-3D'

if case_to_do == 'torus_hd':    
    problemToConfig = '--with-problem=torus9'
    methodGasOrMHD =  '--with-gas=hydro'    
    inputFileList = ['../tst/cylindrical/athinput.torus9_hydro_2D', '../tst/cylindrical/athinput.torus9_hydro_2D_2']
    inputFile = '../tst/cylindrical/athinput.torus9_hydro_2D'
   
    
if case_to_do == 'torus_mhd':    
    problemToConfig = '--with-problem=torus9'
    methodGasOrMHD =  '--with-gas=mhd'
    METHOD = '--with-flux=hlld'   
    # METHOD = '--with-flux=hlle'
    ORDER = '--with-order=2p'

    if mpi:
        inputFile = '../tst/cylindrical/athinput.torus9_MPI'     
    else:
        inputFile = '../tst/cylindrical/athinput.torus9_hydro_2D'


subproc.check_call(['rm', '-f', './bin/*.bin'])

 
if case_to_do == 'torus_hd':
    METHOD = '--with-flux=roe'
else:    
   Integrator = '--with-integrator=vl'

compLev = '0123'


if '0' in compLev:
    subproc.check_call(['make', 'clean'])

    if mpi:
        subproc.check_call([PATH+'./configure', \
        '--with-coord=cylindrical', Integrator, \
        '--enable-fofc', methodGasOrMHD, METHOD, \
        ORDER, '--enable-mpi', problemToConfig])
    else:
        subproc.check_call([PATH+'./configure', \
        '--with-coord=cylindrical', Integrator, \
        '--enable-fofc', methodGasOrMHD, METHOD, \
        ORDER, problemToConfig])

    # subproc.check_call([PATH+'./configure', '--with-coord=cylindrical',Integrator, '--enable-fofc', methodGasOrMHD, METHOD, ORDER, problemToConfig])

#    subproc.check_call([PATH+'./configure', '--with-coord=cylindrical', methodGasOrMHD, METHOD, ORDER, problemToConfig])
    

if '1' in compLev:
    subproc.check_call(['make', 'clean'])

if '2' in compLev:
    if mpi:
        subproc.check_call(['make', 'all', 'MACHINE=macosxmpi'])        
    else:
        subproc.check_call(['make', 'all', 'MACHINE=macosx'])        

if '3' in compLev:
   
    if mpi:
        subproc.check_call(['/Users/dora/Applications/LIBS/openmpi/bin/mpirun', '-np', '4', './athena', '-i', \
                        inputFile],  cwd = './bin' )
    else:
        subproc.check_call(['./athena', '-i', inputFile],  cwd = './bin' )

        
#make all MACHINE=macosxmpi

#/Users/dora/Applications/LIBS/openmpi/bin/mpirun -np 4 ./athena -i
                        #../tst/cylindrical/athinput.torus9_hydro_2D
