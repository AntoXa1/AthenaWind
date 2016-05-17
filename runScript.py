#!/Applications/Canopy64.app/appdata/canopy-1.5.1.2730.macosx-x86_64/Canopy.app/Contents/bin/python
import subprocess as subproc
import sys

# caseToDo = 'torus_hd'
caseToDo = 'torus_mhd'
# caseToDo = 'hkdisk_mhd'
# caseToDo = 'hkdisk_hd'

#./configure --with-coord=cylindrical --with-problem=cylwindrotb
PATH = "/Users/dora/WORK/ECLIPSE_SPACE/AthenaWind"


if caseToDo == '1':    
    problemToConfig = '--with-problem=cylwindrot'
    methodGasOrMHD =  '--with-gas=hydro'
    inputFile = '../tst/cylindrical/athinput.cylwindrot-3D'
    

if caseToDo == '2':    
    problemToConfig = '--with-problem=cylwindrotb'
    methodGasOrMHD =  '--with-gas=mhd'
    inputFile = '../tst/cylindrical/athinput.cylwindrotb-2D'

if caseToDo == '3':    
    problemToConfig = '--with-problem=torus9'

    methodGasOrMHD =  '--with-gas=hydro'
    
    inputFile = '../tst/cylindrical/athinput.torus9_hydro_2D'

if caseToDo == 'hkdisk_mhd':    
    problemToConfig = '--with-problem=hkdisk'
    methodGasOrMHD =  '--with-gas=mhd'
    inputFile = '../tst/cylindrical/athinput.hkdisk-3D'

if caseToDo == 'hkdisk_hd':
    problemToConfig = '--with-problem=hkdisk'
    methodGasOrMHD =  '--with-gas=hydro'
    inputFile = '../tst/cylindrical/athinput.hkdisk-3D'


if caseToDo == 'torus_hd':    
    problemToConfig = '--with-problem=torus9'
    methodGasOrMHD =  '--with-gas=hydro'    
    inputFile = '../tst/cylindrical/athinput.torus9_hydro_2D'

if caseToDo == 'torus_mhd':    
    problemToConfig = '--with-problem=torus9'

    methodGasOrMHD =  '--with-gas=mhd'
    
    inputFile = '../tst/cylindrical/athinput.torus9_hydro_2D'

subproc.check_call(['rm', '-f', './bin/*.bin'])

subproc.check_call(['/Users/dora/WORK/ECLIPSE_SPACE/AthenaWind/./configure', '--with-coord=cylindrical',methodGasOrMHD, problemToConfig])
subproc.check_call(['make', 'clean'])
subproc.check_call(['make', 'all', 'MACHINE=macosx'])
subproc.check_call(['./athena', '-i', inputFile],  cwd = './bin' )
