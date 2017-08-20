#
# main targets: rebuild machine-generated coefficient files
#

default	: 1d 2d 3d
dirs	: 1d.dirs 2d.dirs 3d.dirs

.PHONY	: 1d
1d	:
	-rm -f 1d.cube.size*/*.c 1d.cube.size*/*.h
	cat ../util.maple ../interpolate.maple \
            cube_posns.maple 1d.maple \
           | maple 2>&1 >1d.log

.PHONY	: 2d
2d	:
	-rm -f 2d.cube.size*/*.c 2d.cube.size*/*.h
	cat ../util.maple ../interpolate.maple \
            cube_posns.maple 2d.maple \
           | maple 2>&2 >2d.log

.PHONY	: 3d
3d	:
	-rm -f 3d.cube.size*/*.c 3d.cube.size*/*.h
	cat ../util.maple ../interpolate.maple \
            cube_posns.maple 3d.maple \
           | maple 2>&2 >3d.log

################################################################################

#
# targets to initially create directories for machine-generated coefficients
# these only need to be run if you're setting up a new interpolator,
# i.e. if you don't already have the directories in cvs
# (you still have to add them to cvs by hand)
#

.PHONY	: 1d.dirs
1d.dirs:
	-rm -rf 1d.cube.size*
	mkdir 1d.cube.size2
	mkdir 1d.cube.size3
	mkdir 1d.cube.size4
	mkdir 1d.cube.size5
	mkdir 1d.cube.size6
	mkdir 1d.cube.size7

.PHONY	: 2d.dirs
2d.dirs:
	-rm -rf 2d.cube.size*
	mkdir 2d.cube.size2
	mkdir 2d.cube.size3
	mkdir 2d.cube.size4
	mkdir 2d.cube.size5
	mkdir 2d.cube.size6

.PHONY	: 3d.dirs
3d.dirs:
	-rm -rf 3d.cube.size*
	mkdir 3d.cube.size2
	mkdir 3d.cube.size3
	mkdir 3d.cube.size4
	mkdir 3d.cube.size5
	mkdir 3d.cube.size6
