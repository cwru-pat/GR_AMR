# $Header: /usr/local/svn/cvs-repositories/numrelcvs/AEIThorns/AEILocalInterp/src/Lagrange-maximum-degree/makefile,v 1.1.1.1 2003-07-06 11:16:15 jthorn Exp $

.PHONY	: no-default-target
no-default-target:
	@echo 'there is (deliberately) no default target in this makefile,'
	@echo 'because (re)generating all the coefficients is quite slow :('
	@echo '==> see the makefile to find out about explicit targets'

#
# main targets: rebuild machine-generated coefficient files
#

.PHONY	: 1d
1d	:
	-rm -f 1d.coeffs/*/*.c
	cat ../util.maple ../interpolate.maple \
            ../common/cube_posns.maple fns.maple 1d.maple \
           | maple 2>&1 >1d.log

# with Maple 7, this takes about 15 seconds on a 1.7GHz P4
.PHONY	: 2d
2d	:
	-rm -f 2d.coeffs/*/*.c
	cat ../util.maple ../interpolate.maple \
            ../common/cube_posns.maple fns.maple 2d.maple \
           | maple 2>&1 >2d.log

# with Maple 7, this takes about 8 to 9 minutes on a 1.7GHz P4
.PHONY	: 3d
3d	:
	-rm -f 3d.coeffs/*/*.c
	cat ../util.maple ../interpolate.maple \
            ../common/cube_posns.maple fns.maple 3d.maple \
           | maple 2>&1 >3d.log

################################################################################

#
# targets to initially create directories for machine-generated coefficients
# these only need to be run if you're setting up a new interpolator,
# i.e. if you don't already have the directories in cvs
# (you still have to add them to cvs by hand)
#

.PHONY	: 1d.dirs
1d.dirs:
	-rm -rf 1d.coeffs
	mkdir 1d.coeffs
	cd 1d.coeffs && mkdir 1d.cube.order1.smooth0
	cd 1d.coeffs && mkdir 1d.cube.order2.smooth0
	cd 1d.coeffs && mkdir 1d.cube.order3.smooth0
	cd 1d.coeffs && mkdir 1d.cube.order4.smooth0
	cd 1d.coeffs && mkdir 1d.cube.order5.smooth0
	cd 1d.coeffs && mkdir 1d.cube.order6.smooth0

.PHONY	: 2d.dirs
2d.dirs:
	-rm -rf 2d.coeffs
	mkdir 2d.coeffs
	cd 2d.coeffs && mkdir 2d.cube.order1.smooth0
	cd 2d.coeffs && mkdir 2d.cube.order2.smooth0
	cd 2d.coeffs && mkdir 2d.cube.order3.smooth0
	cd 2d.coeffs && mkdir 2d.cube.order4.smooth0

.PHONY	: 3d.dirs
3d.dirs:
	-rm -rf 3d.coeffs
	mkdir 3d.coeffs
	cd 3d.coeffs && mkdir 3d.cube.order1.smooth0
	cd 3d.coeffs && mkdir 3d.cube.order2.smooth0
	cd 3d.coeffs && mkdir 3d.cube.order3.smooth0
	cd 3d.coeffs && mkdir 3d.cube.order4.smooth0
