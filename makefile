SHELL   = /bin/sh

export computer_os:=   $(shell uname -s)
export computer_arch:= $(shell uname -m)
export computer_name:=$(shell echo $(HOST) | awk '{ sub(/[0-9]+/,""); print }')

WORK_DIR = work
SUB_DIR = subs
TRTSC_DIR = tr_tsc
TSCM_DIR = tsc_m
TSCA_DIR = tsc_a
TSCR_DIR = tsc_r
TV80_DIR = tv80
FIX_DIR = fixup
GLF_DIR = glf23
MMM71_DIR = mmm_7_1
FIXLSC_DIR = fixuplsc
LSC_DIR = lsc
TRDATBUF_DIR = trdatbuf
VERSION = f90lib
EXE = gotsc_$(computer_name)
EXE_LSC = gotlsc_$(computer_name)

NTCCMOD = ${NTCCHOME}/mod
NETCDFMOD = ${NETCDFHOME}/include

ifeq ($(shell uname -n | cut -c-4),cori)
# --- cori start --- #
export NTCCHOME = /global/homes/j/jinchen/project/ntccsvn/files/tshare/LINUX_cori_intel
export NTCCMOD = ${NTCCHOME}/mod
export NTCCINC = ${NTCCHOME}/include
export NTCCLIB = -L${NTCCHOME}/lib
export NETCDFMOD = ${NETCDF_DIR}/include
export NETCDFLIB = ${NETCDF_DIR}/lib

XPLASMA_LIB = -L${NTCCHOME}/lib \
  -lplasma_state -lps_xplasma2 -lplasma_state_kernel \
  -lold_xplasma -lxplasma2 -lxplasma_debug \
  -lgeqdsk_mds -lmds_sub -lmdstransp \
  -lvaxonly -lnscrunch -lfluxav -lr8bloat \
  -lpspline \
  -llsode -llsode_linpack \
  -lcomput -lportlib \
  -lmclib -lsmlib -ltrgraf -lureadsub -ltridiag

#https://www.nersc.gov/nusers/resources/franklin/software/libsci.php
#Libsci includes BLAS, LAPACK, ScaLAPACK, BLACS, and SuperLU routines.
#there is no explicit reference to LibSci in the compiler line.
#It is included in the default user environment,
#and the compiler wrappers perform the linking automatically.
LAPACK_LIB = 
EZCDF=-L${NTCCHOME}/lib -lezcdf #$(NETCDF_FLIB)

LDRTSC = $(NCAR) $(XLIB)
LD = ftn

else
# --- pppl start --- #
TRXPL2PS_LIB = -L${NTCC_ROOT}/lib \
-ltrxplib -ltrread -ltr_getnl -lrp_kernel -lrplot_io -lxdatmgr -lsplitn \
-ltrdatbuf_lib -lechmod_iolib -luflib -lmds_sub -lufhdf -lkey_access \
-L/usr/lib/gcc/x86_64-redhat-linux/4.4.4 -lstdc++

XPLASMA_LIB = -L${NTCC_ROOT}/lib \
  -lplasma_state -lps_xplasma2 -lplasma_state_kernel \
  -lold_xplasma -lxplasma2 -lgeqdsk_mds \
  -lmdstransp -lnscrunch -lfluxav -lr8bloat -lsmlib -lmclib \
  -linterp_sub -lezcdf \
  -lvaxonly -ltridiag -lcomput -lpspline \
  -llsode -llsode_linpack \
  -lportlib \
  -L${MDSPLUS}/lib -lMdsLib \

LDRTSC = \
 -Wl,-rpath -Wl,${NCARG_ROOT} -L${NCARG_ROOT}/lib -lncarg -lncarg_gks -lncarg_c -lncarg_ras -lngmath -lcgm -lz -lsz \
 -Wl,-rpath -Wl,$/usr/lib -L/usr/lib -lfreetype -lX11 \
 -L$(CAIRO_HOME)/lib -lcairo \

ifeq "$(findstring intel,$(F90HOME))" "intel"
NCDIR=$(NETCDF_C_HOME)
NFDIR=$(NETCDF_FORTRAN_HOME)
LAPACK_LIB = -L${F90HOME}/mkl/lib/intel64 \
             -lmkl_core -lmkl_intel_lp64 \
             -lmkl_lapack95_lp64 -lmkl_blas95_lp64 \
             -lmkl_sequential \
             ##-lmkl_ipf -lmkl_em64t -lbmkl_cdft -lguide
EZCDF=-Wl,-rpath -Wl,${PSPLINE_HOME}/lib -L$(PSPLINE_HOME)/lib -lpspline /p/swim/jchen/TSC/INTEL2019/lib/libpspline.a -Wl,-rpath,$(NCDIR)/lib -Wl,-rpath,$(NFDIR)/lib -L$(NCDIR)/lib -L$(NFDIR)/lib -lnetcdf -lnetcdff
#LD = mpif90 -save -xcheck bounds
LD = mpif90 -save
endif

endif

TSC_LIB  = libtrtsc.a libsubs.a libtsc_m.a libglf.a libmmm_7_1.a libtsc_a.a libtsc_r.a libtv80.a libfixupp.a libtrdatbuf.a  
TLSC_LIB = libtrtsc.a libsubs.a libtsc_m.a libglf.a libmmm_7_1.a libtsc_a.a liblsc.a   libtv80.a libfixlsc.a libtrdatbuf.a

LINK_LIB     = $(TSC_LIB) $(TRXPL2PS_LIB) $(XPLASMA_LIB) $(LAPACK_LIB) $(EZCDF) $(LDRTSC) $(HDF5_HOME)/lib/libhdf5.a $(HDF5_HOME)/lib/libhdf5_cpp.a $(HDF5_HOME)/lib/libhdf5_fortran.a
LINK_LIB_LSC = $(TLSC_LIB) $(XPLASMA_LIB) $(LAPACK_LIB) $(EZCDF) $(LDRTSC) $(HDF5_HOME)/lib/libhdf5.a $(HDF5_HOME)/lib/libhdf5_cpp.a $(HDF5_HOME)/lib/libhdf5_fortran.a

tsc: get_svn_version.c
	@echo "Beginning making a new tsc"
	@echo "compile on machine: " ${computer_os} ${computer_arch} ${computer_name}
	@echo "transp lib: " ${NTCCHOME}/lib
	@echo "lapack lib: " ${LAPACK_LIB}
	@echo "LD lib: " ${LINK_LIB}
	@echo "compiler: " ${LD}
	+@ [ -d $@ ] || mkdir -p $(WORK_DIR) 
	@-(cd $(WORK_DIR);rm *.a *.o)
	@cd $(TRTSC_DIR); gmake
	@cd $(TSCM_DIR); gmake
	@cd $(TSCA_DIR); gmake
	@cd $(TSCR_DIR); gmake
	@cd $(SUB_DIR); gmake
	@cd $(GLF_DIR); gmake
	@cd $(MMM71_DIR); gmake
	@cd $(FIX_DIR);  gmake
	@cd $(TV80_DIR); gmake
	@cd $(TRDATBUF_DIR); gmake
	@-(cd $(WORK_DIR);ln -s ../tr_tsc/libtrtsc.a .) 
	@-(cd $(WORK_DIR);ln -s ../tsc_m/libtsc_m.a .) 
	@-(cd $(WORK_DIR);ln -s ../tsc_a/libtsc_a.a .) 
	@-(cd $(WORK_DIR);ln -s ../tsc_r/libtsc_r.a .) 
	@-(cd $(WORK_DIR);ln -s ../fixup/libfixupp.a .)
	@-(cd $(WORK_DIR);ln -s ../tv80/libtv80.a .)
	@-(cd $(WORK_DIR);ln -s ../subs/libsubs.a .); 
	@-(cd $(WORK_DIR);ln -s ../glf23/libglf.a .;)
	@-(cd $(WORK_DIR);ln -s ../mmm_7_1/libmmm_7_1.a .;)
	@-(cd $(WORK_DIR);ln -s ../trdatbuf/libtrdatbuf.a .)
	@-(cd $(WORK_DIR);ar -x libsubs.a gf.o)
	@-(cd $(WORK_DIR);ar -x libtsc_m.a limits.o)
	@-(cd $(WORK_DIR);ar -x libtrtsc.a tsc.o auxheat.o curdrive.o sawtooth.o)
	@-(cd $(WORK_DIR);$(LD) -o $(EXE) -v tsc.o auxheat.o curdrive.o sawtooth.o limits.o gf.o $(LINK_LIB))
	@mv $(WORK_DIR)/$(EXE) $(EXE)
	@echo "Release version - gotsc - is now updated"
	@echo ""

tlsc: 
	@echo "Beginning making tlsc"
	@echo "compile on machine: " ${computer_os} ${computer_arch} ${host_name}
	@echo "transp lib: " ${NTCCHOME}/lib
	@echo "lapack lib: " ${LAPACK_LIB}
	@echo "lD lib: " ${LINK_LIB}
	@echo "compiler: " ${LD}
	+@ [ -d $@ ] || mkdir -p $(WORK_DIR) 
	@-(cd $(WORK_DIR);rm *.a *.o)
	@cd $(TRTSC_DIR); gmake
	@cd $(TSCM_DIR); gmake
	@cd $(TSCA_DIR); gmake
	@cd $(LSC_DIR);  gmake
	@cd $(SUB_DIR); gmake
	@cd $(GLF_DIR); gmake
	@cd $(MMM71_DIR); gmake
	@cd $(FIXLSC_DIR);  gmake
	@cd $(TV80_DIR); gmake
	@cd $(TRDATBUF_DIR); gmake
	@-(cd $(WORK_DIR);ln -s ../tr_tsc/libtrtsc.a .) 
	@-(cd $(WORK_DIR);ln -s ../tsc_m/libtsc_m.a .)
	@-(cd $(WORK_DIR);ln -s ../tsc_a/libtsc_a.a .)
	@-(cd $(WORK_DIR);ln -s ../lsc/liblsc.a .) 
	@-(cd $(WORK_DIR);ln -s ../fixuplsc/libfixlsc.a .)
	@-(cd $(WORK_DIR);ln -s ../tv80/libtv80.a .)
	@-(cd $(WORK_DIR);ln -s ../subs/libsubs.a .)
	@-(cd $(WORK_DIR);ln -s ../glf23/libglf.a .)
	@-(cd $(WORK_DIR);ln -s ../mmm_7_1/libmmm_7_1.a .)
	@-(cd $(WORK_DIR);ln -s ../trdatbuf/libtrdatbuf.a .)
	@-(cd $(WORK_DIR);ar -x libsubs.a gf.o)
	@-(cd $(WORK_DIR);ar -x libtsc_m.a limits.o)
	@-(cd $(WORK_DIR);ar -x libtrtsc.a tsc.o auxheat.o curdrive.o sawtooth.o)
	@-(cd $(WORK_DIR);$(LD) -o $(EXE_LSC) tsc.o auxheat.o curdrive.o sawtooth.o limits.o gf.o $(LINK_LIB_LSC) )
	@mv $(WORK_DIR)/$(EXE_LSC) $(EXE_LSC)
	@echo "Release version - gotlsc - is now updated"
	@echo ""

##
## on every build, record the working copy revision string
##
##
## Then any executable that links in get_svn_version.o will be able
## to call the function get_svn_version() to get a string that
## describes exactly what revision was built.
##

	@echo ""
.PHONY: clean
clean:
	@echo "Beginning clean up directories"
	@rm -rf $(WORK_DIR) 
	@rm -rf $(TRTSC_DIR)/object
	@rm -rf $(TSCM_DIR)/object
	@rm -rf $(TSCA_DIR)/object
	@rm -rf $(TSCR_DIR)/object
	@rm -rf $(SUB_DIR)/object
	@rm -rf $(FIX_DIR)/object
	@rm -rf $(GLF_DIR)/object
	@rm -rf $(MMM71_DIR)/object
	@rm -rf $(TV80_DIR)/object
	@rm -rf $(TRDATBUF_DIR)/object
	@rm -rf $(FIXLSC_DIR)/object
	@rm -rf $(LSC_DIR)/object
	@rm -f $(TRTSC_DIR)/*.a
	@rm -f $(TSCM_DIR)/*.a
	@rm -f $(TSCA_DIR)/*.a
	@rm -f $(TSCR_DIR)/*.a
	@rm -f $(SUB_DIR)/*.a
	@rm -f $(GLF_DIR)/*.a
	@rm -f $(MMM71_DIR)/*.a
	@rm -f $(LSC_DIR)/*.a
	@rm -f $(FIX_DIR)/*.a
	@rm -f $(FIXLSC_DIR)/*.a
	@rm -f $(TV80_DIR)/*.a
	@rm -f $(TRDATBUF_DIR)/*.a
	@echo "working directories cleaned up"
	@echo " "

