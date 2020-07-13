#!/bin/sh

#include "config.h"
MAKEIN = ./make.in
include $(MAKEIN)

# --- Makefile for Quick Program ---
#  				- v 3.00 2019/03/30 Madu Manathunga
#				- v 2.00 2010/10/25 Yipu Miao
#				- v 1.18 2009/09/16 John Faver Exp $ 
#				- Makefile created by mkmf.pl $Id:
#	--------
#	 INDEX
#	--------
#	A. Compiler Setting			! Intel Fortran 9.0+ or GNU Fortran is recommended for single CPU Version
#	B. Make folders
#		! mpif90 is recommended for MPI Multi-CPU Version
#	C. Make Object Files		! Source files --> Object files
#	D. Make Executed files		! Object files --> Executed files
#	E. Self-defined Option		! Make option

#************************************************************************
#                  A. Compiler Settings
# 
#   FC specifies f90 compiler
#   FFLAGS are compliation options
#   LFLAGS are linking flags
#
#************************************************************************

cudaobj    =   $(objfolder)/gpu_write_info.o $(objfolder)/gpu.o $(objfolder)/gpu_type.o \
	$(objfolder)/gpu_get2e.o
cudaxcobj = $(objfolder)/gpu_getxc.o

cudalibxcobj=$(objfolder)/gga_c_am05.o $(objfolder)/gga_c_bcgp.o $(objfolder)/gga_c_bmk.o $(objfolder)/gga_c_cs1.o \
	$(objfolder)/gga_c_ft97.o $(objfolder)/gga_c_gapc.o $(objfolder)/gga_c_gaploc.o $(objfolder)/gga_c_hcth_a.o \
	$(objfolder)/gga_c_lm.o $(objfolder)/gga_c_lyp.o $(objfolder)/gga_c_op_b88.o $(objfolder)/gga_c_op_g96.o \
	$(objfolder)/gga_c_op_pbe.o $(objfolder)/gga_c_op_pw91.o $(objfolder)/gga_c_optc.o $(objfolder)/gga_c_op_xalpha.o \
	$(objfolder)/gga_c_p86.o $(objfolder)/gga_c_pbe.o $(objfolder)/gga_c_pbeloc.o $(objfolder)/gga_c_pw91.o \
	$(objfolder)/gga_c_q2d.o $(objfolder)/gga_c_regtpss.o $(objfolder)/gga_c_revtca.o $(objfolder)/gga_c_scan_e0.o \
	$(objfolder)/gga_c_sg4.o $(objfolder)/gga_c_sogga11.o $(objfolder)/gga_c_tca.o $(objfolder)/gga_c_w94.o \
	$(objfolder)/gga_c_wi.o $(objfolder)/gga_c_wl.o $(objfolder)/gga_c_zpbeint.o $(objfolder)/gga_c_zvpbeint.o \
	$(objfolder)/gga_k_dk.o $(objfolder)/gga_k_exp4.o $(objfolder)/gga_k_meyer.o $(objfolder)/gga_k_ol1.o \
	$(objfolder)/gga_k_ol2.o $(objfolder)/gga_k_pearson.o $(objfolder)/gga_k_tflw.o $(objfolder)/gga_k_thakkar.o \
	$(objfolder)/gga_x_2d_b86.o $(objfolder)/gga_x_2d_b86_mgc.o $(objfolder)/gga_x_2d_b88.o $(objfolder)/gga_x_2d_pbe.o \
	$(objfolder)/gga_x_airy.o $(objfolder)/gga_x_ak13.o $(objfolder)/gga_x_am05.o $(objfolder)/gga_x_b86.o \
	$(objfolder)/gga_x_b88.o $(objfolder)/gga_x_bayesian.o $(objfolder)/gga_x_beefvdw.o $(objfolder)/gga_x_bpccac.o \
	$(objfolder)/gga_x_c09x.o $(objfolder)/gga_x_cap.o $(objfolder)/gga_xc_b97.o $(objfolder)/gga_x_chachiyo.o \
	$(objfolder)/gga_xc_th1.o $(objfolder)/gga_xc_th2.o $(objfolder)/gga_xc_th3.o $(objfolder)/gga_x_dk87.o \
	$(objfolder)/gga_x_eg93.o $(objfolder)/gga_x_ft97.o $(objfolder)/gga_x_g96.o $(objfolder)/gga_x_hcth_a.o \
	$(objfolder)/gga_x_herman.o $(objfolder)/gga_x_hjs_b88_v2.o $(objfolder)/gga_x_hjs.o $(objfolder)/gga_x_htbs.o \
	$(objfolder)/gga_x_kt.o $(objfolder)/gga_x_lag.o $(objfolder)/gga_x_lg93.o $(objfolder)/gga_x_lv_rpw86.o \
	$(objfolder)/gga_x_mpbe.o $(objfolder)/gga_x_n12.o $(objfolder)/gga_x_optx.o $(objfolder)/gga_x_pbea.o \
	$(objfolder)/gga_x_pbe.o $(objfolder)/gga_x_pbeint.o $(objfolder)/gga_x_pbepow.o $(objfolder)/gga_x_pbetrans.o \
	$(objfolder)/gga_x_pw86.o $(objfolder)/gga_x_pw91.o $(objfolder)/gga_x_q2d.o $(objfolder)/gga_x_rge2.o \
	$(objfolder)/gga_x_rpbe.o $(objfolder)/gga_x_sg4.o $(objfolder)/gga_x_sogga11.o $(objfolder)/gga_x_ssb_sw.o \
	$(objfolder)/gga_x_vmt84.o $(objfolder)/gga_x_vmt.o $(objfolder)/gga_x_wc.o $(objfolder)/hyb_gga_xc_wb97.o \
	$(objfolder)/lda_c_1d_csc.o $(objfolder)/lda_c_1d_loos.o $(objfolder)/lda_c_2d_amgb.o $(objfolder)/lda_c_2d_prm.o \
	$(objfolder)/lda_c_chachiyo.o $(objfolder)/lda_c_gk72.o $(objfolder)/lda_c_gombas.o $(objfolder)/lda_c_hl.o \
	$(objfolder)/lda_c_lp96.o $(objfolder)/lda_c_ml1.o $(objfolder)/lda_c_pk09.o $(objfolder)/lda_c_pw.o \
	$(objfolder)/lda_c_pz.o $(objfolder)/lda_c_rc04.o $(objfolder)/lda_c_rpa.o $(objfolder)/lda_c_vwn_1.o \
	$(objfolder)/lda_c_vwn_2.o $(objfolder)/lda_c_vwn_3.o $(objfolder)/lda_c_vwn_4.o $(objfolder)/lda_c_vwn.o \
	$(objfolder)/lda_c_vwn_rpa.o $(objfolder)/lda_c_wigner.o $(objfolder)/lda_k_tf.o $(objfolder)/lda_k_zlp.o \
	$(objfolder)/lda_x_2d.o $(objfolder)/lda_xc_1d_ehwlrg.o $(objfolder)/lda_xc_ksdt.o $(objfolder)/lda_xc_teter93.o \
	$(objfolder)/lda_x.o $(objfolder)/lda_xc_zlp.o $(objfolder)/lda_x_rel.o 
#       $(objfolder)/lda_x_erf.o $(objfolder)/hyb_mgga_xc_wb97mv.o $(objfolder)/hyb_mgga_x_dldf.o $(objfolder)/hyb_mgga_x_m05.o \
	$(objfolder)/mgga_c_b88.o $(objfolder)/mgga_c_bc95.o $(objfolder)/mgga_c_cs.o $(objfolder)/mgga_c_kcis.o \
	$(objfolder)/mgga_c_m05.o $(objfolder)/mgga_c_m06l.o $(objfolder)/mgga_c_m08.o $(objfolder)/mgga_c_pkzb.o \
	$(objfolder)/mgga_c_revscan.o $(objfolder)/mgga_c_revtpss.o $(objfolder)/mgga_c_scan.o $(objfolder)/mgga_c_tpss.o \
	$(objfolder)/mgga_c_tpssloc.o $(objfolder)/mgga_c_vsxc.o $(objfolder)/mgga_k_pc07.o $(objfolder)/mgga_x_br89_explicit.o \
	$(objfolder)/mgga_xc_b97mv.o $(objfolder)/mgga_xc_b98.o $(objfolder)/mgga_xc_cc06.o $(objfolder)/mgga_xc_lp90.o \
	$(objfolder)/mgga_xc_zlp.o $(objfolder)/mgga_x_gvt4.o $(objfolder)/mgga_x_gx.o $(objfolder)/mgga_x_lta.o \
	$(objfolder)/mgga_x_m06l.o $(objfolder)/mgga_x_m08.o $(objfolder)/mgga_x_m11.o $(objfolder)/mgga_x_m11_l.o \
	$(objfolder)/mgga_x_mbeef.o $(objfolder)/mgga_x_mbeefvdw.o $(objfolder)/mgga_x_mk00.o $(objfolder)/mgga_x_mn12.o \
	$(objfolder)/mgga_x_ms.o $(objfolder)/mgga_x_mvs.o $(objfolder)/mgga_x_pbe_gx.o $(objfolder)/mgga_x_pkzb.o \
	$(objfolder)/mgga_x_sa_tpss.o $(objfolder)/mgga_x_scan.o $(objfolder)/mgga_x_tau_hcth.o $(objfolder)/mgga_x_tm.o \
	$(objfolder)/mgga_x_tpss.o $(objfolder)/mgga_x_vt84.o
cublasobj       = $(objfolder)/fortran_thunking.o
cusolverobj     = $(objfolder)/quick_cusolver.o 
#----------------------
# octree files
#----------------------
octobj    = $(objfolder)/grid_packer.o $(objfolder)/octree.o

#----------------------
# quick modules and object files
#----------------------

modobj= $(objfolder)/quick_mpi_module.o $(objfolder)/quick_constants_module.o $(objfolder)/quick_method_module.o \
        $(objfolder)/quick_molspec_module.o $(objfolder)/quick_gaussian_class_module.o $(objfolder)/quick_size_module.o \
        $(objfolder)/quick_amber_interface_module.o $(objfolder)/quick_basis_module.o $(objfolder)/quick_calculated_module.o \
        $(objfolder)/quick_divcon_module.o $(objfolder)/quick_ecp_module.o $(objfolder)/quick_electrondensity_module.o \
        $(objfolder)/quick_files_module.o $(objfolder)/quick_gridpoints_module.o $(objfolder)/quick_mfcc_module.o \
        $(objfolder)/quick_params_module.o $(objfolder)/quick_pb_module.o $(objfolder)/quick_scratch_module.o \
        $(objfolder)/quick_timer_module.o $(objfolder)/quick_scf_module.o $(objfolder)/quick_gradient_module.o \
	$(objfolder)/quick_all_module.o

MAIN = $(srcfolder)/main.o

OBJ =   $(objfolder)/initialize.o $(objfolder)/read_job_and_atom.o $(objfolder)/fmm.o \
        $(objfolder)/getMolSad.o $(objfolder)/getMol.o $(objfolder)/shell.o $(objfolder)/schwarz.o \
        $(objfolder)/quick_one_electron_integral.o $(objfolder)/getEnergy.o $(objfolder)/inidivcon.o \
        $(objfolder)/ecp.o $(objfolder)/hfoperator.o $(objfolder)/nuclear.o \
        $(objfolder)/dft.o $(objfolder)/sedftoperator.o \
        $(objfolder)/scf.o $(objfolder)/uscf.o $(objfolder)/finalize.o $(objfolder)/uhfoperator.o \
        $(objfolder)/udftoperator.o $(objfolder)/usedftoperator.o \
        $(objfolder)/uelectdii.o $(objfolder)/mpi_setup.o $(objfolder)/quick_debug.o \
        $(objfolder)/calMP2.o $(objfolder)/optimize.o $(objfolder)/gradient.o $(objfolder)/hessian.o \
        $(objfolder)/CPHF.o $(objfolder)/frequency.o $(objfolder)/MFCC.o $(objfolder)/basis.o \
        $(objfolder)/fake_amber_interface.o $(objfolder)/scf_operator.o  

SUBS = $(objfolder)/Angles.o $(objfolder)/copyDMat.o $(objfolder)/copySym.o \
	$(objfolder)/degen.o $(objfolder)/denspt.o $(objfolder)/diag.o $(objfolder)/dipole.o \
	$(objfolder)/EffChar.o $(objfolder)/eigvec.o \
	$(objfolder)/ekinetic.o $(objfolder)/findBlock.o $(objfolder)/fmt.o $(objfolder)/getinum.o \
	$(objfolder)/getNum.o $(objfolder)/greedy_distrubute.o $(objfolder)/hrr.o $(objfolder)/iatoi.o \
	$(objfolder)/iatoimp.o $(objfolder)/io.o $(objfolder)/iwhole.o \
	$(objfolder)/lbfgs.o $(objfolder)/Lsolve.o $(objfolder)/matComp.o $(objfolder)/matMul.o \
	$(objfolder)/order.o $(objfolder)/orthog.o $(objfolder)/PriCol.o $(objfolder)/PriSym.o \
	$(objfolder)/PrtAct.o $(objfolder)/PrtDat.o $(objfolder)/PrtErr.o $(objfolder)/PrtLab.o \
	$(objfolder)/PrtMsg.o $(objfolder)/PrtTim.o $(objfolder)/PrtWrn.o $(objfolder)/pteval.o \
	$(objfolder)/quick_open.o $(objfolder)/random.o $(objfolder)/rdinum.o $(objfolder)/rdnml.o \
	$(objfolder)/rdnum.o $(objfolder)/rdword.o $(objfolder)/readPDB.o $(objfolder)/spdfgh.o \
	$(objfolder)/ssw.o $(objfolder)/sum2Mat.o $(objfolder)/transpose.o $(objfolder)/tridi.o \
	$(objfolder)/upcase.o $(objfolder)/vett.o $(objfolder)/whatis.o $(objfolder)/whole.o \
	$(objfolder)/wrtRestart.o $(objfolder)/xnorm.o $(objfolder)/zeroMatrix.o $(objfolder)/zmake.o \
	$(objfolder)/pt2der.o $(objfolder)/sswder.o $(objfolder)/denspt_new_imp.o \
	$(objfolder)/pteval_new_imp.o $(objfolder)/scaMatMul.o

TESTAPI=$(srcfolder)/quick_api_test.o

all: quick quick.cuda
#************************************************************************
# 
#                  B. Make necessary directories
#  
#************************************************************************
makefolders:
	mkdir -p $(objfolder) $(exefolder) $(libfolder)

#************************************************************************
# 
#                 C. Make Object Files
# 
#************************************************************************

#================= common subroutine library ============================
quick_subs:
	cd $(subfolder) && make all

#================= quick module library =================================
quick_modules:
	cd $(modfolder) && make all
#================= octree subroutines   =================================
octree:
	cd $(octfolder) && make all
#============= targets for cuda =========================================
quick_cuda:
	cd $(cudafolder) && make allbutxc
	cd $(cudafolder) && make xc 
		
#================= targets for BLAS =====================================
blas:
	cd $(blasfolder) && make
#==================== libxc cpu library =================================
libxc_cpu:
	cd $(libxcfolder) && make libxc_cpu
#==================== libxc cpu library =================================
libxc_gpu:
	cd $(libxcfolder) && make libxc_gpu
	cd $(libxcfolder)/maple2c_device && make all
#=============== targets for CUBLAS =====================================

$(cublasobj):$(objfolder)/%.o:$(cublasfolder)/%.c
	$(CC) $(CPP_FLAG) -c $< -o $@

#=============== targets for CUSOLVER ===================================

$(cusolverobj):$(objfolder)/%.o:$(cusolverfolder)/%.c
	$(CC) $(CPP_FLAG) -c $< -o $@

#===================== target for general src files =====================

$(OBJ):$(objfolder)/%.o:$(srcfolder)/%.f90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -I$(objfolder) -c $< -o $@
#===================== target for main src files ========================
$(MAIN):$(srcfolder)/%.o:$(srcfolder)/%.f90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -I$(objfolder) -c $< -o $@
#===================== target for compiling api test ====================
$(TESTAPI):$(srcfolder)/%.o:$(srcfolder)/%.f90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -I$(objfolder) -c $< -o $@
#==================== target configuration files ========================
cpconfig:
	cp $(configfolder)/config.h $(srcfolder)/config.h
cpconfig.cuda:
	cp $(configfolder)/config.cuda.h $(srcfolder)/config.h
cpconfig.cuda.SP:
	cp $(configfolder)/config.cuda.SP.h $(srcfolder)/config.h
cpconfig.cuda.MPI:
	cp $(configfolder)/config.cuda.MPI.h $(srcfolder)/config.h
cpconfig.MPI:
	cp $(configfolder)/config.MPI.h $(srcfolder)/config.h

#**********************************************************************
# 
#                 C. Make Executables
# 
#**********************************************************************

ifeq ($(INSTALLER), serial)
install: cpconfig libxc_cpu octree quick_modules quick_subs $(OBJ) $(MAIN) blas 
	$(ARCH) $(ARCHFLAGS) $(libfolder)/libquick.$(LIBEXT) $(objfolder)/*.o
	$(FC) -o $(exefolder)/quick $(MAIN) -L$(libfolder) -lquick -lblas -lxc $(LDFLAGS)
endif

ifeq ($(INSTALLER), cuda)
install: cpconfig.cuda libxc_gpu octree quick_cuda quick_modules quick_subs $(OBJ) $(MAIN) $(cusolverobj) $(cublasobj)
	$(ARCH) $(ARCHFLAGS) $(libfolder)/libquickcu.$(LIBEXT) $(objfolder)/*.o
	$(FC) -o $(exefolder)/quick.cuda $(MAIN) -L$(libfolder) -lquickcu -lxc $(CFLAGS) $(LDFLAGS)
endif

ifeq ($(INSTALLER), mpi)
install: cpconfig.MPI libxc_cpu octree quick_modules quick_subs $(OBJ) $(MAIN) blas
	$(ARCH) $(ARCHFLAGS) $(libfolder)/libquickmpi.$(LIBEXT) $(objfolder)/*.o
	$(FC) -o $(exefolder)/quick.MPI $(MAIN) -L$(libfolder) -lquickmpi -lblas -lxc $(LDFLAGS)
endif

ifeq ($(INSTALLER), cudampi)
install: cpconfig.cuda.MPI libxc_gpu octree quick_cuda quick_modules quick_subs $(OBJ) $(MAIN) blas $(cusolverobj) $(cublasobj)
	$(ARCH) $(ARCHFLAGS) $(libfolder)/libquickcumpi.$(LIBEXT) $(objfolder)/*.o
	$(FC) -o $(exefolder)/quick.cuda.MPI $(MAIN) -L$(libfolder) -lquickcumpi -lblas -lxc $(CFLAGS) $(LDFLAGS)
endif

#quick.cuda.SP: makefolders cpconfig.cuda.SP quick_cuda quick_modules quick_subs quick_pprs $(OBJ) $(cusolverobj) $(cublasobj)
#	$(FC) -o quick.cuda.SP $(OBJ) $(modobj) $(cudaobj) $(SUBS) $(cusolverobj) $(cublasobj) $(CFLAGS) 

quicklib: quick $(TESTAPI)
	$(FC) -o $(exefolder)/testapi.o $(TESTAPI) -I$(objfolder) -L$(libfolder) -lquick -lblas -lxc $(LDFLAGS) 

quickculib: quick.cuda $(TESTAPI)
	$(FC) -o $(exefolder)/testapi.o $(TESTAPI) -I$(objfolder) -L$(libfolder) -lquickcu -lxc $(CFLAGS) $(LDFLAGS) 

quickmpilib: quick.MPI $(TESTAPI)
	$(FC) -DQUAPI_MPIV -o $(exefolder)/testapi.o $(TESTAPI) -I$(objfolder) -L$(libfolder) -lquickmpi -lblas -lxc $(LDFLAGS) 

quickcumpilib: quick.cuda.MPI $(TESTAPI) 
	$(FC) -DQUAPI_MPIV -o $(exefolder)/testapi.o $(TESTAPI) -I$(objfolder) -L$(libfolder) -lquickcumpi -lblas -lxc $(CFLAGS) $(LDFLAGS) 
		
#************************************************************************
# 
#                 D. Self-defined Option
# 
#************************************************************************

# - 1. Clean object files
clean: 
	-rm -f $(objfolder)/* $(srcfolder)/*.o 
	cd $(cudafolder) && make clean
	cd $(subfolder) && make clean
	cd $(blasfolder) && make clean
	cd $(modfolder) && make clean
	cd $(libxcfolder) && make clean
	cd $(libxcfolder)/maple2c_device && make clean	

uninstall: 
	-rm -f $(exefolder)/quick* $(exefolder)/testapi.o $(objfolder)/* $(libfolder)/* $(srcfolder)/*.o 
	cd $(cudafolder) && make clean
	cd $(subfolder) && make clean
	cd $(blasfolder) && make clean
	cd $(modfolder) && make clean
	cd $(libxcfolder) && make clean
	cd $(libxcfolder)/maple2c_device && make clean

# - 2. Make tags for source files
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)

include $(srcfolder)/depend 
