module quick_c_interface
    ! ------------------------------------------------------------------------
    ! gpu_get_oshell
    ! ------------------------------------------------------------------------

    interface
        subroutine gpu_get_oshell_eri(o, ob) &
            BIND(C, NAME="gpu_get_oshell_eri_")
            ! void gpu_get_oshell_eri_(QUICKDouble* o, QUICKDouble* ob);
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: o(*)
            real(kind=c_double) :: ob(*)
        end subroutine gpu_get_oshell_eri
    end interface

    interface
        subroutine gpu_get_oshell_eri_grad(grad) &
            BIND(C, NAME="gpu_get_oshell_eri_grad_")
            ! void gpu_get_oshell_eri_grad_(QUICKDouble* grad)
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: grad(*)
        end subroutine gpu_get_oshell_eri_grad
    end interface

    interface
        subroutine gpu_get_oshell_xc(Eelxc, aelec, belec, o, ob) &
            BIND(C, NAME="gpu_get_oshell_xc_")
            ! void gpu_get_oshell_xc_(QUICKDouble* Eelxc, QUICKDouble* aelec, QUICKDouble* belec, QUICKDouble *o, QUICKDouble *ob);
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: Eelxc
            real(kind=c_double) :: aelec
            real(kind=c_double) :: belec
            real(kind=c_double) :: o(*)
            real(kind=c_double) :: ob(*)
        end subroutine gpu_get_oshell_xc
    end interface

    interface
        subroutine gpu_get_oshell_xcgrad(grad) &
            BIND(C, NAME="gpu_get_oshell_xcgrad_")
            ! void gpu_get_oshell_xcgrad_(QUICKDouble *grad)
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: grad(*)
        end subroutine gpu_get_oshell_xcgrad
    end interface

    ! ------------------------------------------------------------------------
    ! gpu_get_cshell
    ! ------------------------------------------------------------------------

    interface
        subroutine gpu_get_cshell_eri(o) &
            BIND(C, NAME="gpu_get_cshell_eri_")
            ! void gpu_get_cshell_eri_(QUICKDouble* o);
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: o(*)
        end subroutine gpu_get_cshell_eri
    end interface

    interface
        subroutine gpu_get_cshell_eri_grad(grad) &
            BIND(C, NAME="gpu_get_cshell_eri_grad_")
            ! void gpu_get_cshell_eri_grad_(QUICKDouble* grad)
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: grad(*)
        end subroutine gpu_get_cshell_eri_grad
    end interface

    interface
        subroutine gpu_get_cshell_xc(Eelxc, aelec, belec, o) &
            BIND(C, NAME="gpu_get_cshell_xc_")
            ! void gpu_get_cshell_xc_(QUICKDouble* Eelxc, QUICKDouble* aelec, QUICKDouble* belec, QUICKDouble *o)
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: Eelxc
            real(kind=c_double) :: aelec
            real(kind=c_double) :: belec
            real(kind=c_double) :: o(*)
        end subroutine gpu_get_cshell_xc
    end interface

    interface
        subroutine gpu_get_cshell_xcgrad(grad) &
            BIND(C, NAME="gpu_get_cshell_xcgrad_")
            ! void gpu_get_cshell_xcgrad_(QUICKDouble *grad);
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: grad(*)
        end subroutine gpu_get_cshell_xcgrad
    end interface

    ! ------------------------------------------------------------------------
    ! gpu_get_<other>
    ! ------------------------------------------------------------------------

    interface
        subroutine gpu_get_ssw(gridx, gridy, gridz, wtang, rwt, rad3, sswt, weight, gatm, count) &
            BIND(C, NAME="gpu_get_ssw_")
            ! void gpu_get_ssw_(QUICKDouble *gridx, QUICKDouble *gridy, QUICKDouble *gridz, QUICKDouble *wtang, QUICKDouble *rwt, QUICKDouble *rad3, QUICKDouble *sswt, QUICKDouble *weight, int *gatm, int *count){
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double, c_int
            real(kind=c_double) :: gridx(*)
            real(kind=c_double) :: gridy(*)
            real(kind=c_double) :: gridz(*)
            real(kind=c_double) :: wtang(*)
            real(kind=c_double) :: rwt(*)
            real(kind=c_double) :: rad3(*)
            real(kind=c_double) :: sswt(*)
            real(kind=c_double) :: weight(*)
            integer(kind=c_int) :: gatm(*)
            integer(kind=c_int) :: count
        end subroutine gpu_get_ssw
    end interface

    ! ------------------------------------------------------------------------
    ! gpu_upload
    ! ------------------------------------------------------------------------

    interface
        subroutine gpu_upload_method(quick_method, is_oshell, hyb_coeff) &
            BIND(C, NAME="gpu_upload_method_")
            ! void gpu_upload_method_(int* quick_method, bool* is_oshell, double* hyb_coeff);
            ! ATTENTION: this must be adjusted to QUICKDouble
            ! MSOEGTROP ToDo: DataType does not match
            use, intrinsic :: iso_c_binding, only : c_double, c_int, c_bool
            integer(kind=c_int) :: quick_method
            logical(kind=c_int) :: is_oshell
            real(kind=c_double) :: hyb_coeff
        end subroutine gpu_upload_method
    end interface

    interface
        subroutine gpu_upload_xyz(atom_xyz) &
            BIND(C, NAME="gpu_upload_xyz_")
            ! void gpu_upload_xyz_(QUICKDouble* atom_xyz)
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: atom_xyz(*)
        end subroutine gpu_upload_xyz
    end interface

    interface
        subroutine gpu_upload_atom_and_chg(atom, atom_chg) &
            BIND(C, NAME="gpu_upload_atom_and_chg_")
            ! void gpu_upload_atom_and_chg_(int* atom, QUICKDouble* atom_chg);
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double, c_int
            integer(kind=c_int) :: atom(*)
            real(kind=c_double) :: atom_chg(*)
        end subroutine gpu_upload_atom_and_chg
    end interface

    interface
        subroutine gpu_upload_calculated(o, co, vec, dense) &
            BIND(C, NAME="gpu_upload_calculated_")
            ! void gpu_upload_calculated_(QUICKDouble* o, QUICKDouble* co, QUICKDouble* vec, QUICKDouble* dense);
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double, c_int
            real(kind=c_double) :: o(*)
            real(kind=c_double) :: co(*)
            real(kind=c_double) :: vec(*)
            real(kind=c_double) :: dense(*)
        end subroutine gpu_upload_calculated
    end interface

    interface
        subroutine gpu_upload_calculated_beta(ob, denseb) &
            BIND(C, NAME="gpu_upload_calculated_beta_")
            ! void gpu_upload_calculated_beta_(QUICKDouble* ob, QUICKDouble* denseb);
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double, c_int
            real(kind=c_double) :: ob(*)
            real(kind=c_double) :: denseb(*)
        end subroutine gpu_upload_calculated_beta
    end interface

    interface
        subroutine gpu_upload_basis( &
            nshell, nprim, jshell, jbasis, maxcontract, &
            ncontract, itype, aexp, dcoeff, &
            first_basis_function, last_basis_function, first_shell_basis_function, last_shell_basis_function, &
            ncenter, kstart, katom, ktype, kprim, kshell, Ksumtype, &
            Qnumber, Qstart, Qfinal, Qsbasis, Qfbasis, &
            gccoeff, cons, gcexpo, KLMN) &
            BIND(C, NAME="gpu_upload_basis_")
            ! void gpu_upload_basis_(int* nshell, int* nprim, int* jshell, int* jbasis, int* maxcontract,
            !     int* ncontract, int* itype,     QUICKDouble* aexp,      QUICKDouble* dcoeff,
            !     int* first_basis_function, int* last_basis_function, int* first_shell_basis_function, int* last_shell_basis_function,
            !     int* ncenter,   int* kstart,    int* katom,     int* ktype,     int* kprim,  int* kshell, int* Ksumtype,
            !     int* Qnumber,   int* Qstart,    int* Qfinal,    int* Qsbasis,   int* Qfbasis,
            !     QUICKDouble* gccoeff,           QUICKDouble* cons,      QUICKDouble* gcexpo, int* KLMN)
            use, intrinsic :: iso_c_binding, only : c_double, c_int
            integer(kind=c_int) :: nshell
            integer(kind=c_int) :: nprim
            integer(kind=c_int) :: jshell
            integer(kind=c_int) :: jbasis
            integer(kind=c_int) :: maxcontract
            integer(kind=c_int) :: ncontract(*)
            integer(kind=c_int) :: itype(*)
            real(kind=c_double) :: aexp(*)
            real(kind=c_double) :: dcoeff(*)
            integer(kind=c_int) :: first_basis_function(*)
            integer(kind=c_int) :: last_basis_function(*)
            integer(kind=c_int) :: first_shell_basis_function(*)
            integer(kind=c_int) :: last_shell_basis_function(*)
            integer(kind=c_int) :: ncenter(*)
            integer(kind=c_int) :: kstart(*)
            integer(kind=c_int) :: katom(*)
            integer(kind=c_int) :: ktype(*)
            integer(kind=c_int) :: kprim(*)
            integer(kind=c_int) :: kshell(*)
            integer(kind=c_int) :: Ksumtype(*)
            integer(kind=c_int) :: Qnumber(*)
            integer(kind=c_int) :: Qstart(*)
            integer(kind=c_int) :: Qfinal(*)
            integer(kind=c_int) :: Qsbasis(*)
            integer(kind=c_int) :: Qfbasis(*)
            real(kind=c_double) :: gccoeff(*)
            real(kind=c_double) :: cons(*)
            real(kind=c_double) :: gcexpo(*)
            integer(kind=c_int) :: KLMN(*)
        end subroutine gpu_upload_basis
    end interface

    interface
        subroutine gpu_upload_cutoff(cutMatrix, integralCutoff, primLimit, DMCutoff) &
            BIND(C, NAME="gpu_upload_cutoff_")
            ! void gpu_upload_cutoff_(QUICKDouble* cutMatrix, QUICKDouble* integralCutoff,QUICKDouble* primLimit, QUICKDouble* DMCutoff)
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: cutMatrix(*)
            real(kind=c_double) :: integralCutoff
            real(kind=c_double) :: primLimit
            real(kind=c_double) :: DMCutoff
        end subroutine gpu_upload_cutoff
    end interface

    interface
        subroutine gpu_upload_cutoff_matrix(YCutoff, cutPrim) &
            BIND(C, NAME="gpu_upload_cutoff_matrix_")
            ! void gpu_upload_cutoff_matrix_(QUICKDouble* YCutoff,QUICKDouble* cutPrim)
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: YCutoff(*)
            real(kind=c_double) :: cutPrim(*)
        end subroutine gpu_upload_cutoff_matrix
    end interface

    interface
        subroutine gpu_upload_density_matrix(dense) &
            BIND(C, NAME="gpu_upload_density_matrix_")
            ! void gpu_upload_density_matrix_(QUICKDouble* dense);
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: dense(*)
        end subroutine gpu_upload_density_matrix
    end interface

    interface
        subroutine gpu_upload_grad(gradCutoff) &
            BIND(C, NAME="gpu_upload_grad_")
            ! void gpu_upload_grad_(QUICKDouble* gradCutoff);
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: gradCutoff
        end subroutine gpu_upload_grad
    end interface

    interface
        subroutine gpu_upload_libxc(nof_functionals, functional_id, xc_polarization, ierr) &
            BIND(C, NAME="gpu_upload_libxc_")
            ! void gpu_upload_libxc_(int* nof_functionals, int* functional_id, int* xc_polarization, int *ierr)
            use, intrinsic :: iso_c_binding, only : c_int
            integer(kind=c_int) :: nof_functionals
            integer(kind=c_int) :: functional_id(*)
            integer(kind=c_int) :: xc_polarization
            integer(kind=c_int) :: ierr
        end subroutine gpu_upload_libxc
    end interface

    interface
        subroutine gpu_upload_dft_grid(gridxb, gridyb, gridzb, gridb_sswt, gridb_weight, gridb_atm, bin_locator, basf, primf, basf_counter, primf_counter, bin_counter,gridb_count, nbins, nbtotbf, nbtotpf, isg, sigrad2, DMCutoff) &
            BIND(C, NAME="gpu_upload_dft_grid_")
            ! void gpu_upload_dft_grid_(QUICKDouble *gridxb, QUICKDouble *gridyb, QUICKDouble *gridzb, QUICKDouble *gridb_sswt, QUICKDouble *gridb_weight, int *gridb_atm, int *bin_locator, int *basf, int *primf, int *basf_counter, int *primf_counter, int *bin_counter,int *gridb_count, int *nbins, int *nbtotbf, int *nbtotpf, int *isg, QUICKDouble *sigrad2, QUICKDouble *DMCutoff){
            ! ATTENTION: this must be adjusted to QUICKDouble
            use, intrinsic :: iso_c_binding, only : c_double, c_int
            real(kind=c_double) :: gridxb(*)
            real(kind=c_double) :: gridyb(*)
            real(kind=c_double) :: gridzb(*)
            real(kind=c_double) :: gridb_sswt(*)
            real(kind=c_double) :: gridb_weight(*)
            integer(kind=c_int) :: gridb_atm(*)
            integer(kind=c_int) :: bin_locator(*)
            integer(kind=c_int) :: basf(*)
            integer(kind=c_int) :: primf(*)
            integer(kind=c_int) :: basf_counter(*)
            integer(kind=c_int) :: primf_counter(*)
            integer(kind=c_int) :: bin_counter(*)
            integer(kind=c_int) :: gridb_count
            integer(kind=c_int) :: nbins
            integer(kind=c_int) :: nbtotbf
            integer(kind=c_int) :: nbtotpf
            integer(kind=c_int) :: isg
            real(kind=c_double) :: sigrad2(*)
            real(kind=c_double) :: DMCutoff
        end subroutine gpu_upload_dft_grid
    end interface

    ! ------------------------------------------------------------------------
    ! gpu init
    ! ------------------------------------------------------------------------

    interface
        subroutine gpu_set_device(gpu_dev_id, ierr) &
            BIND(C, NAME="gpu_set_device_")
            ! void gpu_set_device_(int* gpu_dev_id, int* ierr);
            use, intrinsic :: iso_c_binding, only : c_int
            integer(kind=c_int) :: gpu_dev_id
            integer(kind=c_int) :: ierr
        end subroutine gpu_set_device
    end interface

    interface
        subroutine gpu_init(ierr) &
            BIND(C, NAME="gpu_init_")
            ! void gpu_init_(int* ierr);
            use, intrinsic :: iso_c_binding, only : c_int
            integer(kind=c_int) :: ierr
        end subroutine gpu_init
    end interface

    interface
        subroutine gpu_setup(natom, nbasis, nElec, imult, molchg, iAtomType) &
            BIND(C, NAME="gpu_setup_")
            ! void gpu_setup_(int* natom, int* nbasis, int* nElec, int* imult, int* molchg, int* iAtomType)
            use, intrinsic :: iso_c_binding, only : c_int
            integer(kind=c_int) :: natom
            integer(kind=c_int) :: nbasis
            integer(kind=c_int) :: nElec
            integer(kind=c_int) :: imult
            integer(kind=c_int) :: molchg
            integer(kind=c_int) :: iAtomType
        end subroutine gpu_setup
    end interface

    interface
        subroutine gpu_reupload_dft_grid() &
            BIND(C, NAME="gpu_reupload_dft_grid_")
        end subroutine gpu_reupload_dft_grid
    end interface

    interface
        subroutine gpu_startup(ierr) &
            BIND(C, NAME="gpu_startup_")
            ! void gpu_startup_(int* ierr)
            use, intrinsic :: iso_c_binding, only : c_int
            integer(kind=c_int) :: ierr
        end subroutine gpu_startup
    end interface

    ! ------------------------------------------------------------------------
    ! gpu shutdown
    ! ------------------------------------------------------------------------

    interface
        subroutine gpu_delete_dft_grid() &
            BIND(C, NAME="gpu_delete_dft_grid")
        end subroutine gpu_delete_dft_grid
    end interface

    interface
        subroutine gpu_delete_dft_dev_grid() &
            BIND(C, NAME="gpu_delete_dft_dev_grid")
        end subroutine gpu_delete_dft_dev_grid
    end interface

    interface
        subroutine gpu_delete_libxc(ierr) &
            BIND(C, NAME="gpu_delete_libxc_")
            ! void gpu_delete_libxc_(int *ierr)
            use, intrinsic :: iso_c_binding, only : c_int
            integer(kind=c_int) :: ierr
        end subroutine gpu_delete_libxc
    end interface

    interface
        subroutine gpu_cleanup() &
            BIND(C, NAME="gpu_cleanup_")
        end subroutine gpu_cleanup
    end interface

    interface
        subroutine gpu_shutdown(ierr) &
            BIND(C, NAME="gpu_shutdown_")
            ! void gpu_shutdown_(int* ierr)
            use, intrinsic :: iso_c_binding, only : c_int
            integer(kind=c_int) :: ierr
        end subroutine gpu_shutdown
    end interface

    ! ------------------------------------------------------------------------
    ! gpu diagnostics
    ! ------------------------------------------------------------------------

    interface
        subroutine gpu_get_device_info(gpu_dev_count, gpu_dev_id, gpu_dev_mem, gpu_num_proc, gpu_core_freq, gpu_dev_name, name_len, majorv, minorv, ierr) &
            BIND(C, NAME="gpu_get_device_info_")
            ! void gpu_get_device_info_(int* gpu_dev_count, int* gpu_dev_id,int* gpu_dev_mem,
            !            int* gpu_num_proc,double* gpu_core_freq,char* gpu_dev_name,int* name_len, int* majorv, int* minorv, int* ierr);
            use, intrinsic :: iso_c_binding, only : c_double, c_int, c_char
            integer(kind=c_int) :: gpu_dev_count
            integer(kind=c_int) :: gpu_dev_id
            integer(kind=c_int) :: gpu_dev_mem
            integer(kind=c_int) :: gpu_num_proc
            real(kind=c_double) :: gpu_core_freq
            character(kind=c_char) :: gpu_dev_name(name_len)
            integer(kind=c_int) :: name_len
            integer(kind=c_int) :: majorv
            integer(kind=c_int) :: minorv
            integer(kind=c_int) :: ierr
        end subroutine gpu_get_device_info
    end interface

    interface
        subroutine cuda_diag(o, x, hold, E, idegen, vec, co, V2, nbasis) &
            BIND(C, NAME="cuda_diag_")
            ! void CUDA_DIAG (double* o, const double* x,double* hold,
            !       const double* E, const double* idegen,
            !       const double* vec, const double* co,
            !       const double* V2, const int* nbasis)
            ! call cuda_diag(quick_scratch%hold, quick_scratch%tmpx,quick_scratch%tmphold,&
            ! quick_scratch%Sminhalf, quick_scratch%IDEGEN1, quick_scratch%hold2,quick_scratch%tmpco, quick_scratch%V, nbasis)
            ! MSOEGTROP ToDo: IDEGEN HAS WRONG TYPE
            use, intrinsic :: iso_c_binding, only : c_double, c_int
            real(kind=c_double) :: o(*)
            real(kind=c_double) :: x(*)
            real(kind=c_double) :: hold(*)
            real(kind=c_double) :: E(*)
            integer(kind=c_int) :: idegen(*)
            real(kind=c_double) :: vec(*)
            real(kind=c_double) :: co(*)
            real(kind=c_double) :: V2(*)
            integer(kind=c_int) :: nbasis
        end subroutine cuda_diag
    end interface

    interface
        subroutine cuda_diag_idegenf(o, x, hold, E, idegen, vec, co, V2, nbasis) &
            BIND(C, NAME="cuda_diag_")
            ! void CUDA_DIAG (double* o, const double* x,double* hold,
            !       const double* E, const double* idegen,
            !       const double* vec, const double* co,
            !       const double* V2, const int* nbasis)
            ! call cuda_diag(quick_scratch%hold, quick_scratch%tmpx,quick_scratch%tmphold,&
            ! quick_scratch%Sminhalf, quick_scratch%IDEGEN1, quick_scratch%hold2,quick_scratch%tmpco, quick_scratch%V, nbasis)
            ! MSOEGTROP ToDo: IDEGEN HAS WRONG TYPE
            use, intrinsic :: iso_c_binding, only : c_double, c_int
            real(kind=c_double) :: o(*)
            real(kind=c_double) :: x(*)
            real(kind=c_double) :: hold(*)
            real(kind=c_double) :: E(*)
            real(kind=c_double) :: idegen(*)
            real(kind=c_double) :: vec(*)
            real(kind=c_double) :: co(*)
            real(kind=c_double) :: V2(*)
            integer(kind=c_int) :: nbasis
        end subroutine cuda_diag_idegenf
    end interface

    ! ------------------------------------------------------------------------
    ! mgpu
    ! ------------------------------------------------------------------------

    interface
        subroutine mgpu_init(mpirank, mpisize, device, ierr) &
            BIND(C, NAME="mgpu_init_")
            ! void mgpu_init_(int *mpirank, int *mpisize, int *device, int* ierr)
            use, intrinsic :: iso_c_binding, only : c_int
            integer(kind=c_int) :: mpirank
            integer(kind=c_int) :: mpisize
            integer(kind=c_int) :: device
            integer(kind=c_int) :: ierr
        end subroutine mgpu_init
    end interface

    interface
        subroutine mgpu_shutdown(ierr) &
            BIND(C, NAME="mgpu_shutdown_")
            ! void mgpu_shutdown_(int* ierr)
            use, intrinsic :: iso_c_binding, only : c_int
            integer(kind=c_int) :: ierr
        end subroutine mgpu_shutdown
    end interface

    interface
        subroutine mgpu_get_device_info(dev_id, gpu_dev_mem, gpu_num_proc,gpu_core_freq,gpu_dev_name,name_len, majorv, minorv) &
            BIND(C, NAME="mgpu_get_device_info_")
            ! void mgpu_get_device_info_(int* dev_id,int* gpu_dev_mem,
            !     int* gpu_num_proc,double* gpu_core_freq,char* gpu_dev_name,int* name_len, int* majorv, int* minorv)
            use, intrinsic :: iso_c_binding, only : c_double, c_int, c_char
            integer(kind=c_int) :: dev_id
            integer(kind=c_int) :: gpu_dev_mem
            integer(kind=c_int) :: gpu_num_proc
            real(kind=c_double) :: gpu_core_freq
            character(kind=c_char) :: gpu_dev_name(*)
            integer(kind=c_int) :: name_len
            integer(kind=c_int) :: majorv
            integer(kind=c_int) :: minorv
        end subroutine mgpu_get_device_info
    end interface

    interface
        subroutine mgpu_query(mpisize, mpirank, mgpu_id, ierr) &
            BIND(C, NAME="mgpu_query_")
            ! void mgpu_query_(int* mpisize, int *mpirank, int *mgpu_id, int* ierr)
            use, intrinsic :: iso_c_binding, only : c_int
            integer(kind=c_int) :: mpisize
            integer(kind=c_int) :: mpirank
            integer(kind=c_int) :: mgpu_id
            integer(kind=c_int) :: ierr
        end subroutine mgpu_query
    end interface

    interface
        subroutine mgpu_get_xclb_time(t_xclb) &
            BIND(C, NAME="mgpu_get_xclb_time_")
            ! void mgpu_get_xclb_time_(double *t_xclb){
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: t_xclb
        end subroutine mgpu_get_xclb_time
    end interface

    interface
        subroutine mgpu_get_xcrb_time(t_xcrb, t_xcpg) &
            BIND(C, NAME="mgpu_get_xcrb_time_")
            ! void mgpu_get_xcrb_time_(double* t_xcrb, double* t_xcpg){
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: t_xcrb
            real(kind=c_double) :: t_xcpg
        end subroutine mgpu_get_xcrb_time
    end interface

    interface
        subroutine mgpu_get_2elb_time(t_2elb) &
            BIND(C, NAME="mgpu_get_2elb_time_")
            ! void mgpu_get_2elb_time_(double* t_2elb){
            use, intrinsic :: iso_c_binding, only : c_double
            real(kind=c_double) :: t_2elb
        end subroutine mgpu_get_2elb_time
    end interface

    ! ------------------------------------------------------------------------
    ! gpack library
    ! ------------------------------------------------------------------------

    interface
        subroutine gpack_initialize() &
            BIND(C, NAME="gpack_initialize_")
        ! void gpack_initialize_(){
        end subroutine gpack_initialize
    end interface

    interface
        subroutine gpack_finalize() &
            BIND(C, NAME="gpack_finalize_")
        ! void gpack_finalize_(){
        end subroutine gpack_finalize
    end interface

    interface
        subroutine gpack_pack_pts(grid_ptx, grid_pty, grid_ptz, grid_atm, grid_sswt, grid_weight, &
            arr_size, natoms, nbasis, maxcontract, DMCutoff, sigrad2, ncontract, aexp, dcoeff, &
            ncenter, itype, xyz, ngpts, nbins, nbtotbf, nbtotpf, toct, tprscrn) &
            BIND(C, NAME="gpack_pack_pts_")
            ! void gpack_pack_pts_(double *grid_ptx, double *grid_pty, double *grid_ptz, int *grid_atm, double *grid_sswt, double *grid_weight, int *arr_size, int *natoms, int *nbasis, int *maxcontract, double *DMCutoff, double *sigrad2, int *ncontract, double *aexp, double *dcoeff, int *ncenter, int *itype, double *xyz, int *ngpts, int *nbins, int *nbtotbf, int *nbtotpf, double *toct, double *tprscrn);
            use, intrinsic :: iso_c_binding, only : c_double, c_int
            real(kind=c_double) :: grid_ptx(*)
            real(kind=c_double) :: grid_pty(*)
            real(kind=c_double) :: grid_ptz(*)
            integer(kind=c_int) :: grid_atm(*)
            real(kind=c_double) :: grid_sswt(*)
            real(kind=c_double) :: grid_weight(*)
            integer(kind=c_int) :: arr_size
            integer(kind=c_int) :: natoms
            integer(kind=c_int) :: nbasis
            integer(kind=c_int) :: maxcontract
            real(kind=c_double) :: DMCutoff
            real(kind=c_double) :: sigrad2(*)
            integer(kind=c_int) :: ncontract(*)
            real(kind=c_double) :: aexp(*)
            real(kind=c_double) :: dcoeff(*)
            integer(kind=c_int) :: ncenter(*)
            integer(kind=c_int) :: itype(*)
            real(kind=c_double) :: xyz(*)
            integer(kind=c_int) :: ngpts
            integer(kind=c_int) :: nbins
            integer(kind=c_int) :: nbtotbf
            integer(kind=c_int) :: nbtotpf
            real(kind=c_double) :: toct
            real(kind=c_double) :: tprscrn
        end subroutine gpack_pack_pts
    end interface

    interface
        subroutine get_gpu_grid_info(gridx, gridy, gridz, ssw, weight, atm, bin_locator, basf, primf, basf_counter, primf_counter, bin_counter) &
            BIND(C, NAME="get_gpu_grid_info_")
            ! void get_gpu_grid_info_(double *gridx, double *gridy, double *gridz, double *ssw, double *weight, int *atm, int *bin_locator, int *basf, int *primf, int *basf_counter, int *primf_counter, int *bin_counter);
            use, intrinsic :: iso_c_binding, only : c_double, c_int
            real(kind=c_double) :: gridx(*)
            real(kind=c_double) :: gridy(*)
            real(kind=c_double) :: gridz(*)
            real(kind=c_double) :: ssw(*)
            real(kind=c_double) :: weight(*)
            integer(kind=c_int) :: atm(*)
            integer(kind=c_int) :: bin_locator(*)
            integer(kind=c_int) :: basf(*)
            integer(kind=c_int) :: primf(*)
            integer(kind=c_int) :: basf_counter(*)
            integer(kind=c_int) :: primf_counter(*)
            integer(kind=c_int) :: bin_counter(*)
        end subroutine get_gpu_grid_info
    end interface

    interface
        subroutine get_cpu_grid_info(gridx, gridy, gridz, ssw, weight, atm, basf, primf, basf_counter, primf_counter, bin_counter) &
            BIND(C, NAME="get_cpu_grid_info_")
            ! void get_cpu_grid_info_(double *gridx, double *gridy, double *gridz, double *ssw, double *weight, int *atm, int *basf, int *primf, int *basf_counter, int *primf_counter, int *bin_counter);
            use, intrinsic :: iso_c_binding, only : c_double, c_int
            real(kind=c_double) :: gridx(*)
            real(kind=c_double) :: gridy(*)
            real(kind=c_double) :: gridz(*)
            real(kind=c_double) :: ssw(*)
            real(kind=c_double) :: weight(*)
            integer(kind=c_int) :: atm(*)
            integer(kind=c_int) :: basf(*)
            integer(kind=c_int) :: primf(*)
            integer(kind=c_int) :: basf_counter(*)
            integer(kind=c_int) :: primf_counter(*)
            integer(kind=c_int) :: bin_counter(*)
        end subroutine get_cpu_grid_info
    end interface

    ! ------------------------------------------------------------------------
    ! to be sorted
    ! ------------------------------------------------------------------------

    interface
        subroutine gpu_aoint(leastIntegralCutoff, maxIntegralCutoff, intNum, intFileName) &
            BIND(C, NAME="gpu_aoint_")
            ! void gpu_aoint_(QUICKDouble* leastIntegralCutoff, QUICKDouble* maxIntegralCutoff, int* intNum, char* intFileName)
            use, intrinsic :: iso_c_binding, only : c_double, c_int, c_char
            real(kind=c_double) :: leastIntegralCutoff
            real(kind=c_double) :: maxIntegralCutoff
            integer(kind=c_int) :: intNum
            character(kind=c_char) :: intFileName(*)
        end subroutine gpu_aoint
    end interface
end module quick_c_interface