REM Best called with
REM CD D:\QuantumChemistry\Quick\QUICK
REM CMD /K build_intel_oneapi.bat
REM EXIT

SET PATH=C:\Windows\system32;C:\Windows;C:\Windows\System32\Wbem;C:\Windows\System32\WindowsPowerShell\v1.0\;C:\Windows\System32\OpenSSH\;C:\bin\Ninja;C:\bin\CMake\bin

call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Auxiliary\Build\vcvars64.bat"
call "T:\Intel\setvars_no_vc.bat"

echo %PATH%

del /S /Q build
mkdir build
cd build

SET CC=cl
SET CXX=cl
SET FC=ifort

cmake .. -G Ninja -DCOMPILER=MANUAL ^
    -DCUDA_NVCC_FLAGS="--ptxas-options=-O0" ^
    -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=TRUE -DBUILD_SHARED_LIBS=TRUE ^
    -DMPI=TRUE -DCUDA=TRUE -DQUICK_USER_ARCH=turing -DMKL_STATIC=TRUE ^
    -DOpenMP_MKL_EXTRA_LIBRARIES=libiomp5md ^
    --trace --log-level=TRACE --debug-find --debug-trycompile 2> ../../LOG_14_VS2019_nvccO0.log

REM C:\bin\ninja\ninja.exe -k 0

REM This has an effect, but teh effect is that it does not configure any more
REM     -DCUDA_NVCC_FLAGS="--ptxas-options=-O0" ^


GOTO :EOF

REM    -DOpenMP_CXX_LIB_NAMES=libiomp5md -DOpenMP_C_LIB_NAMES=libiomp5md -DOpenMP_Fortran_LIB_NAMES=libiomp5md ^
REM    -DOpenMP_CXX_LIB_NAMES=libiomp5md -DOpenMP_C_LIB_NAMES=libiomp5md -DOpenMP_Fortran_LIB_NAMES=libiomp5md ^

===== Export all symbols from library in cmake build =====

See https://blog.kitware.com/create-dlls-on-windows-without-declspec-using-new-cmake-export-all-feature/
ADD -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=TRUE -DBUILD_SHARED_LIBS=TRUE ^

===== Dump all symbols from a library =====

dumpbin /ALL quick.lib > quick.lib.dmp

===== DEBUG Ninja build steps =====

=== Build single file with "-j 1 -v -d keeprsp" reveals command and rsp file
C:\bin\Ninja\ninja.exe -v -d keeprsp src\quick.lib
C:\bin\Ninja\ninja.exe -v -d keeprsp src\quick.exe

/VERBOSE:LIB

===== Maple source files - in case there are bugs in the generated C code =====
https://www.e-cam2020.eu:10443/esl/libxc/-/tree/4.0.1/maple

===== cmake debug options =====
See https://cmake.org/cmake/help/latest/manual/cmake.1.html

===== OPENMP =====

Cmake3.20.3 FindOpenMP is broken for Intel on windows
- there doesn't seem to be a comman dline options required (except for the include)
- it seems to be required to explicitly link libiomp5md
- See https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl/link-line-advisor.html

QUICK FindMKL doesn't use OpenMP_CXX_LIBRARIES for its MKL build tests (fixed there)
CMAKE and intel https://www.scivision.dev/intel-compiler-cmake-make-load-libraries/
The Intel C compiler on Windows is always called ICL
SEE https://community.intel.com/t5/Intel-C-Compiler/What-s-the-difference-between-icx-icl-and-icc-compilers/td-p/1224714
The Intel C++ compiler is also called ICL
See https://software.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/compiler-setup/using-the-command-line/using-compiler-options.html
Amber Windows build notes
http://ambermd.org/pmwiki/pmwiki.php/Main/CMake
ERRORS
Not found T:/Intel/compiler/latest/windows/redist/intel64/compiler/*
Is in     T:/Intel/compiler/latest/windows/redist/intel64_win/compiler/*
COOL TOOL: HRML script to compute proper linker line
T:\Intel\mkl\2021.2.0\documentation\en\common
