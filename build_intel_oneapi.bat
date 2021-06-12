SET PATH=C:\Windows\system32;C:\Windows;C:\Windows\System32\Wbem;C:\Windows\System32\WindowsPowerShell\v1.0\;C:\Windows\System32\OpenSSH\;C:\bin\Ninja;C:\bin\CMake\bin

call T:\Intel\setvars_no_vc.bat

echo %PATH%

cd QUICK
del /S /Q build
mkdir build
cd build


cmake .. -G Ninja -DCOMPILER=INTEL ^
    -DMPI=TRUE -DCUDA=TRUE -DQUICK_USER_ARCH=tesla -DMKL_STATIC=TRUE ^
    -DOpenMP_MKL_EXTRA_LIBRARIES=libiomp5md ^
    --trace --log-level=TRACE --debug-find --debug-trycompile 2> ../../LOG_10_MKL_EXTRA.log

REM    -DOpenMP_CXX_LIB_NAMES=libiomp5md -DOpenMP_C_LIB_NAMES=libiomp5md -DOpenMP_Fortran_LIB_NAMES=libiomp5md ^
REM    -DOpenMP_CXX_LIB_NAMES=libiomp5md -DOpenMP_C_LIB_NAMES=libiomp5md -DOpenMP_Fortran_LIB_NAMES=libiomp5md ^


REM ===== cmake debug options =====
REM See https://cmake.org/cmake/help/latest/manual/cmake.1.html

REM ===== OPENMP =====
REM
REM Cmake3.20.3 FindOpenMP is broken for Intel on windows
REM - there doesn't seem to be a comman dline options required (except for the include)
REM - it seems to be required to explicitly link libiomp5md
REM - See https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl/link-line-advisor.html
REM
REM QUICK FindMKL doesn't use OpenMP_CXX_LIBRARIES for its MKL build tests (fixed there)

REM CMAKE and intel https://www.scivision.dev/intel-compiler-cmake-make-load-libraries/

REM The Intel C compiler on Windows is always called ICL
REM SEE https://community.intel.com/t5/Intel-C-Compiler/What-s-the-difference-between-icx-icl-and-icc-compilers/td-p/1224714
REM The Intel C++ compiler is also called ICL
REM See https://software.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/compiler-setup/using-the-command-line/using-compiler-options.html

REM Amber Windows build notes
REM http://ambermd.org/pmwiki/pmwiki.php/Main/CMake

REM ERRORS
REM Not found T:/Intel/compiler/latest/windows/redist/intel64/compiler/*
REM Is in     T:/Intel/compiler/latest/windows/redist/intel64_win/compiler/*

REM COOL TOOL: HRML script to compute proper linker line
REM T:\Intel\mkl\2021.2.0\documentation\en\common
