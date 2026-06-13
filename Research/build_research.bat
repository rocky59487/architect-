@echo off
REM Research experiments build — NOT part of the engine gate (see Research\README.md).
REM Compiles FrameCore sources once into Research\obj_core\, then links each
REM experiment exe into Research\bin\. Engine sources are read-only inputs.
REM
REM Usage:  build_research.bat [-skipcore] [exp_name]
REM   -skipcore  reuse existing Research\obj_core\*.obj
REM   exp_name   build only that experiment (e.g. exp_sparse_buckling)
setlocal enabledelayedexpansion
set "ROOT=%~dp0.."
set "FS=Plugins\FrameSolver"
set "EIGEN=%EIGEN_DIR%"
if "%EIGEN%"=="" if defined UE_ENGINE_ROOT set "EIGEN=%UE_ENGINE_ROOT%\Engine\Source\ThirdParty\Eigen"
if "%EIGEN%"=="" set "EIGEN=%~dp0..\..\UE_5.7\Engine\Source\ThirdParty\Eigen"
if not exist "%EIGEN%\Eigen" ( echo [research] Eigen include root not found: "%EIGEN%" & exit /b 1 )
set "VSWHERE=C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe"

set "SKIPCORE=0"
set "ONLY="
:parseargs
if "%~1"=="" goto doneargs
if /I "%~1"=="-skipcore" ( set "SKIPCORE=1" ) else ( set "ONLY=%~1" )
shift
goto parseargs
:doneargs

set "VSDIR="
if exist "%VSWHERE%" (
  for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -prerelease -latest -property installationPath 2^>nul`) do set "VSDIR=%%i"
  if "!VSDIR!"=="" for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -property installationPath 2^>nul`) do set "VSDIR=%%i"
)
if "!VSDIR!"=="" ( echo [research] could not locate Visual Studio via vswhere. & exit /b 1 )
for %%d in ("%VSWHERE%") do set "PATH=%%~dpd;!PATH!"
call "!VSDIR!\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 ( echo [research] vcvars64 failed & exit /b 1 )

pushd "%ROOT%"
if not exist "Research\obj_core" mkdir "Research\obj_core"
if not exist "Research\obj_exp"  mkdir "Research\obj_exp"
if not exist "Research\bin"      mkdir "Research\bin"
if not exist "Research\out"      mkdir "Research\out"

set CFLAGS=/nologo /EHsc /std:c++17 /O2 /MD /utf-8 /DEIGEN_MPL2_ONLY /I"%EIGEN%" /I"%FS%\Source\FrameCore\Public" /I"%FS%\Source\FrameCore\Private" /I"Research\common"

if "%SKIPCORE%"=="1" goto exps
echo [research] compiling FrameCore objs...
cl %CFLAGS% /MP /c /Fo:Research\obj_core\ ^
   %FS%\Source\FrameCore\Private\Node.cpp ^
   %FS%\Source\FrameCore\Private\Material.cpp ^
   %FS%\Source\FrameCore\Private\Section.cpp ^
   %FS%\Source\FrameCore\Private\Member.cpp ^
   %FS%\Source\FrameCore\Private\FrameModel.cpp ^
   %FS%\Source\FrameCore\Private\ElementStiffness.cpp ^
   %FS%\Source\FrameCore\Private\BeamColumnElement.cpp ^
   %FS%\Source\FrameCore\Private\MITC4ShellElement.cpp ^
   %FS%\Source\FrameCore\Private\FrameSolver.cpp ^
   %FS%\Source\FrameCore\Private\ElasticAllowable.cpp ^
   %FS%\Source\FrameCore\Private\Grillage.cpp ^
   %FS%\Source\FrameCore\Private\SelfWeight.cpp ^
   %FS%\Source\FrameCore\Private\Combination.cpp ^
   %FS%\Source\FrameCore\Private\InfluenceLine.cpp ^
   %FS%\Source\FrameCore\Private\ModalAnalysis.cpp ^
   %FS%\Source\FrameCore\Private\BucklingAnalysis.cpp ^
   %FS%\Source\FrameCore\Private\ResponseSpectrum.cpp ^
   %FS%\Source\FrameCore\Private\ModalDynamics.cpp ^
   %FS%\Source\FrameCore\Private\Connectivity.cpp ^
   %FS%\Source\FrameCore\Private\Collapse.cpp
if errorlevel 1 ( echo [research] CORE COMPILE FAILED & popd & exit /b 1 )

:exps
set "EXPS=WS_N_incremental\exp_incremental_refactor WS_B_solver\exp_sparse_buckling WS_B_solver\exp_million_dof WS_B_solver\exp_solver_compare WS_B_solver\exp_matrix_free_operator WS_B_solver\exp_matrix_free_pcg WS_B_solver\exp_bsr6_matvec WS_B_solver\exp_framecore_bsr6_matvec WS_B_solver\exp_threaded_element_apply WS_B_solver\exp_pcg_recycling WS_B_solver\exp_line_schwarz_pcg WS_B_solver\exp_parallel_pcg WS_B_solver\exp_full_apply_oracle WS_B_solver\exp_mechanism_guard_oracle WS_C_pdelta\exp_pdelta_convergence WS_D_tensiononly\exp_tension_only WS_H_sizeopt\exp_size_opt WS_I_beso\exp_beso_truss WS_N_incremental\exp_dynamic_inherit"

for %%E in (%EXPS%) do (
  set "SRC=Research\%%E.cpp"
  for %%F in ("%%E") do set "BASE=%%~nF"
  if exist "!SRC!" (
    set "BUILDIT=1"
    if not "%ONLY%"=="" if /I not "!BASE!"=="%ONLY%" set "BUILDIT=0"
    if "!BUILDIT!"=="1" (
      echo [research] building !BASE! ...
      cl %CFLAGS% /Fo:Research\obj_exp\ /Fe:Research\bin\!BASE!.exe "!SRC!" Research\obj_core\*.obj
      if errorlevel 1 ( echo [research] !BASE! COMPILE FAILED & popd & exit /b 1 )
    )
  )
)

echo [research] OK
popd
exit /b 0
