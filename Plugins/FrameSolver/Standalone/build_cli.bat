@echo off
REM Builds frame_cli.exe (stdin/stdout solver) for the OpenSees offline comparison harness.
REM Same FrameCore translation units as build.bat, but main.cpp -> frame_cli.cpp.
setlocal enabledelayedexpansion
set "ROOT=%~dp0.."
set "EIGEN=E:\project\UE_5.7\Engine\Source\ThirdParty\Eigen"
set "VSWHERE=C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe"

set "VSDIR="
if exist "%VSWHERE%" (
  for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -prerelease -latest -property installationPath 2^>nul`) do set "VSDIR=%%i"
  if "!VSDIR!"=="" for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -property installationPath 2^>nul`) do set "VSDIR=%%i"
)
if "!VSDIR!"=="" ( echo [build_cli] could not locate Visual Studio via vswhere. & exit /b 1 )
rem Put vswhere's own folder on PATH so vcvars64's internal bare `vswhere` call resolves quietly.
for %%d in ("%VSWHERE%") do set "PATH=%%~dpd;%PATH%"
call "!VSDIR!\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 ( echo [build_cli] vcvars64 failed & exit /b 1 )

pushd "%ROOT%"
if not exist "Standalone\obj_cli" mkdir "Standalone\obj_cli"

set "GITSHA=unknown"
for /f "usebackq tokens=*" %%g in (`git -C "%ROOT%" rev-parse --short HEAD 2^>nul`) do set "GITSHA=%%g"

cl /nologo /EHsc /std:c++17 /O2 /MD /utf-8 /DEIGEN_MPL2_ONLY /DFRAMECORE_BUILD_SHA=\"!GITSHA!\" ^
   /Fe:Standalone\frame_cli.exe ^
   /Fo:Standalone\obj_cli\ ^
   /I"%EIGEN%" ^
   /I"Source\FrameCore\Public" ^
   /I"Source\FrameCore\Private" ^
   Source\FrameCore\Private\Node.cpp ^
   Source\FrameCore\Private\Material.cpp ^
   Source\FrameCore\Private\Section.cpp ^
   Source\FrameCore\Private\Member.cpp ^
   Source\FrameCore\Private\FrameModel.cpp ^
   Source\FrameCore\Private\ElementStiffness.cpp ^
   Source\FrameCore\Private\BeamColumnElement.cpp ^
   Source\FrameCore\Private\MITC4ShellElement.cpp ^
   Source\FrameCore\Private\FrameSolver.cpp ^
   Source\FrameCore\Private\ElasticAllowable.cpp ^
   Source\FrameCore\Private\Grillage.cpp ^
   Source\FrameCore\Private\SelfWeight.cpp ^
   Source\FrameCore\Private\Combination.cpp ^
   Source\FrameCore\Private\InfluenceLine.cpp ^
   Source\FrameCore\Private\ModalAnalysis.cpp ^
   Source\FrameCore\Private\BucklingAnalysis.cpp ^
   Source\FrameCore\Private\ResponseSpectrum.cpp ^
   Source\FrameCore\Private\ModalDynamics.cpp ^
   Source\FrameCore\Private\Connectivity.cpp ^
   Source\FrameCore\Private\Collapse.cpp ^
   Source\FrameCore\Private\DynamicCollapse.cpp ^
   Source\FrameCore\Private\Reanalysis.cpp ^
   Source\FrameCore\Private\PDeltaAnalysis.cpp ^
   Source\FrameCore\Private\CorotationalAnalysis.cpp ^
   Source\FrameCore\Private\TensionOnly.cpp ^
   Source\FrameCore\Private\SizeOpt.cpp ^
   Source\FrameCore\Private\Topology.cpp ^
   Standalone\frame_cli_core.cpp ^
   Standalone\frame_cli.cpp
if errorlevel 1 ( echo [build_cli] COMPILE FAILED & popd & exit /b 1 )

echo [build_cli] OK -^> Standalone\frame_cli.exe
popd
exit /b 0
