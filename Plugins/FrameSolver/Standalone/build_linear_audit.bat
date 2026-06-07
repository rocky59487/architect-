@echo off
REM Builds and runs the post-F17/F25 deep audit table for FrameCore.
setlocal enabledelayedexpansion
set "ROOT=%~dp0.."
set "EIGEN=E:\project\UE_5.7\Engine\Source\ThirdParty\Eigen"
set "VSWHERE=C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe"

set "VSDIR="
if exist "%VSWHERE%" (
  for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -prerelease -latest -property installationPath 2^>nul`) do set "VSDIR=%%i"
  if "!VSDIR!"=="" for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -property installationPath 2^>nul`) do set "VSDIR=%%i"
)
if "!VSDIR!"=="" (
  echo [build] could not locate Visual Studio via vswhere. Open an "x64 Native Tools" prompt and run cl directly.
  exit /b 1
)
call "!VSDIR!\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 ( echo [build] vcvars64 failed & exit /b 1 )

pushd "%ROOT%"
if not exist "Standalone\obj_linear_audit" mkdir "Standalone\obj_linear_audit"

cl /nologo /EHsc /std:c++17 /O2 /MD /utf-8 /DEIGEN_MPL2_ONLY ^
   /Fe:Standalone\linear_deep_audit.exe ^
   /Fo:Standalone\obj_linear_audit\ ^
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
   Standalone\linear_deep_audit.cpp
if errorlevel 1 ( echo [build] COMPILE FAILED & popd & exit /b 1 )

echo [build] OK -^> Standalone\linear_deep_audit.exe
echo.
"Standalone\linear_deep_audit.exe"
set "RC=%errorlevel%"
popd
exit /b %RC%
