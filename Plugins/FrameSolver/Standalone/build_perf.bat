@echo off
REM Builds frame_perf.exe, an in-process FrameCore performance benchmark.
setlocal enabledelayedexpansion
set "ROOT=%~dp0.."
set "EIGEN=E:\project\UE_5.7\Engine\Source\ThirdParty\Eigen"
set "VSWHERE=C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe"

set "VSDIR="
if exist "%VSWHERE%" (
  for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -prerelease -latest -property installationPath 2^>nul`) do set "VSDIR=%%i"
  if "!VSDIR!"=="" for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -property installationPath 2^>nul`) do set "VSDIR=%%i"
)
if "!VSDIR!"=="" ( echo [build_perf] could not locate Visual Studio via vswhere. & exit /b 1 )
rem Put vswhere's own folder on PATH so vcvars64's internal bare `vswhere` call resolves quietly.
for %%d in ("%VSWHERE%") do set "PATH=%%~dpd;%PATH%"
call "!VSDIR!\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 ( echo [build_perf] vcvars64 failed & exit /b 1 )

pushd "%ROOT%"
if not exist "Standalone\obj_perf" mkdir "Standalone\obj_perf"

cl /nologo /EHsc /std:c++17 /O2 /MD /utf-8 /DEIGEN_MPL2_ONLY ^
   /Fe:Standalone\frame_perf.exe ^
   /Fo:Standalone\obj_perf\ ^
   /I"%EIGEN%" ^
   /I"Source\FrameCore\Public" ^
   /I"Source\FrameCore\Private" ^
   Source\FrameCore\Private\Section.cpp ^
   Source\FrameCore\Private\Member.cpp ^
   Source\FrameCore\Private\FrameModel.cpp ^
   Source\FrameCore\Private\ElementStiffness.cpp ^
   Source\FrameCore\Private\BeamColumnElement.cpp ^
   Source\FrameCore\Private\MITC4ShellElement.cpp ^
   Source\FrameCore\Private\FrameSolver.cpp ^
   Source\FrameCore\Private\ElasticAllowable.cpp ^
   Source\FrameCore\Private\Grillage.cpp ^
   Standalone\frame_perf.cpp
if errorlevel 1 ( echo [build_perf] COMPILE FAILED & popd & exit /b 1 )

echo [build_perf] OK -^> Standalone\frame_perf.exe
popd
exit /b 0
