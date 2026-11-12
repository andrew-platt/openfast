@call "C:\Program Files (x86)\Intel\oneAPI\setvars-vcvarsall.bat" %VS_VER%

for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\compiler\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
@call "C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\env\vars.bat"

set "OverallErrorLevel=0"

echo ==== ifx version
ifx /QV

echo ==== Does RANDOM_NUMBER give the same stream when an array is filled in pieces?
ifx /O2 .github\issue3527\rn_chunk.f90 /exe:rn_chunk.exe
if %ERRORLEVEL% NEQ 0 set "OverallErrorLevel=1"
rn_chunk.exe
if %ERRORLEVEL% NEQ 0 set "OverallErrorLevel=1"

echo ==== TurbSim smoke runs
for %%E in (build\bin\TurbSim*.exe) do (
   echo ---- %%E RanLux
   %%E .github\issue3527\smoke.inp
   if errorlevel 1 set "OverallErrorLevel=1"
   echo ---- %%E intrinsic RNG
   %%E .github\issue3527\smoke_intrinsic.inp
   if errorlevel 1 set "OverallErrorLevel=1"
)
dir .github\issue3527

exit /b %OverallErrorLevel%
