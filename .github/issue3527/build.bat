@call "C:\Program Files (x86)\Intel\oneAPI\setvars-vcvarsall.bat" %VS_VER%

for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\compiler\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
@call "C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\env\vars.bat"

@REM Make the script that generates the git version description ignore dirty
@REM since building the Visual Studio projects modifies files
powershell -command "(Get-Content -Path '.\vs-build\CreateGitVersion.bat') -replace '--dirty', '' | Set-Content -Path '.\vs-build\CreateGitVersion.bat'"

set "OverallErrorLevel=0"

echo "Build TurbSim (Release|x64)"
devenv vs-build/OpenFAST.sln /Build "Release|x64" /Project TurbSim
if %ERRORLEVEL% NEQ 0 set "OverallErrorLevel=1"

echo "Build TurbSim (OpenMP_Release|x64)"
devenv vs-build/OpenFAST.sln /Build "OpenMP_Release|x64" /Project TurbSim
if %ERRORLEVEL% NEQ 0 set "OverallErrorLevel=1"

echo "List executables in build\bin"
dir build\bin

exit /b %OverallErrorLevel%
