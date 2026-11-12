@call "C:\Program Files (x86)\Intel\oneAPI\setvars-vcvarsall.bat" %VS_VER%

for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\compiler\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
@call "C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\env\vars.bat"

@REM Make the script that generates the git version description ignore dirty
@REM since building the Visual Studio projects modifies files
powershell -command "(Get-Content -Path '.\vs-build\CreateGitVersion.bat') -replace '--dirty', '' | Set-Content -Path '.\vs-build\CreateGitVersion.bat'"

@REM The Registry project is a dependency in the solution and may fail to build; it is
@REM not needed for TurbSim, so success is judged by whether the executables exist.
echo "Build TurbSim (Release|x64)"
devenv vs-build/OpenFAST.sln /Build "Release|x64" /Project TurbSim

echo "Build TurbSim (OpenMP_Release|x64)"
devenv vs-build/OpenFAST.sln /Build "OpenMP_Release|x64" /Project TurbSim

cd /d build\bin || exit /b 1
if exist TurbSim_Release.exe ren TurbSim_Release.exe TurbSim.exe
if exist TurbSim_OpenMP_Release.exe ren TurbSim_OpenMP_Release.exe TurbSim_OpenMP.exe

echo "List executables in build\bin"
dir

if not exist TurbSim.exe exit /b 1
if not exist TurbSim_OpenMP.exe exit /b 1
exit /b 0
