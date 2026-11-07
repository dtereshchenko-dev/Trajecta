[CmdletBinding()]
param(
    [string]$Repository,
    [string]$BuildDirectory,
    [string]$QtDirectory = 'C:\Qt\6.9.3\msvc2022_64',
    [string]$IfwDirectory = 'C:\Qt\Tools\QtInstallerFramework\4.11',
    [string]$VcpkgPrefix = 'C:\vcpkg\installed\x64-windows',
    [switch]$StageOnly,
    [switch]$PrepareUpdateRepository
)
$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
if ($PrepareUpdateRepository -and !$StageOnly) {
    throw 'PrepareUpdateRepository requires StageOnly so the generated repository is retained.'
}
$root = [IO.Path]::GetFullPath($PSScriptRoot)
if (!$Repository) { $Repository = Split-Path -Parent $root }
if (!$BuildDirectory) { $BuildDirectory = Join-Path $Repository 'build\Desktop_Qt_6_9_3_MSVC2022_64bit-Release-Trajecta' }

function Require-Path([string]$Path) {
    if (!(Test-Path -LiteralPath $Path)) { throw "Required input missing: $Path" }
}
function Invoke-Checked([string]$Executable, [string[]]$Arguments) {
    & $Executable @Arguments
    if ($LASTEXITCODE -ne 0) { throw "$Executable failed with exit code $LASTEXITCODE" }
}
function Copy-Required([string]$Source, [string]$Destination) {
    Require-Path $Source
    $null = New-Item -ItemType Directory -Force -Path (Split-Path -Parent $Destination)
    Copy-Item -LiteralPath $Source -Destination $Destination -Force
}
function Set-InstallerVersion([string]$Path, [string]$Version, [string]$RootElement) {
    $document = [System.Xml.XmlDocument]::new()
    $document.PreserveWhitespace = $true
    $document.Load($Path)
    $node = $document.SelectSingleNode("/$RootElement/Version")
    if (!$node -or $node.InnerText -ne '@TRAJECTA_VERSION@') {
        throw "Version placeholder missing from $Path"
    }
    $node.InnerText = $Version
    if ($RootElement -eq 'Installer') {
        $title = $document.SelectSingleNode('/Installer/Title')
        if (!$title -or !$title.InnerText.Contains('@TRAJECTA_VERSION@')) {
            throw "Title version placeholder missing from $Path"
        }
        $title.InnerText = $title.InnerText.Replace('@TRAJECTA_VERSION@', $Version)
    }
    $document.Save($Path)
}
function Assert-Local([string]$Path) {
    $full = [IO.Path]::GetFullPath($Path)
    if (!$full.StartsWith($root + '\', [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing operation outside installer directory: $full"
    }
    return $full
}

function Remove-StagedPath([string]$Payload, [string]$Relative) {
    $rootPath = (Resolve-Path -LiteralPath $Payload).Path.TrimEnd('\') + '\'
    $path = Join-Path $rootPath $Relative
    if (!(Test-Path -LiteralPath $path)) { return }
    $resolved = (Resolve-Path -LiteralPath $path).Path
    if (!$resolved.StartsWith($rootPath, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to remove item outside staged payload: $resolved"
    }
    $item = Get-Item -LiteralPath $resolved -Force
    if (($item.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) {
        throw "Refusing linked staged path: $resolved"
    }
    if ($item.PSIsContainer) {
        if (Get-ChildItem -LiteralPath $resolved -Recurse -Force | Where-Object {
            ($_.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0
        } | Select-Object -First 1) { throw "Refusing linked entry under staged path: $resolved" }
        Remove-Item -LiteralPath $resolved -Recurse -Force
    } else {
        Remove-Item -LiteralPath $resolved -Force
    }
}

function Remove-StagedStyle([string]$Payload, [string]$Name) {
    foreach ($relative in @("qml\QtQuick\Controls\$Name", "Qt6QuickControls2$Name.dll",
            "Qt6QuickControls2${Name}StyleImpl.dll")) {
        Remove-StagedPath -Payload $Payload -Relative $relative
    }
}

function Get-RuntimeDependencyScript {
    return @'
cmake_minimum_required(VERSION 3.21)
foreach(required PAYLOAD DUMPBIN REPORT)
    if(NOT DEFINED ${required})
        message(FATAL_ERROR "Missing ${required}")
    endif()
endforeach()
set(CMAKE_GET_RUNTIME_DEPENDENCIES_PLATFORM "windows+pe")
set(CMAKE_GET_RUNTIME_DEPENDENCIES_TOOL "dumpbin")
set(CMAKE_GET_RUNTIME_DEPENDENCIES_COMMAND "${DUMPBIN}")
file(GLOB executables "${PAYLOAD}/*.exe")
file(GLOB_RECURSE libraries "${PAYLOAD}/*.dll")
file(GET_RUNTIME_DEPENDENCIES
    EXECUTABLES ${executables}
    LIBRARIES ${libraries}
    DIRECTORIES "${PAYLOAD}" ${SEARCH_DIRECTORIES}
    RESOLVED_DEPENDENCIES_VAR resolved
    UNRESOLVED_DEPENDENCIES_VAR unresolved
    CONFLICTING_DEPENDENCIES_PREFIX conflicts
    PRE_EXCLUDE_REGEXES "^[Aa][Pp][Ii]-[Mm][Ss]-" "^[Ee][Xx][Tt]-[Mm][Ss]-"
    POST_EXCLUDE_REGEXES ".*[/\\][Ww][Ii][Nn][Dd][Oo][Ww][Ss][/\\]([Ss][Yy][Ss][Tt][Ee][Mm]32|[Ss][Yy][Ss][Ww][Oo][Ww]64)[/\\].*")
if(unresolved)
    message(FATAL_ERROR "Unresolved dependencies: ${unresolved}")
endif()
foreach(name IN LISTS conflicts_FILENAMES)
    set(paths "${conflicts_${name}}")
    list(GET paths 0 first)
    file(SHA256 "${first}" expected)
    foreach(path IN LISTS paths)
        file(SHA256 "${path}" actual)
        if(NOT actual STREQUAL expected)
            message(FATAL_ERROR "Conflicting versions of ${name}: ${paths}")
        endif()
    endforeach()
    list(APPEND resolved "${first}")
endforeach()
file(REAL_PATH "${PAYLOAD}" payload_real)
foreach(dependency IN LISTS resolved)
    if(COPY_DEPENDENCIES)
        get_filename_component(name "${dependency}" NAME)
        if(EXISTS "${PAYLOAD}/${name}")
            file(REAL_PATH "${PAYLOAD}/${name}" destination)
        else()
            set(destination "${payload_real}/${name}")
        endif()
        file(REAL_PATH "${dependency}" source)
        if(NOT source STREQUAL destination)
            file(COPY_FILE "${dependency}" "${PAYLOAD}/${name}" ONLY_IF_DIFFERENT)
        endif()
    else()
        file(REAL_PATH "${dependency}" dependency_real)
        cmake_path(IS_PREFIX payload_real "${dependency_real}" NORMALIZE local)
        if(NOT local)
            message(FATAL_ERROR "Dependency outside payload: ${dependency}")
        endif()
    endif()
endforeach()
list(SORT resolved)
string(REPLACE ";" "\n" resolved "${resolved}")
file(WRITE "${REPORT}" "${resolved}\n")
message(STATUS "Runtime dependency closure verified")
'@
}

function Invoke-PayloadVerification {
    [CmdletBinding()]
    param(
        [Parameter(Mandatory)][string]$Payload,
        [Parameter(Mandatory)][string]$Repository,
        [Parameter(Mandatory)][string]$BuildDirectory,
        [Parameter(Mandatory)][string]$QtDirectory,
        [Parameter(Mandatory)][string]$CMake,
        [Parameter(Mandatory)][string]$Dumpbin,
        [Parameter(Mandatory)][string]$ReportDirectory,
        [Parameter(Mandatory)][string]$DependencyScript
    )
    $ErrorActionPreference = 'Stop'
    Set-StrictMode -Version Latest
    $Payload = [IO.Path]::GetFullPath($Payload)
    $ReportDirectory = [IO.Path]::GetFullPath($ReportDirectory)
    $root = [IO.Path]::GetFullPath($PSScriptRoot) + '\'
    foreach ($path in @($Payload, $ReportDirectory)) {
        if (!$path.StartsWith($root, [StringComparison]::OrdinalIgnoreCase)) { throw "Not inside installer workspace: $path" }
    }
    $null = New-Item -ItemType Directory -Path $ReportDirectory -Force
    $required = @('Trajecta.exe', 'TrajectaCLI.exe', 'StrewnCLI.exe', 'OrbitalCore.dll', 'StrewnCore.dll',
        'gdal.dll', 'libcurl.dll', 'Qt6Core.dll', 'Qt6Network.dll', 'platforms\qwindows.dll', 'tls\qschannelbackend.dll',
        'msvcp140.dll', 'vcruntime140.dll', 'vcruntime140_1.dll', 'vcomp140.dll',
        'share\proj\proj.db', 'share\gdal\grib2_table_4_2_0_0.csv',
        'linux_p1550p2650.440', 'OBSCODE.dat', 'RADCODE.dat',
        'debias\efcc18.dat', 'debias\fcct14.dat', 'weights\vfcc17.rules',
        'earth_orientation\finals2000A.all', 'earth_orientation\eopc04.62-now', 'earth_orientation\eopc04.cnc',
        'earth_orientation\ET-UT.dat')
    foreach ($file in $required) {
        if (!(Test-Path -LiteralPath "$Payload\$file" -PathType Leaf) -or (Get-Item -LiteralPath "$Payload\$file").Length -eq 0) {
            throw "Missing/empty payload file: $file"
        }
    }
    foreach ($relative in @('payload-manifest.json', 'qmltooling',
            'qml\QtQuick\Controls\Fusion', 'qml\QtQuick\Controls\Imagine',
            'qml\QtQuick\Controls\Material', 'qml\QtQuick\Controls\Universal')) {
        if (Test-Path -LiteralPath "$Payload\$relative") {
            throw "Unexpected release payload item: $relative"
        }
    }
    foreach ($file in @('Trajecta.exe', 'TrajectaCLI.exe', 'StrewnCLI.exe', 'OrbitalCore.dll', 'StrewnCore.dll')) {
        $folder = switch ($file) {
            'Trajecta.exe' { 'TrajectaUI' }
            'TrajectaCLI.exe' { 'TrajectaCLI' }
            'OrbitalCore.dll' { 'OrbitalCore' }
            default { 'StrewnField' }
        }
        if ((Get-FileHash "$BuildDirectory\$folder\$file").Hash -ne (Get-FileHash "$Payload\$file").Hash) {
            throw "Stale binary in payload: $file"
        }
    }

    $variables = @('PATH', 'GDAL_DATA', 'PROJ_DATA', 'PROJ_LIB', 'GDAL_DRIVER_PATH', 'QT_PLUGIN_PATH',
        'QML_IMPORT_PATH', 'QML2_IMPORT_PATH', 'QT_QPA_PLATFORM', 'QT_QUICK_BACKEND')
    $saved = @{}
    foreach ($name in $variables) { $saved[$name] = [Environment]::GetEnvironmentVariable($name, 'Process') }
    $test = "$Payload\StrewnCoreTests.exe"
    $offscreen = "$Payload\platforms\qoffscreen.dll"
    if (Test-Path -LiteralPath $offscreen) { throw 'Verification platform plugin already exists in payload' }
    if (Test-Path -LiteralPath $test) { throw 'Verification test executable already exists in payload' }
    $testDlls = @(Get-ChildItem -LiteralPath "$BuildDirectory\StrewnField" -Filter 'gtest*.dll')
    foreach ($dll in $testDlls) {
        if (Test-Path -LiteralPath "$Payload\$($dll.Name)") { throw "Unexpected test DLL in payload: $($dll.Name)" }
    }
    try {
        foreach ($name in $variables) { Remove-Item -LiteralPath "Env:$name" -ErrorAction SilentlyContinue }
        $env:PATH = "$env:SystemRoot\System32;$env:SystemRoot"
        $env:QT_QPA_PLATFORM = 'offscreen'
        $env:QT_QUICK_BACKEND = 'software'
        & $CMake "-DPAYLOAD=$Payload" "-DDUMPBIN=$Dumpbin" "-DREPORT=$ReportDirectory\dependencies-local.txt" `
            -P $DependencyScript
        if ($LASTEXITCODE -ne 0) { throw 'Runtime dependency verification failed' }
        Push-Location $ReportDirectory
        try {
            foreach ($name in @('TrajectaCLI', 'StrewnCLI')) {
                & "$Payload\$name.exe" --help *> "$ReportDirectory\$name-help.log"
                if ($LASTEXITCODE -ne 0) { throw "$name could not start in clean environment" }
            }
            Copy-Item -LiteralPath "$BuildDirectory\StrewnField\StrewnCoreTests.exe" -Destination $test
            $testDlls | Copy-Item -Destination $Payload
            Copy-Item -LiteralPath "$QtDirectory\plugins\platforms\qoffscreen.dll" -Destination $offscreen
            & $test '--gtest_filter=WeatherTest.*:JobTest.MapsAreLibraryOwnedAndReplayable' `
                "--gtest_output=json:$ReportDirectory/raster-tests.json" *> "$ReportDirectory\raster-tests.log"
            if ($LASTEXITCODE -ne 0) { throw "Raster tests failed; see $ReportDirectory\raster-tests.log" }
            $job = Get-Content -LiteralPath "$Repository\StrewnField\examples\cx1-offline.json" -Raw | ConvertFrom-Json
            $job.monte_carlo.runs = 3
            $job.monte_carlo.workers = 2
            $job.output_directory = "$ReportDirectory\smoke-run"
            $job | ConvertTo-Json -Depth 20 | Set-Content -LiteralPath "$ReportDirectory\smoke-job.json" -Encoding UTF8
            & "$Payload\StrewnCLI.exe" "$ReportDirectory\smoke-job.json" *> "$ReportDirectory\strewn-smoke.log"
            if ($LASTEXITCODE -ne 0) { throw "Strewn smoke run failed; see $ReportDirectory\strewn-smoke.log" }
            $maps = Get-Content -LiteralPath "$Repository\StrewnField\examples\maps-fixed.json" -Raw | ConvertFrom-Json
            $maps | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath "$ReportDirectory\smoke-maps.json" -Encoding UTF8
            & "$Payload\StrewnCLI.exe" --render "$ReportDirectory\smoke-run" "$ReportDirectory\smoke-maps" `
                "$ReportDirectory\smoke-maps.json" *> "$ReportDirectory\strewn-render.log"
            if ($LASTEXITCODE -ne 0) { throw "Strewn render failed; see $ReportDirectory\strewn-render.log" }
            foreach ($name in @('field.html', 'density.tif', 'contours.geojson')) {
                if (!(Test-Path -LiteralPath "$ReportDirectory\smoke-maps\$name")) { throw "Missing generated artifact: $name" }
            }
            $process = Start-Process -FilePath "$Payload\Trajecta.exe" -WorkingDirectory $ReportDirectory `
                -WindowStyle Hidden -PassThru -RedirectStandardOutput "$ReportDirectory\ui-stdout.log" `
                -RedirectStandardError "$ReportDirectory\ui-stderr.log"
            $null = $process.Handle
            try {
                if ($process.WaitForExit(8000)) {
                    Get-Content -LiteralPath "$ReportDirectory\ui-stderr.log" | Write-Host
                    Get-Content -LiteralPath "$ReportDirectory\ui-stdout.log" | Write-Host
                    throw "UI exited during startup with code $($process.ExitCode)"
                }
            } finally {
                if (!$process.HasExited) { Stop-Process -Id $process.Id -Force; $process.WaitForExit() }
            }
            if ((Get-Content -LiteralPath "$ReportDirectory\ui-stderr.log" -Raw) -match 'failed to load component|module .* is not installed|Cannot load library') {
                throw "UI resource error; see $ReportDirectory\ui-stderr.log"
            }
        } finally { Pop-Location }
    } finally {
        if (Test-Path -LiteralPath $test) { Remove-Item -LiteralPath $test -Force }
        if (Test-Path -LiteralPath $offscreen) { Remove-Item -LiteralPath $offscreen -Force }
        foreach ($dll in $testDlls) {
            $path = "$Payload\$($dll.Name)"
            if (Test-Path -LiteralPath $path) { Remove-Item -LiteralPath $path -Force }
        }
        foreach ($name in $variables) {
            if ($null -eq $saved[$name]) { Remove-Item -LiteralPath "Env:$name" -ErrorAction SilentlyContinue }
            else { [Environment]::SetEnvironmentVariable($name, $saved[$name], 'Process') }
        }
    }
    Write-Host 'Payload verification passed'
}

$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
Require-Path $vswhere
$vs = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
if (!$vs) { throw 'MSVC x64 toolchain not found' }
& "$vs\Common7\Tools\Launch-VsDevShell.ps1" -Arch amd64 -HostArch amd64 -SkipAutomaticLocation
$env:VSLANG = '1049'
chcp.com 65001 | Out-Null
[Console]::InputEncoding = [Text.UTF8Encoding]::new($false)
[Console]::OutputEncoding = [Text.UTF8Encoding]::new($false)
$cmake = "$vs\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
$dumpbin = (Get-Command dumpbin.exe -ErrorAction Stop).Source
$redist = Join-Path $env:VCToolsRedistDir 'x64'
$binarycreator = Join-Path $IfwDirectory 'bin\binarycreator.exe'
$archivegen = Join-Path $IfwDirectory 'bin\archivegen.exe'
$repogen = Join-Path $IfwDirectory 'bin\repogen.exe'
foreach ($path in @($cmake, $dumpbin, "$QtDirectory\bin\windeployqt.exe", $binarycreator,
        $archivegen, $repogen,
        "$VcpkgPrefix\share\gdal\grib2_table_4_2_0_0.csv", "$VcpkgPrefix\share\proj\proj.db",
        "$root\config\config.xml", "$root\config\icon.ico")) { Require-Path $path }

$binaries = [ordered]@{
    'Trajecta.exe' = "$BuildDirectory\TrajectaUI\Trajecta.exe"
    'TrajectaCLI.exe' = "$BuildDirectory\TrajectaCLI\TrajectaCLI.exe"
    'StrewnCLI.exe' = "$BuildDirectory\StrewnField\StrewnCLI.exe"
    'OrbitalCore.dll' = "$BuildDirectory\OrbitalCore\OrbitalCore.dll"
    'StrewnCore.dll' = "$BuildDirectory\StrewnField\StrewnCore.dll"
}
$requiredData = @('linux_p1550p2650.440', 'OBSCODE.dat', 'RADCODE.dat',
    'debias\efcc18.dat', 'debias\fcct14.dat', 'weights\vfcc17.rules',
    'earth_orientation\finals2000A.all', 'earth_orientation\eopc04.62-now', 'earth_orientation\eopc04.cnc',
    'earth_orientation\ET-UT.dat')
foreach ($path in $binaries.Values) { Require-Path $path }
foreach ($file in $requiredData) { Require-Path "$Repository\data\$file" }
Require-Path "$BuildDirectory\StrewnField\StrewnCoreTests.exe"
$versionFile = "$BuildDirectory\TrajectaUI\trajecta-version.txt"
Require-Path $versionFile
$version = (Get-Content -Raw -LiteralPath $versionFile).Trim()
if ($version -notmatch '^\d+\.\d+\.\d+$') { throw "Invalid application version: $version" }
$executableVersion = (Get-Item -LiteralPath $binaries['Trajecta.exe']).VersionInfo.ProductVersion
if ($executableVersion -ne $version) {
    throw "Application version mismatch: build declares $version, Trajecta.exe reports $executableVersion. Rebuild before packaging."
}
Write-Host "Packaging Trajecta $version"

$stamp = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0, 8)
$stage = Assert-Local "$root\.staging\$stamp"
$package = "$stage\packages\com.trajecta.app"
$payload = "$package\data"
$null = New-Item -ItemType Directory -Path $payload -Force
try {
$dependencyScript = "$stage\runtime-dependencies.cmake"
Get-RuntimeDependencyScript | Set-Content -LiteralPath $dependencyScript -Encoding UTF8
Copy-Item -LiteralPath "$root\packages\com.trajecta.app\meta" -Destination "$package\meta" -Recurse
Copy-Item -LiteralPath "$root\config" -Destination "$stage\config" -Recurse
[IO.File]::AppendAllText("$stage\config\controller.qs", "`n" + [IO.File]::ReadAllText("$root\common\lifecycle.qs"))
Set-InstallerVersion -Path "$stage\config\config.xml" -Version $version -RootElement 'Installer'
Set-InstallerVersion -Path "$package\meta\package.xml" -Version $version -RootElement 'Package'

Write-Host '[1/6] Copying current applications and runtime data'
foreach ($item in $binaries.GetEnumerator()) { Copy-Required $item.Value "$payload\$($item.Key)" }
foreach ($file in $requiredData) { Copy-Required "$Repository\data\$file" "$payload\$file" }
foreach ($directory in @('docs', 'examples')) {
    Require-Path "$Repository\data\$directory"
    Copy-Item -LiteralPath "$Repository\data\$directory" -Destination "$payload\$directory" -Recurse
}
Copy-Required "$Repository\TrajectaUI\THIRD_PARTY.md" "$payload\THIRD_PARTY.md"
Copy-Required "$root\config\icon.ico" "$payload\icon.ico"
Copy-Required "$Repository\LICENSE" "$payload\licenses\GPL-3.0.txt"
Copy-Required "$Repository\licenses\LGPL-3.0.txt" "$payload\licenses\LGPL-3.0.txt"
Copy-Required "$Repository\TrajectaUI\resources\fonts\LICENSE-Inter.txt" "$payload\licenses\Inter-OFL.txt"
Copy-Required "$Repository\TrajectaUI\resources\fonts\LICENSE-JetBrainsMono.txt" "$payload\licenses\JetBrainsMono-OFL.txt"
Copy-Required "$Repository\OrbitalCore\third_party\erfa\LICENSE" "$payload\licenses\ERFA.txt"
foreach ($name in @('leaflet', 'plotly', 'lucide')) {
    Copy-Required "$Repository\StrewnField\output\web\vendor\$name\LICENSE" "$payload\licenses\$name.txt"
}
foreach ($name in @('gdal', 'proj')) {
    $null = New-Item -ItemType Directory -Path "$payload\share\$name" -Force
    Get-ChildItem -LiteralPath "$VcpkgPrefix\share\$name" -File |
        Where-Object { $_.Name -notmatch '\.cmake$|^vcpkg|^usage$' } |
        Copy-Item -Destination "$payload\share\$name"
}
foreach ($name in @('Microsoft.VC143.CRT', 'Microsoft.VC143.OpenMP')) {
    Require-Path "$redist\$name"
    Get-ChildItem -LiteralPath "$redist\$name" -Filter '*.dll' | Copy-Item -Destination $payload
}

Write-Host '[2/6] Deploying Qt/QML, including TLS and rendering support'
$env:PATH = "$QtDirectory\bin;" + $env:PATH
Invoke-Checked "$QtDirectory\bin\windeployqt.exe" @('--release', '--qmldir', "$Repository\TrajectaUI\qml",
    '--verbose', '0',
    '--no-translations', '--no-compiler-runtime', '--skip-plugin-types', 'qmltooling',
    '--openssl-root', $VcpkgPrefix, '--dir', $payload, "$payload\Trajecta.exe")

# The application explicitly selects Qt Quick Controls Basic.
Remove-StagedPath -Payload $payload -Relative 'qmltooling'
foreach ($style in @('Fusion', 'Imagine', 'Material', 'Universal')) {
    Remove-StagedStyle -Payload $payload -Name $style
}

Write-Host '[3/6] Resolving all executable, DLL and Qt plugin dependencies'
# OpenSSL may be loaded dynamically by Qt's TLS backend.
foreach ($name in @('libssl-3-x64.dll', 'libcrypto-3-x64.dll')) {
    Copy-Required "$VcpkgPrefix\bin\$name" "$payload\$name"
}
$search = "$VcpkgPrefix/bin;$QtDirectory/bin;$BuildDirectory/OrbitalCore;$BuildDirectory/StrewnField"
Invoke-Checked $cmake @("-DPAYLOAD=$payload", "-DDUMPBIN=$dumpbin", "-DSEARCH_DIRECTORIES=$search",
    '-DCOPY_DEPENDENCIES=ON', "-DREPORT=$stage\dependencies-source.txt", '-P', $dependencyScript)
$ports = @('eigen3', 'nlohmann-json')
$dllNames = @(Get-ChildItem -LiteralPath $payload -Filter '*.dll' | ForEach-Object Name)
$installed = Split-Path -Parent $VcpkgPrefix
foreach ($list in Get-ChildItem -LiteralPath "$installed\vcpkg\info" -Filter '*_x64-windows.list') {
    foreach ($line in Get-Content -LiteralPath $list.FullName) {
        if ($line -match '^x64-windows/bin/([^/]+\.dll)$' -and $dllNames -contains $Matches[1]) {
            $ports += ($list.BaseName -replace '_[^_]+_x64-windows$', '')
            break
        }
    }
}
foreach ($port in $ports | Sort-Object -Unique) {
    Copy-Required "$VcpkgPrefix\share\$port\copyright" "$payload\licenses\$port.txt"
}

Write-Host '[4/6] Verifying the staged package without developer PATH'
Invoke-PayloadVerification -Payload $payload -Repository $Repository -BuildDirectory $BuildDirectory `
    -QtDirectory $QtDirectory -CMake $cmake -Dumpbin $dumpbin -ReportDirectory "$stage\verification" -DependencyScript $dependencyScript

foreach ($catalog in @(@('OBSCODE.dat', 'obscode.7z'), @('RADCODE.dat', 'radcode.7z'),
        @('linux_p1550p2650.440', 'de440.7z'))) {
    Invoke-Checked $archivegen @("$payload\$($catalog[1])", "$payload\$($catalog[0])")
    Remove-StagedPath -Payload $payload -Relative $catalog[0]
}

$manifest = [ordered]@{
    application_version = $version
    generated_utc = [DateTime]::UtcNow.ToString('o')
    repository = $Repository
    build_directory = $BuildDirectory
    vcpkg_prefix = $VcpkgPrefix
    files = @(Get-ChildItem -LiteralPath $payload -Recurse -File | Sort-Object FullName | ForEach-Object {
        [ordered]@{ path = $_.FullName.Substring($payload.Length + 1).Replace('\', '/');
            bytes = $_.Length; sha256 = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash }
    })
}
$manifest | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath "$stage\verification\payload-manifest.json" -Encoding UTF8
if ($PrepareUpdateRepository) {
    Invoke-Checked $repogen @('-p', "$stage\packages", "$stage\repository")
    Require-Path "$stage\repository\Updates.xml"
    Write-Host "Update repository: $stage\repository"
}
if ($StageOnly) {
    Write-Host "Verified staged package: $payload"
    return
}

Write-Host '[5/6] Building the hybrid installer'
$candidate = "$stage\TrajectaSetup-x86_64.exe"
Invoke-Checked $binarycreator @('--hybrid', '--compression', '5', '-c', "$stage\config\config.xml",
    '-p', "$stage\packages", $candidate)
Require-Path $candidate
if ((Get-Item -LiteralPath $candidate).Length -lt 1MB) { throw 'Installer output is unexpectedly small' }

Write-Host '[6/6] Replacing the current installer and package'
$data = Assert-Local "$root\packages\com.trajecta.app\data"
if (Test-Path -LiteralPath $data) { Remove-Item -LiteralPath $data -Recurse -Force }
$output = Assert-Local "$root\TrajectaSetup-x86_64.exe"
Move-Item -LiteralPath (Assert-Local $payload) -Destination $data
Move-Item -LiteralPath (Assert-Local $candidate) -Destination $output -Force
Write-Host "Installer: $output"
} finally {
    if (!$StageOnly -and (Test-Path -LiteralPath $stage)) {
        Remove-Item -LiteralPath (Assert-Local $stage) -Recurse -Force
    }
}
