[CmdletBinding()]
param([string]$IfwDirectory = 'C:\Qt\Tools\QtInstallerFramework\4.11')

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$base = [IO.Path]::GetFullPath([IO.Path]::GetTempPath()).TrimEnd('\') + '\'
$root = Join-Path $base ('trajecta-ifw-test-' + [guid]::NewGuid().ToString('N'))
$savedLocalAppData = $env:LOCALAPPDATA
$repoRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))

function Run-Checked([string]$Executable, [string[]]$Arguments) {
    & $Executable @Arguments
    if ($LASTEXITCODE -ne 0) { throw "$Executable failed: $LASTEXITCODE" }
}
function Write-Text([string]$Path, [string]$Content) {
    $null = New-Item -ItemType Directory -Path (Split-Path -Parent $Path) -Force
    [IO.File]::WriteAllText($Path, $Content, [Text.UTF8Encoding]::new($false))
}
function Assert-Text([string]$Path, [string]$Expected) {
    if (!(Test-Path -LiteralPath $Path) -or [IO.File]::ReadAllText($Path) -ne $Expected) {
        throw "Unexpected content: $Path"
    }
}

try {
    $null = New-Item -ItemType Directory -Path $root -Force
    $env:LOCALAPPDATA = Join-Path $root 'local'
    $packages = Join-Path $root 'packages'
    $package = Join-Path $packages 'com.trajecta.app'
    $data = Join-Path $package 'data'
    $meta = Join-Path $package 'meta'
    $repository = Join-Path $root 'repository'
    $target = Join-Path $root 'installed with spaces'
    $config = Join-Path $root 'config.xml'
    $source = Get-Content -LiteralPath (Join-Path $repoRoot 'packages/com.trajecta.app/meta/installscript.qs') -Raw
    $start = [regex]::Match($source, 'Component\.prototype\.createOperationsForArchive\s*=\s*function\s*\(\s*archive\s*\)')
    $operations = [regex]::Match($source, 'Component\.prototype\.createOperations\s*=\s*function\s*\(\s*\)')
    if (!$start.Success -or !$operations.Success) { throw 'Component script handlers not found' }
    Write-Text (Join-Path $meta 'installscript.qs') ($source.Substring(0, $operations.Index) + $source.Substring($start.Index))
    $controller = (Get-Content -LiteralPath (Join-Path $repoRoot 'config/controller.qs') -Raw) + "`n" +
        (Get-Content -LiteralPath (Join-Path $repoRoot 'common/lifecycle.qs') -Raw)
    Write-Text (Join-Path $root 'controller.qs') $controller
    $repoUrl = ([uri]$repository).AbsoluteUri
    Write-Text $config "<Installer><Name>Trajecta update smoke</Name><Version>1.0.1</Version><Title>Smoke</Title><Publisher>Test</Publisher><MaintenanceToolName>Uninstall Trajecta</MaintenanceToolName><ControlScript>controller.qs</ControlScript><RemoteRepositories><Repository><Url>$repoUrl</Url><Enabled>1</Enabled></Repository></RemoteRepositories></Installer>"
    $binarycreator = Join-Path $IfwDirectory 'bin/binarycreator.exe'
    $repogen = Join-Path $IfwDirectory 'bin/repogen.exe'
    $archivegen = Join-Path $IfwDirectory 'bin/archivegen.exe'
    foreach ($version in @('1.0.1', '1.0.2')) {
        $essential = if ($version -eq '1.0.1') { 'true' } else { 'false' }
        Write-Text (Join-Path $meta 'package.xml') "<Package><DisplayName>Trajecta</DisplayName><Description>Smoke</Description><Version>$version</Version><ReleaseDate>2026-01-01</ReleaseDate><Default>true</Default><Essential>$essential</Essential><Script>installscript.qs</Script></Package>"
        Write-Text (Join-Path $data 'app.txt') $version
        foreach ($item in @(@('OBSCODE.dat', 'obscode.7z'), @('RADCODE.dat', 'radcode.7z'), @('linux_p1550p2650.440', 'de440.7z'))) {
            $inputFile = Join-Path $root ('input/' + $item[0])
            Write-Text $inputFile $version
            Run-Checked $archivegen @((Join-Path $data $item[1]), $inputFile)
        }
        if ($version -eq '1.0.1') {
            Run-Checked $binarycreator @('--hybrid', '-c', $config, '-p', $packages, (Join-Path $root 'setup.exe'))
            Run-Checked (Join-Path $root 'setup.exe') @('--root', $target, '--default-answer', '--confirm-command', 'install')
            Assert-Text (Join-Path $target 'app.txt') '1.0.1'
            $userCatalog = Join-Path $env:LOCALAPPDATA 'Trajecta/Trajecta/observatories/OBSCODE.dat'
            Assert-Text $userCatalog '1.0.1'
            Write-Text $userCatalog 'user-updated'
            Write-Text (Join-Path $env:LOCALAPPDATA 'Trajecta/Trajecta/work/keep.txt') 'user-work'
            $metadata = Join-Path $target 'components.xml'
            $metadataHash = (Get-FileHash -LiteralPath $metadata -Algorithm SHA256).Hash
            & (Join-Path $root 'setup.exe') '--root' $target '--auto-answer' 'OverwriteTargetDirectory=Yes' '--default-answer' '--confirm-command' 'install'
            if ($LASTEXITCODE -eq 0) { throw 'Repeated installation unexpectedly succeeded.' }
            if ((Get-FileHash -LiteralPath $metadata -Algorithm SHA256).Hash -ne $metadataHash) { throw 'Repeated installation changed metadata.' }
            Assert-Text (Join-Path $target 'app.txt') '1.0.1'
            Assert-Text $userCatalog 'user-updated'
            if (Test-Path -LiteralPath (Join-Path $target 'Uninstall Trajecta.exe.new')) { throw 'Repeated installation left a new maintenance tool.' }
        } else {
            Run-Checked $repogen @('-p', $packages, $repository)
            $check = [System.Diagnostics.Process]::new()
            $check.StartInfo.FileName = Join-Path $target 'Uninstall Trajecta.exe'
            $check.StartInfo.Arguments = 'check-updates'
            $check.StartInfo.UseShellExecute = $false
            $check.StartInfo.RedirectStandardOutput = $true
            $check.StartInfo.RedirectStandardError = $true
            if (!$check.Start()) { throw 'Could not start update check' }
            $output = $check.StandardOutput.ReadToEnd()
            $errors = $check.StandardError.ReadToEnd()
            $check.WaitForExit()
            if ($check.ExitCode -ne 0 -or $output -notmatch '(?s)<updates>.*id="com\.trajecta\.app".*/>.*</updates>') {
                throw "Update check failed: $output $errors"
            }
            & (Join-Path $target 'Uninstall Trajecta.exe') '--default-answer' '--confirm-command' 'update'
            if ($LASTEXITCODE -ne 0) { throw "Maintenance tool failed or requested a restart: $LASTEXITCODE" }
            Assert-Text (Join-Path $target 'app.txt') '1.0.2'
            Assert-Text $userCatalog 'user-updated'
            Assert-Text (Join-Path $env:LOCALAPPDATA 'Trajecta/Trajecta/work/keep.txt') 'user-work'
            Assert-Text (Join-Path $env:LOCALAPPDATA 'Trajecta/Trajecta/ephemerides/linux_p1550p2650.440') '1.0.1'
        }
    }
    Write-Host 'Hybrid update and user data preservation passed.'
} finally {
    $env:LOCALAPPDATA = $savedLocalAppData
    if (([IO.Path]::GetFullPath($root)).StartsWith($base + 'trajecta-ifw-test-', [StringComparison]::OrdinalIgnoreCase) -and (Test-Path -LiteralPath $root)) {
        Remove-Item -LiteralPath $root -Recurse -Force
    }
}

$global:LASTEXITCODE = 0
