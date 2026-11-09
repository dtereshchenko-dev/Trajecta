function Controller() {
    this.removeUserData = false;
    this.targetDir = "";
    installer.installDirectoryChanged.connect(this, Controller.prototype.onTargetDirectoryEdited);
    installer.uninstallationStarted.connect(this, Controller.prototype.onUninstallationStarted);
    installer.uninstallationFinished.connect(this, Controller.prototype.onUninstallationFinished);
}

Controller.prototype.onUninstallationStarted = function () {
    if (
        systemInfo.productType !== "windows" ||
        !installer.isUninstaller() ||
        installer.value("PreserveUserData") === "1"
    )
        return;

    if (installer.isProcessRunning("Trajecta.exe")) {
        if (installer.isCommandLineInstance()) console.log("Close Trajecta before removing it.");
        else QMessageBox.warning("trajecta.running", "Trajecta", "Close Trajecta before removing it.");
        return;
    }

    this.targetDir = installer.value("TargetDir");
    if (installer.isCommandLineInstance()) {
        this.removeUserData = installer.value("RemoveUserData") === "1";
        return;
    }

    var localData = installer.environmentVariable("LOCALAPPDATA");
    var message =
        "Also permanently delete Trajecta settings and user data?\n\n" +
        "This includes saved calculations in the work folder, downloaded data, " +
        "the examples folder (including edited files), and settings in the Windows registry.\n\n" +
        "Work folder: " +
        this.targetDir +
        "/work\n" +
        "Local data: " +
        (localData ? localData + "/Trajecta" : "%LOCALAPPDATA%/Trajecta") +
        "\nRegistry: HKEY_CURRENT_USER\\Software\\Trajecta\n\n" +
        "Choose No to keep your data. This cannot be undone.";
    this.removeUserData =
        QMessageBox.question(
            "trajecta.removeUserData",
            "Remove Trajecta data",
            message,
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        ) === QMessageBox.Yes;
};

Controller.prototype.onUninstallationFinished = function () {
    if (
        systemInfo.productType !== "windows" ||
        !installer.isUninstaller() ||
        installer.value("PreserveUserData") === "1" ||
        installer.status !== QInstaller.Success
    )
        return;

    var failures = [];
    var runCleanup = function (label, body) {
        var script =
            "$ErrorActionPreference = 'Stop'; try { " +
            body +
            " } catch { [Console]::Out.WriteLine($_.ToString()); exit 1 }";
        var result = installer.execute("powershell.exe", [
            "-NoProfile",
            "-NonInteractive",
            "-Command",
            script
        ]);
        if (result.length < 2 || result[1] !== 0)
            failures.push(
                label + ": " + (result.length > 0 && result[0] ? result[0] : "Could not start cleanup.")
            );
    };

    if (this.removeUserData) {
        installer.dropAdminRights();
        // The current user's settings and data are independent of installation cleanup.
        runCleanup(
            "Registry",
            "$path = 'HKCU:\\Software\\Trajecta'; " +
                "if (Test-Path -LiteralPath $path) { Remove-Item -LiteralPath $path -Recurse -Force }"
        );

        var localData = installer.environmentVariable("LOCALAPPDATA");
        if (!localData) {
            failures.push("Local data: LOCALAPPDATA is unavailable.");
        } else {
            var requested = (localData + "/Trajecta").replace(/'/g, "''");
            runCleanup(
                "Local data",
                "$requested = '" +
                    requested +
                    "'; " +
                    "if (-not $env:LOCALAPPDATA -or -not [IO.Path]::IsPathRooted($env:LOCALAPPDATA)) { throw 'Invalid LOCALAPPDATA' }; " +
                    "$path = Join-Path $env:LOCALAPPDATA 'Trajecta'; " +
                    "if (-not [String]::Equals([IO.Path]::GetFullPath($requested), [IO.Path]::GetFullPath($path), [StringComparison]::OrdinalIgnoreCase)) { throw 'User data path mismatch' }; " +
                    "if (Test-Path -LiteralPath $path) { " +
                    "$item = Get-Item -LiteralPath $path -Force; " +
                    "if (($item.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) { throw 'Refusing linked data directory' }; " +
                    "Remove-Item -LiteralPath $path -Recurse -Force }"
            );
        }
    }

    var target = this.targetDir.replace(/'/g, "''");
    var validateTarget =
        "$target = '" +
        target +
        "'; " +
        "if (-not $target -or -not [IO.Path]::IsPathRooted($target)) { throw 'Invalid installation path' }; " +
        "$target = [IO.Path]::GetFullPath($target).TrimEnd([IO.Path]::DirectorySeparatorChar, [IO.Path]::AltDirectorySeparatorChar); " +
        "if ($target.Length -le [IO.Path]::GetPathRoot($target).Length) { throw 'Invalid installation path' }; " +
        "if (Test-Path -LiteralPath $target) { " +
        "$item = Get-Item -LiteralPath $target -Force; " +
        "if (($item.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) { throw 'Refusing linked installation directory' } }; ";
    var removeDirectory =
        "if (Test-Path -LiteralPath $path) { " +
        "$item = Get-Item -LiteralPath $path -Force; " +
        "if (($item.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) { throw ('Refusing linked directory: ' + $path) }; " +
        "if (Get-ChildItem -LiteralPath $path -Recurse -Force | Where-Object { ($_.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0 } | Select-Object -First 1) { throw ('Refusing linked entries: ' + $path) }; " +
        "Remove-Item -LiteralPath $path -Recurse -Force }";

    if (!installer.gainAdminRights()) {
        failures.push("Installation folder: administrative access was not granted.");
    } else {
        if (this.removeUserData) {
            runCleanup(
                "Examples",
                validateTarget + "$path = Join-Path $target 'examples'; " + removeDirectory
            );
            runCleanup(
                "Work folder",
                validateTarget + "$path = Join-Path $target 'work'; " + removeDirectory
            );
        }
        runCleanup(
            "Legacy files",
            validateTarget +
                "$resources = Join-Path $target 'installerResources'; " +
                "$path = Join-Path $resources 'com.orbitalmechanics.app'; " +
                "if (Test-Path -LiteralPath $path) { " +
                "$item = Get-Item -LiteralPath $resources -Force; " +
                "if (($item.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) { throw 'Refusing linked installer resources' }; " +
                "foreach ($name in @('OrbitalMechanics.exe', 'OrbitalMechanicsUI.exe')) { " +
                "$file = Join-Path $target $name; if (Test-Path -LiteralPath $file -PathType Leaf) { Remove-Item -LiteralPath $file -Force } }; " +
                removeDirectory +
                " }"
        );
        runCleanup(
            "Old maintenance file",
            validateTarget +
                "$path = Join-Path $target 'Uninstall Trajecta.exe.new'; " +
                "if (Test-Path -LiteralPath $path -PathType Leaf) { Remove-Item -LiteralPath $path -Force }"
        );
    }

    if (failures.length > 0) {
        var message = "Some items could not be removed:\n\n" + failures.join("\n");
        if (installer.isCommandLineInstance()) console.log(message);
        else QMessageBox.warning("trajecta.removeUserDataFailed", "Trajecta", message);
    }
};
