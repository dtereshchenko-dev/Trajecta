function Component()
{
    component.loaded.connect(this, Component.prototype.onLoaded);
}

Component.prototype.onLoaded = function()
{
    if (installer.isInstaller()) {
        installer.installationFinished.connect(this, Component.prototype.onInstallFinished);

        // Hide unnecessary pages
        installer.setDefaultPageVisible(QInstaller.ComponentSelection, false);
        installer.setDefaultPageVisible(QInstaller.StartMenuSelection, false);

        // Kill running app if needed
        installer.execute("taskkill", ["/IM", "Trajecta.exe", "/F"]);

        // Auto-remove previous installation if exists
        var targetDir = installer.value("TargetDir");
        var maintenanceTool = targetDir + "/Uninstall Trajecta.exe";
        if (installer.fileExists(maintenanceTool)) {
            installer.gainAdminRights();
            installer.execute(maintenanceTool, ["purge", "PreserveUserData=1", "RemoveTargetDir=false"]);
        }

    }
}

Component.prototype.onInstallFinished = function()
{
    if (installer.isInstaller() && installer.status === QInstaller.Success) {
        var targetDir = installer.value("TargetDir");
        installer.execute("cmd", ["/c",
            "for /r \"" + targetDir + "\" %f in (*.tmpUpdate) do @del /q /f \"%f\""]);
        var bootstrapDir = targetDir + "/installerResources/catalog-bootstrap";
        if (installer.fileExists(bootstrapDir))
            installer.execute("cmd", ["/c", "rmdir", "\"" + bootstrapDir + "\""]);
    }

    try {
        if (installer.isInstaller() && installer.status === QInstaller.Success) {
            var page = gui.pageWidgetByObjectName("FinishedPage");
            if (page) {
                var checkBox = page.RunItCheckBox;
                if (checkBox)
                    checkBox.checked = true;
            }
        }
    } catch(e) {}
}

Component.prototype.createOperations = function()
{
    component.createOperations();

    if (systemInfo.productType === "windows") {
        // Start menu shortcut
        component.addOperation("CreateShortcut",
            "@TargetDir@/Trajecta.exe",
            "@StartMenuDir@/Trajecta.lnk",
            "workingDirectory=@TargetDir@",
            "iconPath=@TargetDir@/icon.ico",
            "description=Trajecta");

        // Desktop shortcut
        component.addOperation("CreateShortcut",
            "@TargetDir@/Trajecta.exe",
            "@DesktopDir@/Trajecta.lnk",
            "workingDirectory=@TargetDir@",
            "iconPath=@TargetDir@/icon.ico",
            "description=Trajecta");
    }
}

Component.prototype.createOperationsForArchive = function(archive)
{
    var filename = "";
    var subdirectory = "observatories";
    if (archive.indexOf("obscode.7z") !== -1)
        filename = "OBSCODE.dat";
    else if (archive.indexOf("radcode.7z") !== -1)
        filename = "RADCODE.dat";
    else if (archive.indexOf("de440.7z") !== -1) {
        filename = "linux_p1550p2650.440";
        subdirectory = "ephemerides";
    }
    if (!filename) {
        component.createOperationsForArchive(archive);
        return;
    }

    var local = installer.environmentVariable("LOCALAPPDATA");
    if (!local)
        throw "LOCALAPPDATA is unavailable; cannot install runtime data.";
    var directory = local + "/Trajecta/Trajecta/" + subdirectory;
    if (!installer.fileExists(directory + "/" + filename)) {
        var staging = "@TargetDir@/installerResources/catalog-bootstrap";
        component.addOperation("Extract", archive, staging);
        component.addOperation("Mkdir", directory, "UNDOOPERATION", "");
        component.addOperation("Copy", staging + "/" + filename, directory + "/" + filename, "UNDOOPERATION", "");
        component.addOperation("Delete", staging + "/" + filename, "UNDOOPERATION", "");
    }
}
