function Component() {
    component.loaded.connect(this, Component.prototype.onLoaded);
}

Component.prototype.onLoaded = function () {
    installer.installationFinished.connect(this, Component.prototype.onInstallFinished);
    if (!installer.isCommandLineInstance()) {
        var finishedPage = gui.pageByObjectName("FinishedPage");
        if (finishedPage) finishedPage.entered.connect(this, Component.prototype.onFinishedPageEntered);
    }
    if (installer.isInstaller()) {

        // Hide unnecessary pages
        installer.setDefaultPageVisible(QInstaller.ComponentSelection, false);
        installer.setDefaultPageVisible(QInstaller.StartMenuSelection, false);
    }
};

Component.prototype.onFinishedPageEntered = function () {
    if (!installer.isUpdater() || installer.status !== QInstaller.Success) return;
    var page = gui.pageByObjectName("FinishedPage");
    if (page && page.RunItCheckBox) {
        page.RunItCheckBox.checked = true;
        gui.clickButton(buttons.FinishButton, 200);
    }
};

Component.prototype.onInstallFinished = function () {
    if (installer.status === QInstaller.Success) {
        var local = installer.environmentVariable("LOCALAPPDATA");
        var pckBootstrap = local + "/Trajecta/Trajecta/earth_orientation/.installer-bootstrap";
        if (local && installer.fileExists(pckBootstrap))
            installer.execute("cmd", ["/c", "rmdir", '"' + pckBootstrap + '"']);
    }
    if (installer.isInstaller() && installer.status === QInstaller.Success) {
        var targetDir = installer.value("TargetDir");
        installer.execute("cmd", ["/c", 'for /r "' + targetDir + '" %f in (*.tmpUpdate) do @del /q /f "%f"']);
        var bootstrapDir = targetDir + "/installerResources/catalog-bootstrap";
        if (installer.fileExists(bootstrapDir))
            installer.execute("cmd", ["/c", "rmdir", '"' + bootstrapDir + '"']);
    }

    try {
        if (installer.isInstaller() && installer.status === QInstaller.Success) {
            var page = gui.pageWidgetByObjectName("FinishedPage");
            if (page) {
                var checkBox = page.RunItCheckBox;
                if (checkBox) checkBox.checked = true;
            }
        }
    } catch (e) {}
};

Component.prototype.createOperations = function () {
    component.createOperations();

    if (systemInfo.productType === "windows") {
        // Start menu shortcut
        component.addOperation(
            "CreateShortcut",
            "@TargetDir@/Trajecta.exe",
            "@StartMenuDir@/Trajecta.lnk",
            "workingDirectory=@TargetDir@",
            "iconPath=@TargetDir@/icon.ico",
            "description=Trajecta"
        );

        // Desktop shortcut
        component.addOperation(
            "CreateShortcut",
            "@TargetDir@/Trajecta.exe",
            "@DesktopDir@/Trajecta.lnk",
            "workingDirectory=@TargetDir@",
            "iconPath=@TargetDir@/icon.ico",
            "description=Trajecta"
        );
    }
};

Component.prototype.createOperationsForArchive = function (archive) {
    var filename = "";
    var subdirectory = "observatories";
    if (archive.indexOf("obscode.7z") !== -1) filename = "OBSCODE.dat";
    else if (archive.indexOf("radcode.7z") !== -1) filename = "RADCODE.dat";
    else if (archive.indexOf("de440.7z") !== -1) {
        filename = "linux_p1550p2650.440";
        subdirectory = "ephemerides";
    }
    else if (archive.indexOf("earth_historic.7z") !== -1) {
        filename = "earth_historic.bpc";
        subdirectory = "earth_orientation";
    }
    else if (archive.indexOf("earth_latest.7z") !== -1) {
        filename = "earth_latest.bpc";
        subdirectory = "earth_orientation";
    }
    else if (archive.indexOf("earth_predict.7z") !== -1) {
        filename = "earth_predict.bpc";
        subdirectory = "earth_orientation";
    }
    else if (archive.indexOf("eopc04-data.7z") !== -1) {
        filename = "eopc04.62-now";
        subdirectory = "earth_orientation";
    }
    else if (archive.indexOf("eopc04-corrections.7z") !== -1) {
        filename = "eopc04.cnc";
        subdirectory = "earth_orientation";
    }
    else if (archive.indexOf("historical-time.7z") !== -1) {
        filename = "ET-UT.dat";
        subdirectory = "earth_orientation";
    }
    if (!filename) {
        component.createOperationsForArchive(archive);
        return;
    }

    var local = installer.environmentVariable("LOCALAPPDATA");
    if (!local) throw "LOCALAPPDATA is unavailable; cannot install runtime data.";
    var directory = local + "/Trajecta/Trajecta/" + subdirectory;
    if (!installer.fileExists(directory + "/" + filename)) {
        var staging = subdirectory === "earth_orientation"
            ? directory + "/.installer-bootstrap"
            : "@TargetDir@/installerResources/catalog-bootstrap";
        component.addOperation("Extract", archive, staging);
        component.addOperation("Mkdir", directory, "UNDOOPERATION", "");
        component.addOperation(
            "Copy",
            staging + "/" + filename,
            directory + "/" + filename,
            "UNDOOPERATION",
            ""
        );
        component.addOperation("Delete", staging + "/" + filename, "UNDOOPERATION", "");
    }
};
