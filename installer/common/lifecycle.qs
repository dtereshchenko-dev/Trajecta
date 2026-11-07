var lastOfferedTarget = "";

Controller.prototype.checkExistingInstallation = function(target, onTargetPage)
{
    if (!installer.isInstaller())
        return false;
    target = target || installer.value("TargetDir");
    var name = systemInfo.productType === "windows" ? "/Uninstall Trajecta.exe" : "/TrajectaMaintenance";
    var tool = target + name;
    if (!installer.fileExists(tool))
        return false;

    var message = "Trajecta is already installed in this folder. Open the maintenance tool to update or remove it?\n\n" +
        "Choose No to select another installation folder. The existing installation and user data will not be removed.";
    var answer = QMessageBox.question("trajecta.alreadyInstalled", "Trajecta", message,
        QMessageBox.Yes | QMessageBox.No, QMessageBox.No);
    if (answer === QMessageBox.Yes) {
        if (installer.executeDetached(tool, [], target)) {
            gui.rejectWithoutPrompt();
            return true;
        }
        QMessageBox.warning("trajecta.maintenanceStartFailed", "Trajecta",
            "Could not open the maintenance tool. The existing installation was not changed.");
    }
    if (!onTargetPage)
        gui.clickButton(buttons.BackButton);
    return true;
};

Controller.prototype.TargetDirectoryPageCallback = function()
{
    if (installer.isInstaller())
        Controller.prototype.onTargetDirectoryEdited(installer.value("TargetDir"));
};

Controller.prototype.onTargetDirectoryEdited = function(target)
{
    if (!installer.isInstaller() || installer.isCommandLineInstance())
        return;
    if (target === lastOfferedTarget)
        return;
    var name = systemInfo.productType === "windows" ? "/Uninstall Trajecta.exe" : "/TrajectaMaintenance";
    if (!installer.fileExists(target + name)) {
        lastOfferedTarget = "";
        return;
    }
    lastOfferedTarget = target;
    Controller.prototype.checkExistingInstallation(target, true);
};

Controller.prototype.ComponentSelectionPageCallback = function()
{
    Controller.prototype.checkExistingInstallation("", false);
};

Controller.prototype.ReadyForInstallationPageCallback = function()
{
    if (installer.isInstaller()) {
        Controller.prototype.checkExistingInstallation("", false);
        return;
    }
    var process = systemInfo.productType === "windows" ? "Trajecta.exe" : "Trajecta";
    if (installer.isProcessRunning(process)) {
        QMessageBox.warning("trajecta.runningMaintenance", "Trajecta",
            "Close Trajecta before updating or removing it, then try again.");
        gui.rejectWithoutPrompt();
    }
};
