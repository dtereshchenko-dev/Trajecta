function Controller() {
    this.removeUserData = false;
    installer.installDirectoryChanged.connect(this, Controller.prototype.onTargetDirectoryEdited);
    installer.uninstallationStarted.connect(this, Controller.prototype.onUninstallationStarted);
    installer.uninstallationFinished.connect(this, Controller.prototype.onUninstallationFinished);
}

Controller.prototype.onUninstallationStarted = function () {
    if (!installer.isUninstaller() || installer.value("PreserveUserData") === "1") return;
    if (installer.isProcessRunning("Trajecta")) {
        if (installer.isCommandLineInstance()) console.log("Close Trajecta before removing it.");
        else QMessageBox.warning("trajecta.running", "Trajecta", "Close Trajecta before removing it.");
        return;
    }
    if (installer.isCommandLineInstance()) {
        this.removeUserData = installer.value("RemoveUserData") === "1";
    } else {
        var message =
            "Also permanently delete this user's Trajecta settings, downloaded data and saved calculations?\n\n" +
            "Choose No to keep them. This cannot be undone.";
        this.removeUserData =
            QMessageBox.question(
                "trajecta.removeUserData",
                "Remove Trajecta data",
                message,
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            ) === QMessageBox.Yes;
    }
};

Controller.prototype.onUninstallationFinished = function () {
    if (
        !installer.isUninstaller() ||
        installer.value("PreserveUserData") === "1" ||
        installer.status !== QInstaller.Success ||
        !this.removeUserData
    )
        return;
    var script =
        'set -eu; data=${XDG_DATA_HOME:-"$HOME/.local/share"}; config=${XDG_CONFIG_HOME:-"$HOME/.config"}; ' +
        'case "$data" in /*) ;; *) exit 1;; esac; case "$config" in /*) ;; *) exit 1;; esac; ' +
        '[ "$data" != / ] && [ "$config" != / ] && [ ! -L "$data/Trajecta" ] && [ ! -L "$config/Trajecta" ]; ' +
        'rm -rf -- "$data/Trajecta"; rm -f -- "$config/Trajecta/Trajecta.conf"; ' +
        'rmdir -- "$config/Trajecta" 2>/dev/null || true';
    var result = installer.execute("/bin/sh", ["-c", script]);
    if (result.length < 2 || result[1] !== 0) {
        if (installer.isCommandLineInstance()) console.log("Could not remove Trajecta user data.");
        else
            QMessageBox.warning(
                "trajecta.removeUserDataFailed",
                "Trajecta",
                "Could not remove Trajecta user data."
            );
    }
};
