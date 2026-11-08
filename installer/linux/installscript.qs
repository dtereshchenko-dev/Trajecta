function Component() {
    component.addStopProcessForUpdateRequest("Trajecta");
    component.loaded.connect(this, Component.prototype.onLoaded);
}

Component.prototype.onLoaded = function () {
    if (!installer.isCommandLineInstance()) {
        var finishedPage = gui.pageByObjectName("FinishedPage");
        if (finishedPage) finishedPage.entered.connect(this, Component.prototype.onFinishedPageEntered);
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

Component.prototype.createOperations = function () {
    component.createOperations();
    component.addOperation(
        "Execute",
        "/bin/sh",
        "@TargetDir@/usr/bin/install-data.sh",
        "@TargetDir@/usr/share/trajecta-seed"
    );
    component.addOperation("Mkdir", "@HomeDir@/.local/share/applications");
    component.addOperation(
        "CreateDesktopEntry",
        "@HomeDir@/.local/share/applications/trajecta.desktop",
        "Type=Application\nName=Trajecta\nComment=Orbit determination and impact analysis\n" +
            'Exec="@TargetDir@/usr/bin/trajecta"\nIcon=@TargetDir@/usr/share/icons/hicolor/scalable/apps/trajecta.svg\n' +
            "Categories=Science;Astronomy;\nTerminal=false\n"
    );
};
