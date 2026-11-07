function Component()
{
    component.addStopProcessForUpdateRequest("Trajecta");
}

Component.prototype.createOperations = function()
{
    component.createOperations();
    component.addOperation("Execute", "/bin/sh", "@TargetDir@/usr/bin/install-data.sh",
        "@TargetDir@/usr/share/trajecta-seed");
    component.addOperation("Mkdir", "@HomeDir@/.local/share/applications");
    component.addOperation("CreateDesktopEntry", "@HomeDir@/.local/share/applications/trajecta.desktop",
        "Type=Application\nName=Trajecta\nComment=Orbit determination and impact analysis\n" +
        "Exec=\"@TargetDir@/usr/bin/trajecta\"\nIcon=@TargetDir@/usr/share/icons/hicolor/scalable/apps/trajecta.svg\n" +
        "Categories=Science;Astronomy;\nTerminal=false\n");
};
