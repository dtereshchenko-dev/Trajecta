pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI

BaseDialog {
    id: dialog
    required property DataSourcesModel sources
    title: qsTr("Data sources")
    acceptText: qsTr("Close")
    width: Math.min(650, parent ? parent.width - 40 : 650)
    height: Math.min(730, parent ? parent.height - 40 : 730)
    function sourceGroup(id) {
        if (id === "planetary" || id === "asteroids") return qsTr("Ephemerides")
        if (id === "observatories" || id === "radar") return qsTr("Observatories")
        if (id.indexOf("earth_pck_") === 0 || id === "c04" || id === "delta_t") return qsTr("Earth orientation and time")
        if (id === "weights" || id.indexOf("debias_") === 0) return qsTr("Astrometric weights and corrections")
        return qsTr("Other data")
    }
    onAboutToShow: { sources.refresh(); sources.earthPckUpdate.check() }
    contentItem: ColumnLayout {
        spacing: 8
        RowLayout {
            Layout.fillWidth: true
            CustomText {
                text: dialog.sources.checking ? qsTr("Checking sources...") : qsTr("Runtime data")
                secondary: true
                font.pixelSize: 11
                Layout.fillWidth: true
            }
            CustomButton { iconName: "rotate-ccw"; flat: true; accessibleLabel: qsTr("Check again"); enabled: !dialog.sources.checking; onClicked: { dialog.sources.refresh(); dialog.sources.observatoryUpdate.check(); dialog.sources.earthPckUpdate.check() } }
        }
        CustomText {
            visible: dialog.sources.restartRequired
            text: qsTr("Data files changed. Restart the application to reload cached sources.")
            color: AppTheme.amber
            font.pixelSize: 10
            wrapMode: Text.Wrap
            elide: Text.ElideNone
            Layout.fillWidth: true
        }
        ScrollView {
            id: sourceList
            Layout.fillWidth: true
            Layout.fillHeight: true
            clip: true
            contentWidth: availableWidth
            ScrollBar.vertical: CustomScrollBar {
                parent: sourceList
                x: sourceList.width - width
                y: sourceList.topPadding
                height: sourceList.availableHeight
                policy: ScrollBar.AlwaysOn
            }
            ColumnLayout {
                width: Math.max(0, sourceList.availableWidth - 16)
                spacing: 0
                Repeater {
                    model: dialog.sources.entries
                    delegate: ColumnLayout {
                        id: row
                        required property var modelData
                        required property int index
                        readonly property bool pck: row.modelData.id.indexOf("earth_pck_") === 0
                        readonly property string groupName: dialog.sourceGroup(row.modelData.id)
                        readonly property bool firstInGroup: index === 0 || groupName !== dialog.sourceGroup((dialog.sources.entries[index - 1] || {}).id || "")
                        readonly property var updateState: (pck ? dialog.sources.earthPckUpdate.states[row.modelData.id] : dialog.sources.observatoryUpdate.states[row.modelData.id]) || ({})
                        Layout.fillWidth: true
                        spacing: 4
                        CustomText {
                            visible: row.firstInGroup
                            text: row.groupName
                            color: AppTheme.accent
                            font.pixelSize: 12
                            font.weight: Font.DemiBold
                            Layout.fillWidth: true
                            Layout.topMargin: row.index === 0 ? 8 : 20
                            Layout.bottomMargin: 2
                        }
                        RowLayout {
                            Layout.fillWidth: true
                            Layout.topMargin: 12
                            CustomText { text: row.modelData.title; font.pixelSize: 12; font.weight: Font.DemiBold; Layout.fillWidth: true; elide: Text.ElideRight }
                            CustomText { text: row.modelData.status; color: row.modelData.valid && row.modelData.onDisk ? AppTheme.accent : AppTheme.amber; font.pixelSize: 10 }
                        }
                        CustomText { text: row.modelData.detail; secondary: true; font.pixelSize: 10; wrapMode: Text.Wrap; elide: Text.ElideNone; Layout.fillWidth: true }
                        CustomText { text: row.modelData.impact; secondary: true; font.pixelSize: 10; wrapMode: Text.Wrap; elide: Text.ElideNone; Layout.fillWidth: true }
                        RowLayout {
                            visible: !!row.updateState.available || !!row.updateState.checking || !!row.updateState.downloading
                            Layout.fillWidth: true
                            Layout.topMargin: visible ? 5 : 0
                            CustomButton {
                                objectName: "update-source-" + row.modelData.id
                                visible: !!row.updateState.available && !row.updateState.downloading
                                enabled: !dialog.sources.checking
                                text: qsTr("Update")
                                iconName: "download"
                                onClicked: row.pck ? dialog.sources.earthPckUpdate.install(row.modelData.id) : dialog.sources.observatoryUpdate.install(row.modelData.id)
                            }
                            CustomText { visible: !!row.updateState.checking; text: qsTr("Checking for updates..."); secondary: true; font.pixelSize: 10 }
                            CustomText { visible: !!row.updateState.downloading; text: qsTr("Downloading..."); secondary: true; font.pixelSize: 10 }
                        }
                        CustomText {
                            visible: !!row.updateState.message
                            text: row.updateState.message || ""
                            secondary: true
                            font.pixelSize: 10
                            wrapMode: Text.Wrap
                            elide: Text.ElideNone
                            Layout.fillWidth: true
                        }
                        RowLayout {
                            visible: row.modelData.id === "asteroids" && (dialog.sources.sb441Download.running || !row.modelData.onDisk || (!row.modelData.active && !row.modelData.valid))
                            Layout.fillWidth: true
                            Layout.topMargin: visible ? 5 : 0
                            CustomButton {
                                objectName: "download-sb441"
                                visible: !dialog.sources.sb441Download.running
                                enabled: !dialog.sources.checking
                                text: qsTr("Download SB441 (616 MiB)")
                                iconName: "download"
                                onClicked: dialog.sources.sb441Download.start()
                            }
                            CustomProgressBar { visible: dialog.sources.sb441Download.running; value: dialog.sources.sb441Download.progress; Layout.fillWidth: true }
                            CustomButton { visible: dialog.sources.sb441Download.running; text: qsTr("Cancel"); onClicked: dialog.sources.sb441Download.cancel() }
                        }
                        CustomText {
                            visible: row.modelData.id === "asteroids" && dialog.sources.sb441Download.message.length > 0
                            text: dialog.sources.sb441Download.message
                            secondary: true
                            font.pixelSize: 10
                            wrapMode: Text.Wrap
                            elide: Text.ElideNone
                            Layout.fillWidth: true
                        }
                        Rectangle { Layout.fillWidth: true; Layout.topMargin: 12; implicitHeight: 1; color: AppTheme.line }
                    }
                }
            }
        }
    }
}
