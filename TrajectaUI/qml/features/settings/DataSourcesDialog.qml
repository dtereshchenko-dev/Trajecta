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
    onAboutToShow: sources.refresh()
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
            CustomButton { iconName: "rotate-ccw"; flat: true; accessibleLabel: qsTr("Check again"); enabled: !dialog.sources.checking; onClicked: dialog.sources.refresh() }
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
                        Layout.fillWidth: true
                        spacing: 4
                        RowLayout {
                            Layout.fillWidth: true
                            Layout.topMargin: 12
                            CustomText { text: row.modelData.title; font.pixelSize: 12; font.weight: Font.DemiBold; Layout.fillWidth: true; elide: Text.ElideRight }
                            CustomText { text: row.modelData.status; color: row.modelData.valid && row.modelData.onDisk ? AppTheme.accent : AppTheme.amber; font.pixelSize: 10 }
                        }
                        CustomText { text: row.modelData.detail; secondary: true; font.pixelSize: 10; wrapMode: Text.Wrap; elide: Text.ElideNone; Layout.fillWidth: true }
                        CustomText { text: row.modelData.impact; secondary: true; font.pixelSize: 10; wrapMode: Text.Wrap; elide: Text.ElideNone; Layout.fillWidth: true }
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
