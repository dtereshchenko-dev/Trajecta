pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI

ColumnLayout {
    id: list
    required property RecentFiles recentFiles
    signal fileRequested(int index)
    spacing: 0
    CustomText { text: qsTr("Recent"); font.pixelSize: 14; font.weight: Font.DemiBold; Layout.bottomMargin: 12 }
    CustomText { visible: list.recentFiles.count === 0; text: qsTr("No recent files"); secondary: true; Layout.bottomMargin: 6 }
    Repeater {
        model: list.recentFiles.files
        AbstractButton {
            id: row
            required property string modelData
            required property int index
            objectName: "recent-observation-" + index
            Layout.fillWidth: true
            implicitHeight: 59
            hoverEnabled: true
            leftPadding: 0
            rightPadding: 8
            Accessible.name: list.recentFiles.fileName(index)
            InteractionCursor {}
            onClicked: list.fileRequested(index)
            Accessible.description: list.recentFiles.displayPath(index)
            background: Rectangle {
                color: row.hovered ? AppTheme.rowHover : "transparent"
                border.width: row.activeFocus ? 1 : 0
                border.color: AppTheme.accent
                Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line }
            }
            contentItem: RowLayout {
                spacing: 12
                CustomIcon { name: "file-text"; size: 20; color: AppTheme.muted }
                ColumnLayout {
                    Layout.fillWidth: true
                    spacing: 3
                    CustomText { text: list.recentFiles.fileName(row.index); Layout.fillWidth: true; font.pixelSize: 12 }
                    CustomText { text: list.recentFiles.displayPath(row.index); Layout.fillWidth: true; font.pixelSize: 10; secondary: true; elide: Text.ElideMiddle }
                }
                CustomIcon { name: "arrow-up-right"; size: 16; color: AppTheme.muted }
            }
        }
    }
}
