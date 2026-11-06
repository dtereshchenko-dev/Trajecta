import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../shared/UI.js" as UI

Rectangle {
    id: bar
    required property ApplicationSession session
    property string pageName: qsTr("Solution")
    property bool compact: false
    color: AppTheme.surface
    implicitHeight: UI.commandBarHeight
    RowLayout {
        anchors.fill: parent
        anchors.leftMargin: bar.compact ? 20 : 28
        anchors.rightMargin: bar.compact ? 20 : 28
        spacing: 12
        CustomText { text: qsTr("Workspace"); secondary: true; font.pixelSize: 12 }
        CustomIcon { name: "chevron-right"; size: 12 }
        CustomText { text: bar.pageName; font.pixelSize: 12 }
        Item { Layout.fillWidth: true }
        CustomText { text: bar.session.inputName; visible: text.length > 0; secondary: true; font.pixelSize: 10; Layout.maximumWidth: bar.width * 0.20 }
        RowLayout {
            spacing: 6
            CustomButton { iconName: "folder"; flat: true; accessibleLabel: qsTr("Open working folder"); enabled: bar.session.solution.available; onClicked: bar.session.openDirectory() }
            Rectangle { implicitWidth: 1; implicitHeight: 20; color: AppTheme.line; Layout.leftMargin: 7; Layout.rightMargin: 7 }
            CustomButton {
                objectName: "run-button"
                text: bar.session.busy ? qsTr("Cancel") : bar.session.solution.available ? qsTr("Re-solve") : qsTr("Solve")
                iconName: bar.session.busy ? "square" : "play"
                primary: !bar.session.busy
                enabled: bar.session.busy ? !bar.session.cancelling : bar.session.canRun
                onClicked: bar.session.busy ? bar.session.cancel() : bar.session.run()
            }
        }
    }
    Rectangle { anchors.bottom: parent.bottom; height: 1; width: parent.width; color: AppTheme.line }
}
