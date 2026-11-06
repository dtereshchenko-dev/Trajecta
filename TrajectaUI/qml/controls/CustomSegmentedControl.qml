pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI

Rectangle {
    id: control
    property var options: []
    property var disabledIndices: []
    property int currentIndex: 0
    signal activated(int index)
    implicitWidth: buttons.implicitWidth + 6
    implicitHeight: 31
    color: AppTheme.well
    border.color: AppTheme.line
    radius: 4
    RowLayout {
        id: buttons
        anchors.centerIn: parent
        spacing: 1
        Repeater {
            model: control.options
            Button {
                id: option
                required property string modelData
                required property int index
                objectName: control.objectName + "-option-" + index
                enabled: control.disabledIndices.indexOf(index) < 0
                opacity: enabled ? 1 : 0.45
                implicitHeight: 25
                implicitWidth: implicitContentWidth + 20
                padding: 0
                hoverEnabled: true
                InteractionCursor {}
                contentItem: CustomText { text: option.modelData; font.pixelSize: 11; horizontalAlignment: Text.AlignHCenter; color: control.currentIndex === option.index ? AppTheme.ink : AppTheme.muted }
                background: Rectangle { radius: 3; color: control.currentIndex === option.index ? AppTheme.surface : option.hovered ? AppTheme.hover : "transparent"; border.color: option.activeFocus ? AppTheme.accent : "transparent" }
                onClicked: control.activated(index)
            }
        }
    }
}
