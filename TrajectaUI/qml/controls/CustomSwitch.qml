import QtQuick
import QtQuick.Controls.Basic
import TrajectaUI

Switch {
    id: control
    implicitWidth: 38
    implicitHeight: 33
    padding: 3
    opacity: enabled ? 1 : 0.45
    InteractionCursor {}
    indicator: Rectangle {
        x: 3
        y: (control.height - height) / 2
        width: 32
        height: 19
        radius: 10
        color: control.checked ? AppTheme.accent : AppTheme.controlLine
        border.color: control.activeFocus ? AppTheme.ink : "transparent"
        Rectangle {
            x: control.checked ? 16 : 3
            y: 3
            width: 13
            height: 13
            radius: 7
            color: AppTheme.surface
            Behavior on x { NumberAnimation { duration: 100 } }
        }
    }
}
