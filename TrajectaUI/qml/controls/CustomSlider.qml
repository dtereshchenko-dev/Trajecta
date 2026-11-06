import QtQuick
import QtQuick.Controls.Basic
import TrajectaUI

Slider {
    id: control
    implicitWidth: 200
    implicitHeight: 33
    padding: 8
    opacity: enabled ? 1 : 0.4
    InteractionCursor {}
    background: Rectangle {
        x: control.leftPadding
        y: (control.height - height) / 2
        width: control.availableWidth
        height: 3
        radius: 1
        color: AppTheme.line
        Rectangle { width: control.visualPosition * parent.width; height: parent.height; color: AppTheme.accent; radius: 1 }
    }
    handle: Rectangle {
        x: control.leftPadding + control.visualPosition * (control.availableWidth - width)
        y: (control.height - height) / 2
        width: 12
        height: 12
        radius: 6
        color: control.pressed ? AppTheme.accentStrong : AppTheme.accent
        border.width: control.activeFocus ? 2 : 0
        border.color: AppTheme.ink
    }
}
