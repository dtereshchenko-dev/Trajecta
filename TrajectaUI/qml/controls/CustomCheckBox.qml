import QtQuick
import QtQuick.Controls.Basic
import TrajectaUI
import "../shared/UI.js" as UI

CheckBox {
    id: control
    implicitHeight: UI.controlHeight
    spacing: 8
    padding: 2
    InteractionCursor {}
    indicator: Rectangle {
        x: control.leftPadding
        y: (control.height - height) / 2
        width: 16
        height: 16
        radius: 3
        color: control.checked ? AppTheme.accent : AppTheme.surface
        border.color: control.activeFocus || control.checked ? AppTheme.accent : AppTheme.muted
        CustomText { anchors.centerIn: parent; text: control.checked ? "✓" : ""; color: AppTheme.accentText }
    }
    contentItem: CustomText {
        text: control.text
        leftPadding: 24
    }
}
