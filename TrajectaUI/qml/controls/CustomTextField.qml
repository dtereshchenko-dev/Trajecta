import QtQuick
import QtQuick.Controls.Basic
import TrajectaUI
import "../shared/UI.js" as UI

TextField {
    id: control
    implicitHeight: UI.controlHeight
    padding: 9
    color: AppTheme.ink
    placeholderTextColor: AppTheme.muted
    selectionColor: AppTheme.accent
    selectedTextColor: AppTheme.accentText
    font.family: UI.font
    font.pixelSize: 13
    selectByMouse: true
    InteractionCursor { cursorShape: Qt.IBeamCursor }
    background: Rectangle {
        color: AppTheme.surface
        radius: UI.radius
        border.color: control.activeFocus ? AppTheme.accent : AppTheme.line
    }
}
