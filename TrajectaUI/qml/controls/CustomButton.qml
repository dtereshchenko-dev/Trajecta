import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI
import "../shared/UI.js" as UI

Button {
    id: control
    property string iconName: ""
    property string accessibleLabel: text
    property bool primary: false
    readonly property color foreground: !enabled ? AppTheme.muted : primary ? AppTheme.accentText : checked ? AppTheme.accent : flat && !hovered ? AppTheme.muted : AppTheme.ink
    implicitHeight: UI.controlHeight
    implicitWidth: text.length ? implicitContentWidth + 24 : 31
    leftPadding: text.length ? 12 : 7
    rightPadding: leftPadding
    topPadding: 7
    bottomPadding: 7
    spacing: 8
    font.family: UI.font
    font.pixelSize: 12
    font.weight: Font.Medium
    icon.source: iconName ? UI.getSVG(iconName) : ""
    icon.width: text.length ? 15 : 17
    icon.height: text.length ? 15 : 17
    icon.color: foreground
    palette.buttonText: foreground
    palette.windowText: foreground
    palette.highlightedText: foreground
    palette.brightText: foreground
    palette.highlight: AppTheme.accent
    hoverEnabled: true
    InteractionCursor {}
    opacity: enabled ? 1 : 0.45
    contentItem: RowLayout {
        spacing: control.text.length && control.iconName.length ? 8 : 0
        CustomIcon { visible: control.iconName.length > 0; name: control.iconName; color: control.foreground; size: control.text.length ? 15 : 17 }
        CustomText { visible: control.text.length > 0; text: control.text; font: control.font; color: control.foreground; horizontalAlignment: Text.AlignHCenter; Layout.fillWidth: true }
    }
    background: Rectangle {
        radius: UI.radius
        color: control.primary ? control.hovered ? AppTheme.accentStrong : AppTheme.accent : control.checked ? AppTheme.accentSoft : control.hovered ? AppTheme.nav : control.flat ? "transparent" : AppTheme.surface
        border.color: control.activeFocus ? AppTheme.accent : control.flat ? "transparent" : control.primary ? AppTheme.accent : control.hovered ? AppTheme.controlLineHover : AppTheme.controlLine
    }
    Accessible.name: accessibleLabel
}
