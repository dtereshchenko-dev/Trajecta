import QtQuick
import QtQuick.Controls.Basic
import TrajectaUI
import "../shared/UI.js" as UI

Button {
    id: tab
    implicitWidth: implicitContentWidth
    implicitHeight: 46
    padding: 0
    leftPadding: 0
    rightPadding: 0
    topPadding: 0
    bottomPadding: 0
    hoverEnabled: true
    InteractionCursor {}
    font.family: UI.font
    font.pixelSize: 12
    font.weight: checked ? Font.DemiBold : Font.Normal
    contentItem: CustomText { text: tab.text; font: tab.font; color: tab.checked ? AppTheme.accent : tab.hovered ? AppTheme.ink : AppTheme.muted }
    background: Item {
        Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 2; color: tab.checked || tab.activeFocus ? AppTheme.accent : "transparent" }
    }
}
