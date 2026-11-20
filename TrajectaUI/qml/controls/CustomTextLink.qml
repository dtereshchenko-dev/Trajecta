import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI

Button {
    id: link
    property string iconName: "chevron-right"
    property bool alwaysUnderline: false
    padding: 0
    implicitHeight: 24
    implicitWidth: label.implicitWidth
    hoverEnabled: true
    InteractionCursor {}
    contentItem: RowLayout {
        id: label
        spacing: 6
        CustomText { text: link.text; font.pixelSize: 11; color: AppTheme.accent; font.underline: link.alwaysUnderline || link.hovered || link.activeFocus }
        CustomIcon { name: link.iconName; color: AppTheme.accent; size: 13 }
    }
    background: Item {}
}
