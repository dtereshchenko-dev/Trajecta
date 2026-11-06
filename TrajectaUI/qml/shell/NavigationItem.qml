import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI

Button {
    id: control
    property string iconName: ""
    property string count: ""
    property bool warning: false
    property int labelSize: 13
    property bool compact: false
    implicitHeight: 39
    leftPadding: compact ? 0 : 12
    rightPadding: compact ? 0 : 12
    topPadding: 10
    bottomPadding: 10
    hoverEnabled: true
    InteractionCursor {}
    contentItem: Item {
        implicitWidth: glyph.size + (control.compact ? 0 : 11 + label.implicitWidth + (counter.visible ? 11 + counter.implicitWidth : 0) + (warningIcon.visible ? 11 + warningIcon.size : 0))
        implicitHeight: Math.max(glyph.size, label.implicitHeight)
        CustomIcon {
            id: glyph
            x: control.compact ? (parent.width - width) / 2 : 0
            anchors.verticalCenter: parent.verticalCenter
            name: control.iconName
            color: control.checked ? AppTheme.accent : AppTheme.muted
        }
        CustomText {
            id: label
            objectName: "navigation-label"
            visible: !control.compact
            anchors.left: glyph.right
            anchors.leftMargin: 11
            anchors.right: counter.visible ? counter.left : warningIcon.visible ? warningIcon.left : parent.right
            anchors.rightMargin: counter.visible || warningIcon.visible ? 11 : 0
            anchors.verticalCenter: parent.verticalCenter
            text: control.text
            font.pixelSize: control.labelSize
            color: control.checked ? AppTheme.accent : control.hovered ? AppTheme.ink : AppTheme.muted
            font.weight: control.checked ? Font.DemiBold : Font.Normal
        }
        CustomText {
            id: counter
            objectName: "navigation-count"
            text: control.count
            visible: !control.compact && text.length > 0
            anchors.right: warningIcon.visible ? warningIcon.left : parent.right
            anchors.rightMargin: warningIcon.visible ? 11 : 0
            anchors.verticalCenter: parent.verticalCenter
            secondary: true
            numeric: true
            font.pixelSize: 11
        }
        CustomIcon {
            id: warningIcon
            objectName: "navigation-warning"
            visible: control.warning
            x: parent.width - width
            y: control.compact ? 0 : (parent.height - height) / 2
            name: "triangle-alert"
            color: AppTheme.amber
            size: control.compact ? 12 : 15
        }
    }
    Accessible.name: control.warning ? text + ": " + qsTr("Attention required") : text
    background: Rectangle { radius: 4; color: control.checked ? AppTheme.accentTint : control.hovered ? AppTheme.hover : "transparent"; border.color: control.activeFocus ? AppTheme.accent : "transparent" }
}
