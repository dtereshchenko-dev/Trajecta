pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

Button {
    id: tile
    required property string themeKey
    property string description: ""
    readonly property var colors: UI.themes[themeKey]
    objectName: "theme-" + themeKey
    implicitHeight: 152
    implicitWidth: 200
    padding: 1
    checkable: true
    Accessible.name: text
    InteractionCursor {}
    background: Rectangle { color: AppTheme.surface; radius: 6; border.width: tile.checked || tile.activeFocus ? 2 : 1; border.color: tile.checked || tile.activeFocus ? AppTheme.accent : tile.hovered ? AppTheme.controlLineHover : AppTheme.controlLine }
    contentItem: Column {
        spacing: 0
        Item {
            width: parent.width
            height: 76
            Rectangle { anchors.fill: parent; anchors.margins: 1; radius: 4; color: tile.colors.surface }
            Rectangle {
                x: 1; y: 1; width: 33; height: 75; color: tile.colors.nav
                Column {
                    x: 6; y: 9; spacing: 5
                    Repeater {
                        model: 4
                        Rectangle { required property int index; width: 21; height: 5; radius: 2; color: index === 0 ? tile.colors.accent : tile.colors.line }
                    }
                }
            }
            Rectangle { x: 34; height: 76; width: 1; color: tile.colors.line }
            Column {
                x: 46; y: 10; width: parent.width - 58; spacing: 6
                Rectangle { width: parent.width * 0.54; height: 7; radius: 2; color: tile.colors.ink; opacity: 0.8 }
                Rectangle { width: parent.width * 0.8; height: 5; radius: 2; color: tile.colors.muted }
            }
            Row {
                x: 46; y: 55; spacing: 5
                Repeater {
                    model: [tile.colors.accent, tile.colors.blue, tile.colors.red, tile.colors.amber]
                    Rectangle { required property string modelData; width: 9; height: 9; radius: 5; color: modelData }
                }
            }
            Rectangle { x: parent.width - 43; y: 54; width: 32; height: 11; radius: 2; color: tile.colors.nav }
        }
        Rectangle { width: parent.width; height: 1; color: AppTheme.line }
        RowLayout {
            x: 12; width: parent.width - 24; height: 73; spacing: 8
            ColumnLayout {
                spacing: 3
                Layout.fillWidth: true
                CustomText { text: tile.text; font.pixelSize: 12; font.weight: Font.DemiBold; Layout.fillWidth: true }
                CustomText { text: tile.description; font.pixelSize: 10; secondary: true; wrapMode: Text.Wrap; elide: Text.ElideNone; Layout.fillWidth: true }
            }
            CustomIcon { name: "check"; size: 16; color: AppTheme.accent; opacity: tile.checked ? 1 : 0 }
        }
    }
}
