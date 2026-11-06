pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

SettingsPane {
    id: pane
    required property ApplicationSession session
    title: qsTr("Appearance")
    GridLayout {
        Layout.fillWidth: true
        columns: width >= 558 ? 3 : 2
        columnSpacing: 12
        rowSpacing: 12
        Repeater {
            model: [
                { key: "graphite", name: "Graphite", note: qsTr("Neutral dark with a green accent.") },
                { key: "observatory", name: "Observatory", note: qsTr("Dark blue with a turquoise accent.") },
                { key: "slate", name: "Slate", note: qsTr("Soft dark grey, low contrast.") },
                { key: "nocturne", name: "Nocturne", note: qsTr("Near-black for dark rooms and OLED.") },
                { key: "redlight", name: "Red Light", note: qsTr("Warm red with reduced blue.") },
                { key: "dim", name: "Dim", note: qsTr("Muted light without pure white.") },
                { key: "light", name: "Light", note: qsTr("Original light palette.") }
            ]
            ThemeTile {
                required property var modelData
                Layout.fillWidth: true
                Layout.preferredWidth: 178
                themeKey: modelData.key
                text: modelData.name
                description: modelData.note
                checked: pane.session.theme === themeKey
                onClicked: pane.session.theme = themeKey
            }
        }
    }
}
