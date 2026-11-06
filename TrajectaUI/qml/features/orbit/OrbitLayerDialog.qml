import QtQuick
import QtQuick.Layouts
import TrajectaUI

BaseDialog {
    id: dialog
    title: qsTr("Orbit display")
    property bool planets: true
    property bool grid: true
    property bool labels: true
    property bool perihelion: true
    contentItem: ColumnLayout {
        spacing: 12
        CustomCheckBox { objectName: "orbit-layer-planets"; text: qsTr("Planetary orbits"); checked: dialog.planets; onToggled: dialog.planets = checked }
        CustomCheckBox { objectName: "orbit-layer-grid"; text: qsTr("Ecliptic grid"); checked: dialog.grid; onToggled: dialog.grid = checked }
        CustomCheckBox { objectName: "orbit-layer-labels"; text: qsTr("Body labels"); checked: dialog.labels; onToggled: dialog.labels = checked }
        CustomCheckBox { objectName: "orbit-layer-perihelion"; text: qsTr("Perihelion"); checked: dialog.perihelion; onToggled: dialog.perihelion = checked }
    }
}
