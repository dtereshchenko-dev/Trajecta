pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

SettingsPane {
    id: pane
    required property CalculationSettingsModel settings
    title: "Trajecta"
    Repeater {
        model: [
            [qsTr("Version"), pane.settings.applicationVersion],
            ["Qt", pane.settings.qtVersion],
            [qsTr("Numerical integrator"), "RA15 / Gauss-Radau"],
            [qsTr("Least squares"), "Levenberg-Marquardt"],
            [qsTr("Solution time scale"), "TDB"]
        ]
        SettingRow {
            id: row
            required property var modelData
            label: modelData[0]
            CustomText { text: row.modelData[1]; numeric: true; font.pixelSize: 12 }
        }
    }
}
