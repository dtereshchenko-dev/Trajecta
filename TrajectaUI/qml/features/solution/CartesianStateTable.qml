pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

ColumnLayout {
    id: table
    objectName: "cartesian-state-table"
    required property var values
    required property var sigmas
    spacing: 12
    SectionHeading { title: qsTr("State vector"); Layout.fillWidth: true }
    ParameterTable {
        Layout.fillWidth: true
        rowPrefix: "state"
        parameterTitle: qsTr("Component")
        showSigma: true
        rows: table.values.map((value, index) => ({
            label: ["X", "Y", "Z", "Vx", "Vy", "Vz"][index],
            value: UI.number(value, 13),
            sigma: UI.sigma(table.sigmas[index]),
            unit: index < 3 ? "AU" : "AU/day"
        }))
    }
}
