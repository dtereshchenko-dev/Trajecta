pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

ColumnLayout {
    id: table
    objectName: "cartesian-state-table"
    required property var values
    spacing: 12
    SectionHeading { title: qsTr("State vector"); Layout.fillWidth: true }
    ParameterTable {
        Layout.fillWidth: true
        rowPrefix: "state"
        parameterTitle: qsTr("Component")
        rows: table.values.map((value, index) => ({
            label: ["X", "Y", "Z", "Vx", "Vy", "Vz"][index],
            value: UI.number(value, 13),
            unit: index < 3 ? "AU" : "AU/day"
        }))
    }
}
