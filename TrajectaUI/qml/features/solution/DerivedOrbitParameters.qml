pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

ColumnLayout {
    id: root
    objectName: "derived-orbit-parameters"
    required property OrbitModel orbit
    spacing: 12
    SectionHeading { title: qsTr("Derived quantities"); Layout.fillWidth: true }
    ParameterTable {
        Layout.fillWidth: true
        rowPrefix: "derived"
        rows: [
            {label: "q", value: root.orbit.perigee, unit: "AU", digits: 10},
            {label: "Q", value: root.orbit.eccentricity < 1 ? root.orbit.apogee : NaN, unit: "AU", digits: 10},
            {label: "1/a", value: root.orbit.reciprocalSemiMajorAxis, unit: "1/AU", digits: 10},
            {label: qsTr("Period"), value: root.orbit.eccentricity < 1 ? root.orbit.orbitalPeriod : NaN, unit: qsTr("days"), digits: 4},
            {label: "n", value: root.orbit.meanMotion, unit: "deg/day", digits: 8},
            {label: qsTr("True anomaly"), value: root.orbit.trueAnomaly, unit: "deg", digits: 8}
        ].map(row => ({label: row.label, value: UI.number(row.value, row.digits), unit: row.unit}))
    }
}
