pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

ColumnLayout {
    id: root
    objectName: "derived-orbit-parameters"
    required property OrbitModel orbit
    property bool hasSigmas: false
    spacing: 12
    SectionHeading { title: qsTr("Derived quantities"); Layout.fillWidth: true }
    ParameterTable {
        Layout.fillWidth: true
        rowPrefix: "derived"
        showSigma: true
        rows: [
            {label: "q", value: root.orbit.perigee, sigma: root.orbit.sigmaQ, unit: "AU", digits: 10},
            {label: "Q", value: root.orbit.eccentricity < 1 ? root.orbit.apogee : NaN, sigma: root.orbit.sigmaApogee, unit: "AU", digits: 10},
            {label: "1/a", value: root.orbit.reciprocalSemiMajorAxis, sigma: root.orbit.sigmaReciprocalSemiMajorAxis, unit: "1/AU", digits: 10},
            {label: qsTr("Period"), value: root.orbit.eccentricity < 1 ? root.orbit.orbitalPeriod : NaN, sigma: root.orbit.sigmaPeriod, unit: qsTr("days"), digits: 4},
            {label: "n", value: root.orbit.meanMotion, sigma: root.orbit.sigmaMeanMotion, unit: "deg/day", digits: 8},
            {label: qsTr("True anomaly"), value: root.orbit.trueAnomaly, sigma: root.orbit.sigmaNu, unit: "deg", digits: 8}
        ].map(row => ({label: row.label, value: UI.number(row.value, row.digits), sigma: root.hasSigmas && Number.isFinite(row.value) ? UI.sigma(row.sigma) : "\u2014", unit: row.unit}))
    }
}
