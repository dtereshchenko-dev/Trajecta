pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

ColumnLayout {
    id: root
    objectName: "estimated-size"
    required property SolutionViewModel solution
    readonly property var typical: solution.typicalSize
    visible: solution.sizeEstimates.length > 0
    spacing: 12
    SectionHeading {
        title: qsTr("Estimated size")
        Layout.fillWidth: true
    }
    ParameterTable {
        Layout.fillWidth: true
        rowPrefix: "size"
        parameterTitle: qsTr("Geometric albedo")
        showSigma: true
        sigmaTitle: ""
        rows: root.solution.sizeEstimates.map(row => ({label: row.label, value: row.value, sigma: "", unit: row.unit}))
    }
    InlineValue {
        objectName: "typical-size"
        label: qsTr("Typical, albedo %1").arg(root.typical.label || "")
        value: root.typical.value || "\u2014"
        unit: root.typical.unit || ""
    }
}
