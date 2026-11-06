pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

ColumnLayout {
    id: root
    objectName: "estimated-size"
    required property SolutionViewModel solution
    readonly property var typical: solution.typicalSize
    visible: solution.sizeEstimates.length > 0
    spacing: 12
    SectionHeading {
        title: qsTr("Estimated size")
        detail: "H " + UI.number(root.solution.absoluteMagnitude, 3) + " mag \u00B7 "
            + qsTr("%1 magnitudes").arg(root.solution.magnitudeCount)
        Layout.fillWidth: true
    }
    ParameterTable {
        Layout.fillWidth: true
        rowPrefix: "size"
        parameterTitle: qsTr("Geometric albedo")
        parameterShare: 0.72
        unitShare: 0.06
        rows: root.solution.sizeEstimates
    }
    InlineValue {
        objectName: "typical-size"
        label: qsTr("Typical, albedo %1").arg(root.typical.label || "")
        value: root.typical.value || "\u2014"
        unit: root.typical.unit || ""
    }
}
