pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

ColumnLayout {
    id: panel
    required property MovViewModel mov
    required property MovCalculationForm form
    spacing: 20
    CustomText { text: qsTr("Source data"); font.pixelSize: 11; font.weight: Font.DemiBold }
    KeyValueList {
        Layout.fillWidth: true
        rows: [[qsTr("Source"), panel.mov.sourceName || qsTr("Not available")],
               [qsTr("Observations"), String(panel.mov.observationCount)],
               [qsTr("Nominal / covariance"), panel.mov.externalSource ? qsTr("Not supplied") : qsTr("Automatic")],
               [qsTr("Dynamics"), "RA15"], [qsTr("Range domain"), qsTr("Automatic admissible region")]]
    }
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line }
    CustomText { text: qsTr("Sampling"); font.pixelSize: 11; font.weight: Font.DemiBold }
    KeyValueList {
        Layout.fillWidth: true
        rows: [[qsTr("Grid nodes"), panel.form.gridValid ? String(panel.form.gridNodes) : "\u2014"],
               [qsTr("Preliminary grid"), "40 x 40"], [qsTr("Probability measure"), qsTr("Weighted MOV integral")],
               [qsTr("Maximum horizon"), qsTr("30 days")]]
    }
}
