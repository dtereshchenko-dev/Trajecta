import QtQuick
import QtQuick.Layouts
import TrajectaUI

Rectangle {
    id: hint
    property var observation: ({})
    property point anchorPoint
    visible: observation.sourceRow !== undefined
    width: Math.min(280, parent.width - 16)
    height: content.implicitHeight + 20
    x: Math.max(8, Math.min(anchorPoint.x + 12, parent.width - width - 8))
    y: Math.max(8, Math.min(anchorPoint.y - height - 12, parent.height - height - 8))
    color: AppTheme.surface
    border.color: AppTheme.controlLine
    radius: 4
    z: 100
    ColumnLayout {
        id: content
        anchors { left: parent.left; right: parent.right; top: parent.top; margins: 10 }
        spacing: 6
        CustomText {
            Layout.fillWidth: true
            text: qsTr("Observation %1 · %2").arg(hint.observation.number ?? "").arg(hint.observation.station ?? "")
            font.pixelSize: 12
            font.weight: Font.DemiBold
        }
        CustomText { Layout.fillWidth: true; text: (hint.observation.epoch ?? "") + " UTC"; numeric: true; secondary: true; font.pixelSize: 11 }
        RowLayout {
            Layout.fillWidth: true
            CustomText { text: "Δα cos δ"; secondary: true; font.pixelSize: 11; Layout.fillWidth: true }
            InlineValue { value: hint.observation.raResidual ?? ""; unit: "arcsec" }
        }
        RowLayout {
            Layout.fillWidth: true
            CustomText { text: "Δδ"; secondary: true; font.pixelSize: 11; Layout.fillWidth: true }
            InlineValue { value: hint.observation.decResidual ?? ""; unit: "arcsec" }
        }
        CustomText {
            text: hint.observation.included ? qsTr("Included in fit") : qsTr("Excluded from fit")
            color: hint.observation.included ? AppTheme.accent : AppTheme.red
            font.pixelSize: 11
        }
    }
}
