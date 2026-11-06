import QtQuick
import QtQuick.Layouts
import TrajectaUI

Rectangle {
    id: strip
    required property ObservationTableModel observations
    readonly property var selection: observations.selection
    signal openRequested()
    implicitHeight: 40
    color: AppTheme.accentSoft
    Rectangle { width: 2; height: parent.height; color: AppTheme.accent }
    RowLayout {
        anchors.fill: parent
        anchors.leftMargin: 15
        anchors.rightMargin: 15
        spacing: 16
        CustomText { text: strip.selection.station || qsTr("No observation selected"); font.pixelSize: 11; font.weight: Font.DemiBold }
        CustomText { text: strip.selection.epoch ? strip.selection.epoch + " UTC" : ""; font.pixelSize: 11 }
        CustomText { text: strip.selection.raResidual === undefined ? "" : "Δα cos δ " + strip.selection.raResidual + "″"; numeric: true; font.pixelSize: 10 }
        CustomText { text: strip.selection.decResidual === undefined ? "" : "Δδ " + strip.selection.decResidual + "″"; numeric: true; font.pixelSize: 10 }
        Item { Layout.fillWidth: true }
        CustomTextLink { objectName: "observation-details-link"; enabled: strip.observations.selectedSourceRow >= 0; text: strip.selection.number ? qsTr("Observation %1").arg(strip.selection.number) : qsTr("Observation"); iconName: "arrow-up-right"; onClicked: strip.openRequested() }
    }
}
