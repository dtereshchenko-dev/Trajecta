import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

ColumnLayout {
    id: summary
    required property SolutionViewModel solution
    spacing: 12
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line }
    RowLayout {
        Layout.fillWidth: true
        Layout.topMargin: 3
        Layout.minimumHeight: 30
        CustomText { text: qsTr("Earth MOID"); secondary: true; font.pixelSize: 11; Layout.fillWidth: true }
        InlineValue { objectName: "earth-moid"; value: UI.groupedNumber(summary.solution.earthMoidKm, 1); valueSize: 13; unit: "km" }
    }
    CustomTextLink { objectName: "all-moids-link"; text: qsTr("All planets"); implicitHeight: 16; onClicked: details.open() }
    MoidDetailsDialog { id: details; solution: summary.solution }
}
