pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

ColumnLayout {
    id: view
    required property MovViewModel mov
    readonly property var result: mov.result
    spacing: 20
    RowLayout {
        Layout.fillWidth: true
        CustomComboBox { objectName: "mov-runs"; Layout.fillWidth: true; model: view.mov.runs; textRole: "label"; currentIndex: view.mov.selectedRun; onActivated: view.mov.selectedRun = currentIndex; Accessible.name: qsTr("MOV runs") }
        CustomButton { iconName: "file-text"; accessibleLabel: qsTr("Open MOV HTML report"); enabled: !!view.result.report; onClicked: view.mov.openReport() }
        CustomButton { iconName: "folder"; accessibleLabel: qsTr("Open MOV archive"); enabled: !!view.result.archive; onClicked: view.mov.openArchive() }
    }
    CustomText { visible: !view.mov.runs.length && !view.mov.busy; text: qsTr("No MOV runs"); secondary: true; font.pixelSize: 12 }
    KeyValueList {
        visible: view.mov.runs.length > 0
        Layout.fillWidth: true
        rows: [[qsTr("Sampling"), view.result.sampling || ""], [qsTr("Fitted nodes"), String(view.result.fits ?? "")],
               [qsTr("Weighted nodes"), String(view.result.selected ?? "")], [qsTr("Crossing nodes"), String(view.result.hits ?? "")],
               [qsTr("Failed fits / propagations"), String(view.result.failures ?? "")],
               [qsTr("Weighted P"), view.result.probability === view.result.probabilityUpper ? String(view.result.probability || "") : view.result.probability + " to " + view.result.probabilityUpper],
               [qsTr("Time"), view.result.seconds === undefined ? "" : Number(view.result.seconds).toFixed(2) + " s"]]
    }
    CustomText { visible: !!view.result.cancelled; text: qsTr("Cancelled; no final probability"); color: AppTheme.amber; font.pixelSize: 11 }
    CustomText { visible: !!view.result.warnings; text: view.result.warnings || ""; color: AppTheme.amber; Layout.fillWidth: true; wrapMode: Text.Wrap; font.pixelSize: 11 }
}
