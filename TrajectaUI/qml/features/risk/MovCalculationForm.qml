pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import QtQuick.Dialogs
import TrajectaUI

ColumnLayout {
    id: form
    objectName: "mov-form"
    required property MovViewModel mov
    property int sampling: 0
    property real rangeNodes: 80
    property real rateNodes: 80
    property real workers: 15
    property int observationSelection: 0
    property int target: 0
    property int boundary: 1
    property string startUtc: ""
    property string endUtc: ""
    readonly property bool gridValid: Number.isInteger(rangeNodes) && rangeNodes >= 8 && rangeNodes <= 250
        && Number.isInteger(rateNodes) && rateNodes >= 8 && rateNodes <= 250
        && Number.isInteger(workers) && workers >= 1 && workers <= 64
    readonly property int gridNodes: gridValid ? rangeNodes * rateNodes : 0
    readonly property var options: ({startUtc, endUtc, sampling, rangeNodes, rateNodes, workers, observationSelection, target, boundary})
    function initialize() {
        lastSeed = mov.seedUtc
        if (!mov.seedUtc.length) return
        startUtc = mov.seedUtc
        const date = new Date(startUtc.replace(" ", "T") + "Z")
        if (Number.isFinite(date.getTime())) endUtc = new Date(date.getTime() + 2 * 86400000).toISOString().replace("T", " ").replace("Z", "")
    }
    Component.onCompleted: initialize()
    Connections { target: form.mov; function onChanged() { if (form.mov.seedUtc !== form.lastSeed) { form.lastSeed = form.mov.seedUtc; form.initialize() } } }
    property string lastSeed: ""
    spacing: 20
    CustomText { text: "MOV"; font.pixelSize: 18; font.weight: Font.DemiBold }
    FormSection { title: qsTr("Observations") }
    RowLayout {
        Layout.fillWidth: true
        CustomTextField { Layout.fillWidth: true; text: form.mov.sourceName; readOnly: true; Accessible.name: qsTr("MOV observation source") }
        CustomButton { iconName: "folder-open"; accessibleLabel: qsTr("Open MOV observations"); onClicked: fileDialog.open() }
        CustomButton { visible: form.mov.externalSource; iconName: "rotate-ccw"; accessibleLabel: qsTr("Use current solution observations"); onClicked: form.mov.useAcceptedSolution() }
    }
    FormField {
        label: qsTr("Observation set")
        CustomComboBox {
            objectName: "mov-observation-set"
            Layout.fillWidth: true
            model: form.mov.externalSource ? [qsTr("All observations")] : [qsTr("Accepted fit mask"), qsTr("All observations")]
            currentIndex: form.mov.externalSource ? 0 : form.observationSelection
            Accessible.name: qsTr("MOV observation set")
            onActivated: form.observationSelection = currentIndex
        }
    }
    FormSection { title: qsTr("Calculation interval") }
    GridLayout {
        Layout.fillWidth: true
        columns: width < 460 ? 1 : 2
        columnSpacing: 20
        rowSpacing: 18
        FormField { label: qsTr("Start, UTC"); CustomTextField { objectName: "mov-start"; Layout.fillWidth: true; text: form.startUtc; placeholderText: "YYYY-MM-DD HH:MM:SS"; onTextEdited: form.startUtc = text } }
        FormField { label: qsTr("End, UTC"); CustomTextField { objectName: "mov-end"; Layout.fillWidth: true; text: form.endUtc; placeholderText: "YYYY-MM-DD HH:MM:SS"; onTextEdited: form.endUtc = text } }
        FormField {
            label: qsTr("Target")
            CustomComboBox { objectName: "mov-target"; Layout.fillWidth: true; model: [qsTr("Earth"), qsTr("Moon"), qsTr("Mars")]; currentIndex: form.target; onActivated: form.target = currentIndex; Accessible.name: qsTr("MOV target body") }
        }
        FormField {
            label: qsTr("Boundary")
            CustomComboBox { objectName: "mov-boundary"; Layout.fillWidth: true; model: form.target === 0 ? [qsTr("Surface"), qsTr("Atmospheric entry \u00b7 100 km")] : [qsTr("Surface")]; currentIndex: form.target === 0 ? form.boundary : 0; onActivated: form.boundary = currentIndex; Accessible.name: qsTr("MOV boundary") }
        }
    }
    FormSection { title: qsTr("Sampling") }
    FormField {
        label: qsTr("Method")
        CustomComboBox { objectName: "mov-sampling"; Layout.fillWidth: true; model: [qsTr("Automatic"), qsTr("Admissible-region grid"), qsTr("Cobweb")]; currentIndex: form.sampling; onActivated: form.sampling = currentIndex; Accessible.name: qsTr("MOV sampling") }
    }
    GridLayout {
        Layout.fillWidth: true
        columns: width < 460 ? 1 : 2
        columnSpacing: 20
        rowSpacing: 18
        FormField { label: form.sampling === 2 ? qsTr("Directions") : qsTr("Range nodes / directions"); CustomNumberField { objectName: "mov-range-nodes"; Layout.fillWidth: true; value: form.rangeNodes; minimumValue: 8; onValueEdited: value => form.rangeNodes = value; Accessible.name: qsTr("MOV range nodes") } }
        FormField { label: form.sampling === 2 ? qsTr("Radial nodes") : qsTr("Range-rate / radial nodes"); CustomNumberField { objectName: "mov-rate-nodes"; Layout.fillWidth: true; value: form.rateNodes; minimumValue: 8; onValueEdited: value => form.rateNodes = value; Accessible.name: qsTr("MOV rate nodes") } }
        FormField { label: qsTr("Workers"); CustomNumberField { objectName: "mov-workers"; Layout.fillWidth: true; value: form.workers; minimumValue: 1; onValueEdited: value => form.workers = value; Accessible.name: qsTr("MOV workers") } }
        FormField { label: qsTr("Selection"); CustomTextField { Layout.fillWidth: true; text: "Chi <= 5; H <= 34.5"; readOnly: true } }
    }
    CustomText { visible: !form.gridValid; Layout.fillWidth: true; wrapMode: Text.Wrap; text: qsTr("Node counts: 8 to 250. Workers: 1 to 64. Whole numbers only."); font.pixelSize: 11; color: AppTheme.amber }
    FileDialog { id: fileDialog; title: qsTr("MOV observations"); nameFilters: [qsTr("Observations (*.xml *.psv *.obs *.txt *.dat)"), qsTr("All files (*)")]; onAccepted: form.mov.selectFile(selectedFile) }
}
