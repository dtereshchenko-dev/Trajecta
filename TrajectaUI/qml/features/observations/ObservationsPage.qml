import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

BasePage {
    id: page
    objectName: "observations-page"
    required property ObservationTableModel observations
    property bool canRecalculate: false
    property bool inspectorOpen: false
    readonly property bool inlineInspector: width >= 1250
    readonly property real pageMargin: UI.pageMargin(width)
    signal importRequested()
    signal recalculateRequested()
    title: qsTr("Observations and weights")
    padding: 0
    function openDetails() {
        if (observations.selectedSourceRow < 0) return
        if (inlineInspector) inspectorOpen = true
        else details.open()
    }
    header: Rectangle {
        implicitHeight: 48
        color: AppTheme.base
        RowLayout {
            anchors.fill: parent
            anchors.leftMargin: page.pageMargin
            anchors.rightMargin: page.pageMargin
            spacing: 14
            CustomText { text: page.title; font.pixelSize: 14; font.weight: Font.DemiBold; Layout.fillWidth: true }
            CustomText { text: page.observations.sourceFormat; secondary: true; font.pixelSize: 11 }
            CustomButton {
                objectName: "observations-inspector-toggle"
                iconName: "info"
                flat: true
                checked: page.inlineInspector && page.inspectorOpen
                enabled: page.observations.selectedSourceRow >= 0
                accessibleLabel: qsTr("Observation details")
                onClicked: {
                    if (page.inlineInspector && page.inspectorOpen) page.inspectorOpen = false
                    else page.openDetails()
                }
            }
        }
        Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line }
    }
    contentItem: Item {
        ColumnLayout {
            anchors.fill: parent
            anchors.margins: page.pageMargin
            spacing: 18
            CustomSegmentedControl {
                objectName: "observations-type"
                visible: page.observations.radarCount > 0
                options: [qsTr("Optical (%1)").arg(page.observations.totalCount - page.observations.radarCount),
                    qsTr("Radar (%1)").arg(page.observations.radarCount)]
                currentIndex: page.observations.radarMode ? 1 : 0
                onActivated: index => page.observations.radarMode = index === 1
            }
            ObservationFilters { observations: page.observations; Layout.fillWidth: true; visible: page.observations.totalCount > 0 }
            RowLayout {
                Layout.fillWidth: true
                Layout.fillHeight: true
                spacing: 24
                Item {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    Layout.minimumWidth: 300
                    ObservationsTable {
                        id: table
                        objectName: "observations-table"
                        anchors.fill: parent
                        observations: page.observations
                        visible: page.observations.count > 0
                        onObservationActivated: page.openDetails()
                    }
                    ColumnLayout {
                        anchors.centerIn: parent
                        spacing: 16
                        visible: page.observations.count === 0
                        CustomIcon { Layout.alignment: Qt.AlignHCenter; name: "table-2"; size: 28; color: AppTheme.muted }
                        CustomText { objectName: "observations-empty"; text: page.observations.totalCount ? qsTr("No matching observations") : qsTr("No observations"); Layout.alignment: Qt.AlignHCenter }
                        CustomButton {
                            objectName: "observations-empty-action"
                            Layout.alignment: Qt.AlignHCenter
                            text: page.observations.totalCount ? qsTr("Clear filters") : qsTr("Open observations")
                            iconName: page.observations.totalCount ? "x" : "folder-open"
                            onClicked: { if (page.observations.totalCount) page.observations.clearFilters(); else page.importRequested() }
                        }
                    }
                }
                Rectangle { implicitWidth: 1; color: AppTheme.line; Layout.fillHeight: true; visible: inspector.visible }
                ObservationInspector {
                    id: inspector
                    objectName: "observations-inspector"
                    observations: page.observations
                    visible: page.inlineInspector && page.inspectorOpen && page.observations.selectedSourceRow >= 0
                    Layout.minimumWidth: 300
                    Layout.preferredWidth: 300
                    Layout.maximumWidth: 300
                    Layout.fillHeight: true
                    onCloseRequested: page.inspectorOpen = false
                }
            }
            RowLayout {
                Layout.fillWidth: true
                visible: page.observations.pendingCount > 0 || page.observations.manualCount > 0
                spacing: 10
                CustomText {
                    objectName: "observations-pending-count"
                    text: page.observations.pendingCount ? qsTr("Pending changes: %1").arg(page.observations.pendingCount)
                        : qsTr("Manual choices: %1").arg(page.observations.manualCount)
                    color: page.observations.pendingCount ? AppTheme.amber : AppTheme.muted
                    font.pixelSize: 11
                    Layout.fillWidth: true
                }
                CustomButton {
                    objectName: "observations-reset-manual"
                    text: qsTr("Reset manual selection")
                    enabled: page.observations.editingEnabled && page.observations.manualCount > 0
                    onClicked: page.observations.resetManualSelection()
                }
                CustomButton {
                    objectName: "observations-discard"
                    text: qsTr("Discard")
                    enabled: page.observations.editingEnabled && page.observations.pendingCount > 0
                    onClicked: page.observations.discardEdits()
                }
                CustomButton {
                    objectName: "observations-recalculate"
                    text: qsTr("Re-solve")
                    iconName: "play"
                    primary: true
                    enabled: page.canRecalculate && page.observations.pendingCount > 0
                    onClicked: page.recalculateRequested()
                }
            }
            RowLayout {
                Layout.fillWidth: true
                visible: page.observations.totalCount > 0
                spacing: 18
                PlotLegendItem { label: qsTr("Included: %1").arg(page.observations.stationResiduals.included); color: AppTheme.accent }
                PlotLegendItem { label: qsTr("Excluded: %1").arg(page.observations.stationResiduals.excluded); color: AppTheme.red }
                Item { Layout.fillWidth: true }
                CustomText { text: qsTr("Current solution"); secondary: true; font.pixelSize: 10 }
            }
        }
    }
    ObservationDetailsDialog { id: details; objectName: "observations-page-details"; observations: page.observations }
}
