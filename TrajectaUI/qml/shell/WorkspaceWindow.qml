import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI

ApplicationWindow {
    id: window
    required property ApplicationSession session
    property int currentPage: 0
    property bool navigationCollapsed: false
    property bool closingAfterCancel: false
    width: 1440
    height: 900
    minimumWidth: 1024
    minimumHeight: 720
    visible: true
    title: "Trajecta"
    color: AppTheme.base
    Binding { target: AppTheme; property: "name"; value: window.session.theme }
    DialogBackdrop {
        objectName: "dialog-backdrop"
        anchors.fill: parent
        active: DialogState.activeFor(window)
        RowLayout {
            anchors.fill: parent
            spacing: 0
            NavigationRail {
                Layout.fillHeight: true
                currentPage: window.currentPage
                session: window.session
                compact: window.width <= 1200 || window.navigationCollapsed
                collapsible: window.width > 1200
                onToggleCollapsed: window.navigationCollapsed = !window.navigationCollapsed
                onNavigate: function(page) { window.currentPage = page }
                onSettingsRequested: settingsDialog.open()
                onDataSourcesRequested: dataSourcesDialog.open()
                onImportRequested: importDialog.open()
            }
            ColumnLayout {
                Layout.fillWidth: true
                Layout.fillHeight: true
                spacing: 0
                CommandBar {
                    session: window.session
                    compact: window.width <= 1200
                    pageName: [qsTr("Solution"), qsTr("Observations"), qsTr("Orbit 3D"), qsTr("Ephemerides"), qsTr("Risk Monitor"), qsTr("Close Approaches")][window.currentPage]
                    Layout.fillWidth: true
                }
                RowLayout {
                    visible: window.session.error.length > 0
                    Layout.fillWidth: true
                    Layout.margins: visible ? 12 : 0
                    CustomText { text: window.session.error; color: AppTheme.red; wrapMode: Text.Wrap; Layout.fillWidth: true }
                    CustomButton { text: qsTr("Dismiss"); onClicked: window.session.dismissError() }
                }
                SolutionHeader { visible: window.currentPage !== 2; Layout.fillWidth: true; solution: window.session.solution; viewportWidth: window.width }
                StackLayout {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    currentIndex: window.currentPage
                    SolutionPage {
                        solution: window.session.solution
                        canRecalculate: window.session.canRun
                        viewportWidth: window.width
                        onGaiaRequested: window.session.downloadGaiaFpr()
                        onObservationsRequested: function(station) {
                            window.session.solution.observations.showStation(station)
                            window.currentPage = 1
                        }
                        onImportRequested: importDialog.open()
                    }
                    ObservationsPage {
                        observations: window.session.solution.observations
                        canRecalculate: window.session.canRun
                        onImportRequested: importDialog.open()
                        onRecalculateRequested: window.session.run()
                    }
                    OrbitPage {
                        solution: window.session.solution
                        suspended: window.session.busy || DialogState.activeFor(window)
                        onImportRequested: importDialog.open()
                    }
                    EphemeridesPage {
                        ephemerides: window.session.ephemerides
                        onImportRequested: importDialog.open()
                    }
                    RiskMonitorPage {
                        monitor: window.session.risk
                        mov: window.session.mov
                        atmosphere: window.session.atmosphere
                        solution: window.session.solution
                        fitBusy: window.session.busy
                        onImportRequested: importDialog.open()
                    }
                    CloseApproachesPage {
                        solution: window.session.solution
                        approaches: window.session.approaches
                        fitBusy: window.session.busy
                        onImportRequested: importDialog.open()
                    }
                }
                Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line }
                Rectangle {
                    Layout.fillWidth: true
                    implicitHeight: 30
                    color: AppTheme.surface
                    RowLayout {
                        anchors.fill: parent
                        anchors.leftMargin: 16
                        anchors.rightMargin: 16
                        spacing: 7
                        Rectangle { width: 6; height: 6; radius: 3; color: window.session.busy ? AppTheme.amber : window.session.error.length || window.session.dataSources.selectedUnavailableCount > 0 ? AppTheme.red : !window.session.dataSources.checked || window.session.dataSources.unavailableCount > 0 || window.session.dataSources.restartRequired ? AppTheme.amber : AppTheme.accent }
                        CustomText {
                            text: window.session.busy ? window.session.progress : window.session.error.length ? qsTr("Calculation failed")
                                : !window.session.dataSources.checked ? qsTr("Checking data sources...") : window.session.dataSources.readinessMessage
                            secondary: true
                            font.pixelSize: 10
                            elide: Text.ElideRight
                            Layout.fillWidth: true
                        }
                        CustomTextLink { visible: !window.session.busy && window.session.dataSources.checked && (window.session.dataSources.unavailableCount > 0 || window.session.dataSources.restartRequired); text: qsTr("Data sources"); onClicked: dataSourcesDialog.open() }
                        CustomTextLink {
                            visible: window.session.applicationUpdate.available
                            enabled: !window.session.busy
                            text: qsTr("Update available: %1").arg(window.session.applicationUpdate.version)
                            onClicked: updateDialog.open()
                        }
                        CustomText { text: window.session.solution.available ? qsTr("Results saved to working folder") : ""; secondary: true; font.pixelSize: 10 }
                    }
                }
            }
        }
    }
    ObservationDropArea {
        anchors.fill: parent
        z: 100
        session: window.session
        enabled: !window.session.busy && !DialogState.activeFor(window)
    }
    ImportObservationsDialog { id: importDialog; session: window.session }
    CalculationDialog { id: calculation; session: window.session }
    ObservationReviewDialog { id: observationReview; session: window.session }
    CalculationSettingsDialog { id: settingsDialog; session: window.session }
    DataSourcesDialog { id: dataSourcesDialog; sources: window.session.dataSources }
    BaseDialog {
        id: updateDialog
        title: qsTr("Update Trajecta")
        acceptText: qsTr("Update")
        acceptPrimary: true
        autoAccept: false
        cancelText: qsTr("Later")
        contentItem: CustomText {
            text: window.session.applicationUpdate.error.length ? window.session.applicationUpdate.error
                : qsTr("Trajecta will close and the update tool will open.")
            wrapMode: Text.Wrap
            elide: Text.ElideNone
        }
        onAcceptRequested: {
            if (window.session.applicationUpdate.openUpdater()) {
                updateDialog.close()
                window.close()
            }
        }
    }
    BaseDialog {
        id: closeDialog
        title: qsTr("Calculation in progress")
        acceptText: qsTr("Cancel and exit")
        acceptPrimary: true
        cancelText: qsTr("Keep calculating")
        contentItem: CustomText { text: qsTr("Cancel the calculation and close when it stops?"); wrapMode: Text.Wrap; elide: Text.ElideNone }
        onAccepted: {
            window.closingAfterCancel = true
            window.session.cancel()
        }
    }
    onClosing: function(event) {
        if (session.busy) {
            event.accepted = false
            if (!closingAfterCancel) closeDialog.open()
        }
    }
    Connections {
        target: window.session
        function onRunStarted() {
            importDialog.close()
            settingsDialog.close()
            calculation.open()
        }
        function onObservationReviewReady() {
            calculation.close()
            observationReview.open()
        }
        function onChanged() { if (window.closingAfterCancel && !window.session.busy) window.close() }
    }
}
