pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

BasePage {
    id: page
    objectName: "risk-monitor-page"
    required property RiskMonitorViewModel monitor
    required property MovViewModel mov
    required property AtmosphereViewModel atmosphere
    required property SolutionViewModel solution
    property bool fitBusy: false
    property int currentView: 0
    property bool showingAtmosphere: false
    property bool showingMov: false
    onShowingAtmosphereChanged: { if (showingAtmosphere) showingMov = false }
    onShowingMovChanged: { if (showingMov) showingAtmosphere = false }
    readonly property bool jobBusy: monitor.busy || atmosphere.busy || mov.busy
    readonly property real pageMargin: UI.pageMargin(width)
    signal importRequested()
    title: qsTr("Risk Monitor")
    padding: 0
    header: Rectangle {
        color: AppTheme.base
        implicitHeight: 56
        RowLayout {
            anchors.fill: parent
            anchors.leftMargin: page.pageMargin
            anchors.rightMargin: page.pageMargin
            spacing: 16
            CustomText { text: page.title; font.pixelSize: 14; font.weight: Font.DemiBold; Layout.fillWidth: true }
            CustomSegmentedControl {
                objectName: "risk-view"
                options: [qsTr("New calculation"), qsTr("Runs")]
                currentIndex: page.currentView
                onActivated: index => page.currentView = index
            }
            CustomButton {
                objectName: "risk-run"
                iconName: page.jobBusy ? "square" : "play"
                text: page.jobBusy ? qsTr("Cancel") : qsTr("Run")
                primary: !page.jobBusy
                enabled: page.mov.busy ? !page.mov.cancelling : page.atmosphere.busy ? !page.atmosphere.cancelling : page.monitor.busy ? !page.monitor.cancelling
                    : page.showingMov ? page.mov.canRun && movForm.gridValid : page.showingAtmosphere ? page.atmosphere.canRun : page.monitor.canRun
                accessibleLabel: text
                onClicked: {
                    if (page.mov.busy) page.mov.cancel()
                    else if (page.atmosphere.busy) page.atmosphere.cancel()
                    else if (page.monitor.busy) page.monitor.cancel()
                    else if (page.showingAtmosphere) page.atmosphere.run()
                    else if (page.showingMov) page.mov.run(movForm.options)
                    else page.monitor.run()
                }
            }
        }
        Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line }
    }
    contentItem: ScrollView {
        id: scroll
        objectName: "risk-scroll"
        clip: true
        contentWidth: availableWidth
        ScrollBar.horizontal.policy: ScrollBar.AlwaysOff
        ScrollBar.vertical: CustomScrollBar {
            objectName: "risk-scrollbar"
            parent: scroll
            x: scroll.width - width
            y: scroll.topPadding
            height: scroll.availableHeight
        }
        ColumnLayout {
            width: scroll.availableWidth
            spacing: 0
            ColumnLayout {
                visible: page.showingMov && (page.mov.busy || page.mov.error.length > 0)
                Layout.fillWidth: true
                Layout.margins: visible ? page.pageMargin : 0
                CustomText { Layout.fillWidth: true; text: page.mov.error.length ? page.mov.error : page.mov.progress; color: page.mov.error.length ? AppTheme.red : AppTheme.muted; wrapMode: Text.Wrap; font.pixelSize: 11 }
                CustomProgressBar { Layout.fillWidth: true; visible: page.mov.busy; indeterminate: true }
            }
            RowLayout {
                visible: page.currentView === 0
                Layout.fillWidth: true
                Layout.leftMargin: page.pageMargin
                Layout.rightMargin: page.pageMargin
                Layout.topMargin: 20
                spacing: 0
                Repeater {
                    model: [{method: 0, name: qsTr("Covariance MC"), icon: "activity"}, {method: 1, name: qsTr("Observation MC"), icon: "scan-line"}, {method: 3, name: "MOV", icon: "table-2"}, {method: 2, name: qsTr("Atmospheric MC"), icon: "cloud-sun"}]
                    CustomTab {
                        id: methodTab
                        required property var modelData
                        objectName: "risk-method-" + modelData.method
                        text: modelData.name
                        checked: modelData.method === 3 ? page.showingMov : modelData.method === 2 ? page.showingAtmosphere : !page.showingAtmosphere && !page.showingMov && page.monitor.method === modelData.method
                        implicitWidth: methodContent.implicitWidth + 28
                        implicitHeight: 46
                        contentItem: RowLayout {
                            id: methodContent
                            spacing: 10
                            CustomIcon { name: methodTab.modelData.icon; size: 16; color: methodTab.checked ? AppTheme.accent : AppTheme.muted }
                            CustomText { text: methodTab.text; font.pixelSize: 11; color: methodTab.checked ? AppTheme.accent : AppTheme.muted }
                        }
                        leftPadding: 14
                        rightPadding: 14
                        background: Rectangle {
                            color: methodTab.checked ? AppTheme.accentTint : "transparent"
                            Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 2; color: methodTab.checked || methodTab.activeFocus ? AppTheme.accent : "transparent" }
                        }
                        onClicked: {
                            page.showingAtmosphere = modelData.method === 2
                            page.showingMov = modelData.method === 3
                            if (modelData.method < 2) page.monitor.method = modelData.method
                        }
                    }
                }
                Item { Layout.fillWidth: true }
            }
            Rectangle {
                visible: page.currentView === 0
                Layout.fillWidth: true
                Layout.leftMargin: page.pageMargin
                Layout.rightMargin: page.pageMargin
                implicitHeight: 1
                color: AppTheme.line
            }
            GridLayout {
                id: workbench
                visible: page.currentView === 0
                readonly property bool compact: width < 900
                Layout.fillWidth: true
                Layout.leftMargin: page.pageMargin
                Layout.rightMargin: page.pageMargin
                Layout.topMargin: 24
                Layout.bottomMargin: page.pageMargin
                columns: compact ? 1 : 3
                columnSpacing: 28
                rowSpacing: 24
                RiskCalculationForm {
                    visible: !page.showingAtmosphere && !page.showingMov
                    monitor: page.monitor
                    enabled: !page.monitor.busy
                    Layout.fillWidth: true
                    Layout.preferredWidth: 740
                    Layout.alignment: Qt.AlignTop
                }
                AtmosphereForm {
                    atmosphere: page.atmosphere
                    visible: page.showingAtmosphere
                    enabled: !page.atmosphere.busy
                    Layout.fillWidth: true
                    Layout.preferredWidth: 740
                    Layout.alignment: Qt.AlignTop
                }
                MovCalculationForm {
                    id: movForm
                    mov: page.mov
                    enabled: !page.fitBusy
                    visible: page.showingMov
                    Layout.fillWidth: true
                    Layout.preferredWidth: 740
                    Layout.alignment: Qt.AlignTop
                }
                Rectangle {
                    visible: !workbench.compact
                    Layout.preferredWidth: 1
                    Layout.fillHeight: true
                    color: AppTheme.line
                }
                RiskSourcePanel {
                    visible: !page.showingAtmosphere && !page.showingMov
                    monitor: page.monitor
                    solution: page.solution
                    fitBusy: page.fitBusy
                    Layout.fillWidth: true
                    Layout.preferredWidth: 280
                    Layout.alignment: Qt.AlignTop
                    onImportRequested: page.importRequested()
                }
                AtmosphereSourcePanel {
                    atmosphere: page.atmosphere
                    visible: page.showingAtmosphere
                    Layout.fillWidth: true
                    Layout.preferredWidth: 280
                    Layout.alignment: Qt.AlignTop
                }
                MovSourcePanel {
                    mov: page.mov
                    form: movForm
                    visible: page.showingMov
                    Layout.fillWidth: true
                    Layout.preferredWidth: 280
                    Layout.alignment: Qt.AlignTop
                }
            }
            RiskRunView {
                objectName: "risk-runs"
                monitor: page.monitor
                visible: page.currentView === 1 && !page.showingAtmosphere && !page.showingMov
                Layout.fillWidth: true
                Layout.margins: page.pageMargin
                tableHeight: Math.max(300, page.height - 340)
            }
            AtmosphereRunView {
                atmosphere: page.atmosphere
                visible: page.currentView === 1 && page.showingAtmosphere
                Layout.fillWidth: true
                Layout.margins: page.pageMargin
            }
            MovRunView {
                mov: page.mov
                visible: page.currentView === 1 && page.showingMov
                Layout.fillWidth: true
                Layout.margins: page.pageMargin
            }
        }
    }
    Connections {
        target: page.mov
        function onRunStarted() { page.showingMov = true; page.currentView = 1 }
        function onRunFailed() { if (page.showingMov) page.currentView = 0 }
    }
    Connections {
        target: page.monitor
        function onRunStarted() { page.showingMov = false; page.showingAtmosphere = false; page.currentView = 1 }
    }
    Connections {
        target: page.atmosphere
        function onRunStarted() { page.showingAtmosphere = true; page.currentView = 1 }
    }
}
