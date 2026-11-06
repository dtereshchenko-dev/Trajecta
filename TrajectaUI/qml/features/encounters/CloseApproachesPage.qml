pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

BasePage {
    id: page
    objectName: "close-approaches-page"
    required property SolutionViewModel solution
    required property CloseApproachesViewModel approaches
    property bool fitBusy: false
    property int currentView: 0
    property real maximumDistanceAu: 0.2
    property bool automaticDistance: true
    readonly property real distanceLimit: automaticDistance ? 0 : maximumDistanceAu > 0 ? maximumDistanceAu : NaN
    readonly property string formError: approaches.available
        ? approaches.validationError(interval.startUtc, interval.endUtc, distanceLimit)
        : qsTr("A converged N-body solution is required.")
    readonly property int distanceUnitIndex: approaches.table.distanceUnit
    readonly property int velocityUnitIndex: approaches.table.velocityUnit
    readonly property var distanceUnits: ["au", "km", "LD", "mi", "Re"]
    readonly property var velocityUnits: ["km/s", "mi/h", "au/d"]
    readonly property string distanceUnit: distanceUnits[distanceUnitIndex]
    readonly property string velocityUnit: velocityUnits[velocityUnitIndex]
    readonly property real pageMargin: UI.pageMargin(width)
    signal importRequested()
    title: qsTr("Close Approaches")
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
                objectName: "approaches-view"
                options: [qsTr("New calculation"), qsTr("Events")]
                currentIndex: page.currentView
                onActivated: index => page.currentView = index
            }
            CustomButton {
                objectName: "approaches-run"
                iconName: page.approaches.busy ? "square" : "search"
                text: page.approaches.busy ? qsTr("Cancel") : qsTr("Search")
                primary: !page.approaches.busy
                enabled: page.approaches.busy ? !page.approaches.cancelling : page.approaches.canRun && !page.formError.length
                accessibleLabel: text
                onClicked: page.approaches.busy ? page.approaches.cancel()
                    : page.approaches.run(interval.startUtc, interval.endUtc, page.distanceLimit)
            }
        }
        Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line }
    }
    contentItem: ScrollView {
        id: scroll
        objectName: "approaches-scroll"
        contentWidth: availableWidth
        clip: true
        ScrollBar.horizontal.policy: ScrollBar.AlwaysOff
        ScrollBar.vertical: CustomScrollBar {
            objectName: "approaches-scrollbar"
            parent: scroll
            x: scroll.width - width
            y: scroll.topPadding
            height: scroll.availableHeight
        }
        ColumnLayout {
            width: scroll.availableWidth
            spacing: 0
            ColumnLayout {
                Layout.fillWidth: true
                Layout.leftMargin: page.pageMargin
                Layout.rightMargin: page.pageMargin
                Layout.topMargin: visible ? 16 : 0
                visible: page.approaches.busy || page.approaches.error.length > 0
                CustomText {
                    Layout.fillWidth: true
                    text: page.approaches.busy ? page.approaches.progress : page.approaches.error
                    color: page.approaches.error.length ? AppTheme.red : AppTheme.muted
                    wrapMode: Text.Wrap
                    font.pixelSize: 11
                }
                CustomProgressBar { Layout.fillWidth: true; visible: page.approaches.busy; indeterminate: true }
            }
            GridLayout {
                id: workbench
                visible: page.currentView === 0
                Layout.fillWidth: true
                Layout.margins: page.pageMargin
                columns: width < 900 ? 1 : 3
                columnSpacing: 28
                rowSpacing: 28
                ColumnLayout {
                    Layout.fillWidth: true
                    Layout.preferredWidth: 740
                    Layout.alignment: Qt.AlignTop
                    spacing: 20
                    RowLayout {
                        Layout.fillWidth: true
                        CustomText { text: qsTr("Nominal encounters"); font.pixelSize: 18; font.weight: Font.DemiBold; Layout.fillWidth: true }
                        CustomText { text: qsTr("N-body"); color: AppTheme.accent; font.pixelSize: 10 }
                    }
                    RowLayout {
                        Layout.fillWidth: true
                        Layout.topMargin: -10
                        spacing: 9
                        CustomIcon { name: "orbit"; color: AppTheme.muted; size: 16 }
                        CustomText { text: qsTr("All planets + Moon"); secondary: true; font.pixelSize: 11 }
                    }
                    FormSection { title: qsTr("Search interval") }
                    SearchIntervalForm {
                        id: interval
                        objectName: "approaches-interval"
                        fieldPrefix: "approaches"
                        seedUtc: page.solution.calculationStartUtc
                        sourceIdentity: page.solution.designation
                        Layout.fillWidth: true
                        enabled: !page.approaches.busy
                    }
                    FormSection { title: qsTr("Encounter filter") }
                    GridLayout {
                        Layout.fillWidth: true
                        columns: width < 460 ? 1 : 2
                        columnSpacing: 20
                        enabled: !page.approaches.busy
                        FormField {
                            label: qsTr("Reporting distance")
                            CustomComboBox {
                                objectName: "approaches-distance-mode"
                                Layout.fillWidth: true
                                model: [qsTr("Automatic by body"), qsTr("Custom maximum")]
                                currentIndex: page.automaticDistance ? 0 : 1
                                onActivated: page.automaticDistance = currentIndex === 0
                                Accessible.name: qsTr("Reporting distance policy")
                            }
                        }
                        FormField {
                            visible: !page.automaticDistance
                            label: qsTr("Maximum distance, AU")
                            CustomNumberField {
                                objectName: "approaches-distance"
                                Layout.fillWidth: true
                                value: page.maximumDistanceAu
                                minimumValue: 0
                                exclusiveMinimum: true
                                Accessible.name: qsTr("Maximum body-centred distance in AU")
                                onValueEdited: value => page.maximumDistanceAu = value
                            }
                        }
                        FormField {
                            visible: !page.automaticDistance
                            label: qsTr("Distance, km")
                            CustomTextField {
                                objectName: "approaches-distance-km"
                                Layout.fillWidth: true
                                readOnly: true
                                text: page.maximumDistanceAu > 0 ? UI.groupedNumber(page.maximumDistanceAu * 149597870.7, 0) : ""
                                font.family: UI.monoFont
                                Accessible.name: qsTr("Maximum distance in kilometres")
                            }
                        }
                    }
                    KeyValueList {
                        visible: page.automaticDistance
                        Layout.fillWidth: true
                        rows: [[qsTr("Mercury"), "0.10 au"], [qsTr("Venus"), "0.15 au"],
                               [qsTr("Earth"), "0.50 au"], [qsTr("Mars"), "0.10 au"],
                               [qsTr("Moon"), "0.01 au"], [qsTr("Jupiter / Saturn"), "2.00 au"],
                               [qsTr("Uranus / Neptune"), "2.00 au"]]
                    }
                    CustomText {
                        Layout.fillWidth: true
                        visible: page.formError.length > 0
                        text: page.formError
                        color: AppTheme.amber
                        font.pixelSize: 11
                        wrapMode: Text.Wrap
                    }
                    FormSection { title: qsTr("Output") }
                    RowLayout {
                        Layout.fillWidth: true
                        spacing: 12
                        CustomIcon { name: "table-2"; color: AppTheme.blue; size: 18 }
                        ColumnLayout {
                            Layout.fillWidth: true
                            spacing: 7
                            CustomText { text: qsTr("Encounter catalogue"); font.pixelSize: 11; font.weight: Font.DemiBold; Layout.fillWidth: true }
                            CustomText { text: qsTr("TDB / Body / Distance bounds / Velocities"); font.pixelSize: 10; secondary: true; Layout.fillWidth: true; wrapMode: Text.Wrap }
                        }
                    }
                }
                Rectangle { visible: workbench.columns === 3; Layout.preferredWidth: 1; Layout.fillHeight: true; color: AppTheme.line }
                ColumnLayout {
                    Layout.fillWidth: true
                    Layout.preferredWidth: 280
                    Layout.alignment: Qt.AlignTop
                    spacing: 20
                    CustomText { text: qsTr("Source data"); font.pixelSize: 11; font.weight: Font.DemiBold }
                    KeyValueList {
                        Layout.fillWidth: true
                        rows: [[qsTr("Nominal"), page.solution.available ? page.solution.objectTitle : qsTr("Not available")],
                               [qsTr("Solution"), page.solution.status],
                               [qsTr("Fit epoch, UTC"), page.solution.fitEpochUtc || "\u2014"],
                               [qsTr("Non-gravitational"), page.solution.available ? page.solution.nonGravEnabled ? page.solution.nonGravModel : qsTr("Off") : "\u2014"],
                               [qsTr("Covariance"), page.solution.covariance.dimension >= 6 ? qsTr("Available") : qsTr("Not available")],
                               [qsTr("Distance bounds, 3\u03c3"), qsTr("Local linear")],
                               [qsTr("Time uncertainty, 3\u03c3"), qsTr("Local linear")]]
                    }
                    CustomButton {
                        visible: !page.solution.available
                        text: qsTr("Open observations")
                        iconName: "folder-open"
                        enabled: !page.fitBusy
                        onClicked: page.importRequested()
                    }
                    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line }
                    CustomText { text: qsTr("Search scope"); font.pixelSize: 11; font.weight: Font.DemiBold }
                    KeyValueList {
                        Layout.fillWidth: true
                        rows: [[qsTr("Bodies"), qsTr("8 planets + Moon")],
                               [qsTr("Trajectory"), qsTr("Nominal")],
                               [qsTr("Distance origin"), qsTr("Body centre")],
                               [qsTr("Impact probability"), qsTr("Not calculated")]]
                    }
                }
            }
            ColumnLayout {
                visible: page.currentView === 1
                Layout.fillWidth: true
                Layout.margins: page.pageMargin
                spacing: 20
                SectionHeading { title: qsTr("Encounters"); detail: page.approaches.hasResult ? qsTr("%1 events").arg(page.approaches.table.count) : qsTr("Not calculated"); Layout.fillWidth: true }
                RowLayout {
                    Layout.fillWidth: true
                    CustomText {
                        Layout.fillWidth: true
                        text: page.approaches.summary
                        secondary: true
                        font.pixelSize: 11
                        wrapMode: Text.Wrap
                    }
                    CustomButton {
                        iconName: "folder-open"
                        accessibleLabel: qsTr("Open encounter archive")
                        enabled: page.approaches.archivePath.length > 0
                        onClicked: page.approaches.openArchive()
                    }
                }
                GridLayout {
                    Layout.fillWidth: true
                    Layout.maximumWidth: 500
                    columns: width < 400 ? 1 : 2
                    columnSpacing: 20
                    rowSpacing: 16
                    FormField {
                        label: qsTr("Distance units")
                        CustomComboBox {
                            objectName: "approaches-distance-unit"
                            Layout.fillWidth: true
                            model: page.distanceUnits
                            currentIndex: page.distanceUnitIndex
                            Accessible.name: qsTr("Encounter distance units")
                            onActivated: page.approaches.table.distanceUnit = currentIndex
                        }
                    }
                    FormField {
                        label: qsTr("Velocity units")
                        CustomComboBox {
                            objectName: "approaches-velocity-unit"
                            Layout.fillWidth: true
                            model: page.velocityUnits
                            currentIndex: page.velocityUnitIndex
                            Accessible.name: qsTr("Encounter velocity units")
                            onActivated: page.approaches.table.velocityUnit = currentIndex
                        }
                    }
                }
                CloseApproachTable {
                    objectName: "approaches-events"
                    Layout.fillWidth: true
                    Layout.preferredHeight: Math.max(240, page.height - 330)
                    results: page.approaches.table
                    emptyText: page.approaches.busy ? qsTr("Searching encounters...")
                        : page.approaches.hasResult ? qsTr("No encounters within the reporting limits") : qsTr("No calculated encounters")
                }
            }
        }
    }
    Connections {
        target: page.approaches
        function onCompleted() { page.currentView = 1 }
    }
}
