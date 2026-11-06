pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

ColumnLayout {
    id: panel
    required property RiskMonitorViewModel monitor
    required property SolutionViewModel solution
    property bool fitBusy: false
    signal importRequested()
    spacing: 20
    CustomText { text: qsTr("Source data"); font.pixelSize: 11; font.weight: Font.DemiBold }
    RowLayout {
        visible: panel.solution.available
        Layout.fillWidth: true
        spacing: 8
        Rectangle { implicitWidth: 6; implicitHeight: 6; radius: 3; color: panel.solution.converged ? AppTheme.accent : AppTheme.amber }
        CustomText {
            Layout.fillWidth: true
            text: panel.solution.status
            color: panel.solution.converged ? AppTheme.accent : AppTheme.amber
            font.pixelSize: 11
            wrapMode: Text.Wrap
        }
    }
    KeyValueList {
        visible: panel.monitor.available
        Layout.fillWidth: true
        rows: [[qsTr("Nominal"), panel.solution.objectTitle],
               [qsTr("Full covariance"), panel.monitor.covarianceSize],
               [qsTr("Fit epoch, MJD"), panel.monitor.fitEpoch],
               [qsTr("Included observations"), String(panel.monitor.observationCount)],
               [qsTr("Optical / radar"), (panel.monitor.observationCount - panel.monitor.radarCount) + " / " + panel.monitor.radarCount],
               [qsTr("Non-gravitational"), panel.solution.nonGravEnabled ? panel.solution.nonGravModel : qsTr("Off")],
               [qsTr("LOV"), qsTr("Not required")],
               [qsTr("Date scale"), "UTC"],
               [qsTr("Interval"), panel.monitor.windowDuration || qsTr("Not set")]]
    }
    CustomText {
        objectName: "risk-source-error"
        Layout.fillWidth: true
        visible: text.length > 0
        text: panel.monitor.sourceError
        color: AppTheme.amber
        font.pixelSize: 11
        wrapMode: Text.Wrap
        elide: Text.ElideNone
    }
    CustomButton {
        visible: !panel.solution.available
        text: qsTr("Open observations")
        iconName: "folder-open"
        enabled: !panel.fitBusy
        onClicked: panel.importRequested()
    }
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line }
    CustomText { text: qsTr("Route"); font.pixelSize: 11; font.weight: Font.DemiBold }
    Repeater {
        model: panel.monitor.method === 0
            ? [qsTr("Covariance samples"), qsTr("Interval propagation"), panel.monitor.intervalMode === 1 ? qsTr("Encounter catalogue") : qsTr("Boundary intersections"), qsTr("Archive and HTML")]
            : [qsTr("Synthetic observations"), qsTr("Refit and propagation"), panel.monitor.intervalMode === 1 ? qsTr("Encounter catalogue") : qsTr("Boundary intersections"), qsTr("Archive and HTML")]
        RowLayout {
            id: routeRow
            required property string modelData
            required property int index
            Layout.fillWidth: true
            spacing: 16
            CustomText { text: "0" + (routeRow.index + 1); numeric: true; color: AppTheme.accent; font.pixelSize: 10 }
            CustomText { text: routeRow.modelData; secondary: true; font.pixelSize: 11; Layout.fillWidth: true; wrapMode: Text.Wrap; elide: Text.ElideNone }
        }
    }
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line }
    CustomText { text: qsTr("Outcomes"); font.pixelSize: 11; font.weight: Font.DemiBold }
    CustomText {
        text: qsTr("Boundary crossing / Miss / Calculation failure")
        secondary: true
        font.pixelSize: 11
        Layout.fillWidth: true
        wrapMode: Text.Wrap
        elide: Text.ElideNone
    }
}
