import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

ColumnLayout {
    id: summary
    required property SolutionViewModel solution
    property real viewportWidth: 1440
    signal observationsRequested(string station)
    spacing: 0
    GridLayout {
        id: metrics
        Layout.fillWidth: true
        Layout.bottomMargin: 21
        columns: width < 1000 ? 3 : 6
        columnSpacing: summary.viewportWidth <= 1200 ? 26 : summary.viewportWidth >= 1700 ? 60 : 38
        rowSpacing: 18
        MetricValue { objectName: "solution-optical-rms"; label: summary.solution.observations.radarCount ? qsTr("Optical RMS") : qsTr("RMS"); value: UI.number(summary.solution.rms, 5); context: qsTr("Unweighted · arcsec"); valueSize: summary.viewportWidth <= 1200 ? 18 : 21 }
        MetricValue { objectName: "solution-normalized-rms"; label: qsTr("Normalized RMS"); value: UI.number(summary.solution.normalizedRms, 5); context: qsTr("Weighted · dimensionless"); valueSize: summary.viewportWidth <= 1200 ? 18 : 21 }
        MetricValue { label: qsTr("Observations"); value: String(summary.solution.included); suffix: "/ " + summary.solution.total; context: qsTr("%1 stations").arg(summary.solution.stationCount); valueSize: summary.viewportWidth <= 1200 ? 18 : 21 }
        MetricValue { label: qsTr("Observation arc"); value: UI.number(summary.solution.arcDays, 2); suffix: qsTr("days"); context: summary.solution.arcRange; valueSize: summary.viewportWidth <= 1200 ? 18 : 21 }
        Item { visible: metrics.columns === 6; Layout.fillWidth: true }
        MetricValue { label: qsTr("Solution epoch · TDB"); value: summary.solution.epochMjd; context: "MJD · ECLM J2000"; valueSize: 17 }
    }
    Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line; Layout.bottomMargin: 24 }
    RowLayout {
        id: columns
        Layout.fillWidth: true
        spacing: summary.viewportWidth <= 1200 ? 22 : summary.viewportWidth >= 1700 ? 40 : 30
        readonly property real orbitWidth: Math.max(summary.viewportWidth <= 1200 ? 280 : 300, (width - spacing) / (summary.viewportWidth <= 1200 ? 2.3 : summary.viewportWidth >= 1700 ? 2.8 : 2.6))
        ResidualSection {
            solution: summary.solution
            viewportWidth: summary.viewportWidth
            Layout.fillWidth: true
            Layout.alignment: Qt.AlignTop
        }
        Item {
            Layout.preferredWidth: columns.orbitWidth
            Layout.minimumWidth: columns.orbitWidth
            Layout.preferredHeight: elements.implicitHeight
            Layout.alignment: Qt.AlignTop
            Rectangle { width: 1; height: parent.height; color: AppTheme.line }
            OrbitElementsTable {
                id: elements
                x: summary.viewportWidth <= 1200 ? 20 : 28
                width: parent.width - x
                solution: summary.solution
            }
        }
    }
    SelectedObservationStrip {
        Layout.fillWidth: true
        Layout.topMargin: 20
        observations: summary.solution.observations
        onOpenRequested: observationDetails.open()
    }
    StationDiagnostics {
        Layout.fillWidth: true
        Layout.topMargin: summary.viewportWidth >= 1700 ? 32 : 27
        statistics: summary.solution.observations.stationResiduals
        onObservationsRequested: function(station) { summary.observationsRequested(station) }
    }
    ObservationDetailsDialog { id: observationDetails; observations: summary.solution.observations }
}
