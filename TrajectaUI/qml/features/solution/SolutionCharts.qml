import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI

Flickable {
    id: root
    objectName: "solution-charts"
    required property DiagnosticsModel diagnostics
    property ObservationTableModel observations: null
    property real pageMargin: 24
    property var stats: ({count: 0, includedCount: 0, excludedCount: 0, bins: [], excludedBins: [], histogramMax: 1})
    readonly property bool stacked: width < 1000
    clip: true
    contentWidth: width
    contentHeight: content.height + 56
    boundsBehavior: Flickable.StopAtBounds
    ScrollBar.vertical: CustomScrollBar {}

    function updateCharts() {
        if (!diagnostics) return
        stats = diagnostics.fillResidualCharts(ra.includedSeries, dec.includedSeries, ra.excludedSeries, dec.excludedSeries)
        var span = Math.max(1, (stats.jdMax - stats.timeOriginMjd) * 1.02)
        var limit = Math.ceil(Math.max(stats.draMax, stats.ddecMax) * 1.15 * 2) / 2
        ra.timeOriginMjd = dec.timeOriginMjd = stats.timeOriginMjd
        ra.timeOriginLabel = dec.timeOriginLabel = stats.timeOriginLabel
        ra.bounds = dec.bounds = Qt.rect(0, -limit, span, 2 * limit)
        ra.refreshSelection()
        dec.refreshSelection()
    }
    Component.onCompleted: updateCharts()
    Connections {
        target: root.diagnostics
        function onChanged() { root.updateCharts() }
    }
    ColumnLayout {
        id: content
        x: root.pageMargin
        y: 24
        width: Math.max(1, root.width - 2 * x)
        height: Math.max(implicitHeight, root.height - 56)
        spacing: 20
        RowLayout {
            Layout.fillWidth: true
            spacing: 16
            CustomText { text: qsTr("Fit diagnostics"); font.pixelSize: 18; font.weight: Font.DemiBold; Layout.fillWidth: true }
            CustomText { text: qsTr("%1 optical observations").arg(root.stats.count); secondary: true; font.pixelSize: 11 }
        }
        GridLayout {
            visible: root.stats.count > 0
            Layout.fillWidth: true
            Layout.fillHeight: true
            columns: root.stacked ? 1 : 2
            uniformCellWidths: true
            uniformCellHeights: true
            columnSpacing: 24
            rowSpacing: 24
            ColumnLayout {
                Layout.fillWidth: true
                Layout.fillHeight: true
                Layout.minimumWidth: 0
                spacing: 10
                CustomText { text: qsTr("Right ascension"); font.pixelSize: 14; font.weight: Font.DemiBold }
                Rectangle {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    Layout.minimumHeight: 232
                    Layout.preferredHeight: 232
                    color: AppTheme.surface; radius: 4; border.color: AppTheme.line
                    ResidualTimeChart {
                        id: ra
                        objectName: "ra-time-chart"
                        anchors.fill: parent; anchors.margins: 1
                        observations: root.observations
                    }
                }
            }
            ColumnLayout {
                Layout.fillWidth: true
                Layout.fillHeight: true
                Layout.minimumWidth: 0
                spacing: 10
                CustomText { text: qsTr("Declination"); font.pixelSize: 14; font.weight: Font.DemiBold }
                Rectangle {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    Layout.minimumHeight: 232
                    Layout.preferredHeight: 232
                    color: AppTheme.surface; radius: 4; border.color: AppTheme.line
                    ResidualTimeChart {
                        id: dec
                        objectName: "dec-time-chart"
                        anchors.fill: parent; anchors.margins: 1
                        observations: root.observations
                        declination: true
                    }
                }
            }
            ColumnLayout {
                Layout.columnSpan: root.stacked ? 1 : 2
                Layout.fillWidth: true
                Layout.fillHeight: true
                Layout.minimumWidth: 0
                spacing: 10
                CustomText { text: qsTr("Residual distribution"); font.pixelSize: 14; font.weight: Font.DemiBold }
                Rectangle {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    Layout.minimumHeight: 232
                    Layout.preferredHeight: 232
                    color: AppTheme.surface; radius: 4; border.color: AppTheme.line
                    ResidualHistogram {
                        objectName: "residual-histogram"
                        anchors.fill: parent; anchors.margins: 1
                        includedBins: root.stats.bins
                        excludedBins: root.stats.excludedBins
                        maximum: root.stats.histogramMax
                    }
                }
            }
        }
        RowLayout {
            visible: root.stats.count > 0
            Layout.fillWidth: true
            spacing: 20
            PlotLegendItem { label: qsTr("Included: %1").arg(root.stats.includedCount); color: AppTheme.accent }
            PlotLegendItem { label: qsTr("Excluded: %1").arg(root.stats.excludedCount); color: AppTheme.red }
            Item { Layout.fillWidth: true }
        }
        CustomText {
            objectName: "charts-empty"
            visible: root.stats.count === 0
            Layout.fillWidth: true
            Layout.topMargin: 60
            text: qsTr("No optical residuals available")
            horizontalAlignment: Text.AlignHCenter
            secondary: true
        }
    }
}
