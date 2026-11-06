import QtQuick
import QtCharts
import TrajectaUI

ChartView {
    id: chart
    property var includedBins: []
    property var excludedBins: []
    property real maximum: 1
    property int hoveredBin: -1
    readonly property real binWidth: maximum / Math.max(1, includedBins.length)
    antialiasing: true
    backgroundColor: AppTheme.surface
    plotAreaColor: AppTheme.surface
    plotArea: Qt.rect(0, 0, 0, 0)
    backgroundRoundness: 0
    legend.visible: false
    margins { left: 0; right: 0; top: 0; bottom: 0 }
    clip: true

    function rebuild() {
        if (!includedEdge || !excludedBottom || !excludedTop) return
        hoveredBin = -1
        includedEdge.clear()
        excludedBottom.clear()
        excludedTop.clear()
        for (var i = 0; i < includedBins.length; ++i) {
            var x0 = i * binWidth + binWidth * 0.025
            var x1 = (i + 1) * binWidth - binWidth * 0.025
            var low = includedBins[i]
            var high = low + (excludedBins[i] || 0)
            includedEdge.append(x0, 0); includedEdge.append(x0, low)
            includedEdge.append(x1, low); includedEdge.append(x1, 0)
            excludedBottom.append(x0, 0); excludedBottom.append(x0, low)
            excludedBottom.append(x1, low); excludedBottom.append(x1, 0)
            excludedTop.append(x0, 0); excludedTop.append(x0, high)
            excludedTop.append(x1, high); excludedTop.append(x1, 0)
        }
    }
    onIncludedBinsChanged: Qt.callLater(rebuild)
    onExcludedBinsChanged: Qt.callLater(rebuild)
    onMaximumChanged: Qt.callLater(rebuild)
    Component.onCompleted: rebuild()
    DiagnosticAxis {
        id: horizontal
        min: 0; max: chart.maximum
        titleText: qsTr("Residual magnitude, arcsec")
        labelFormat: chart.maximum < 1 ? "%.2f" : "%.1f"
    }
    DiagnosticAxis {
        id: vertical
        min: 0
        max: {
            var peak = 0
            for (var i = 0; i < chart.includedBins.length; ++i)
                peak = Math.max(peak, chart.includedBins[i] + (chart.excludedBins[i] || 0))
            return Math.max(1, Math.ceil(peak * 1.08))
        }
        labelFormat: "%.0f"
        titleText: qsTr("Observations")
    }
    AreaSeries {
        axisX: horizontal; axisY: vertical
        color: Qt.alpha(AppTheme.accent, 0.7); borderColor: "transparent"
        upperSeries: LineSeries { id: includedEdge; objectName: "histogram-included" }
    }
    AreaSeries {
        axisX: horizontal; axisY: vertical
        color: AppTheme.red; borderColor: "transparent"
        lowerSeries: LineSeries { id: excludedBottom }
        upperSeries: LineSeries { id: excludedTop; objectName: "histogram-excluded" }
    }
    MouseArea {
        anchors.fill: parent
        acceptedButtons: Qt.NoButton
        hoverEnabled: true
        onExited: chart.hoveredBin = -1
        onPositionChanged: function(mouse) {
            var point = chart.mapToValue(Qt.point(mouse.x, mouse.y))
            var bin = Math.floor(point.x / chart.binWidth)
            chart.hoveredBin = bin >= 0 && bin < chart.includedBins.length && point.y >= 0
                && point.y <= chart.includedBins[bin] + (chart.excludedBins[bin] || 0) ? bin : -1
        }
    }
    Rectangle {
        visible: chart.hoveredBin >= 0
        x: Math.max(8, chart.width - width - 24); y: 8
        width: hint.implicitWidth + 20; height: hint.implicitHeight + 16
        color: AppTheme.surface; border.color: AppTheme.controlLine; radius: 4
        CustomText {
            id: hint
            anchors.centerIn: parent
            font.pixelSize: 11
            text: chart.hoveredBin < 0 ? "" : (chart.hoveredBin * chart.binWidth).toFixed(3)
                + " - " + ((chart.hoveredBin + 1) * chart.binWidth).toFixed(3) + " arcsec\n"
                + qsTr("Included: %1 | Excluded: %2").arg(chart.includedBins[chart.hoveredBin]).arg(chart.excludedBins[chart.hoveredBin])
        }
    }
}
