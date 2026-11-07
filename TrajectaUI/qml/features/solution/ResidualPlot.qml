import QtQuick
import QtCharts
import TrajectaUI
import "../../shared/UI.js" as UI

ChartView {
    id: chart
    objectName: "residual-plot"
    property ObservationTableModel observations: null
    property bool timeMode: false
    property int timeComponent: -1
    property real timeOriginMjd: 0
    property string timeOriginLabel: ""
    property bool externalData: false
    property bool uniformColor: false
    property string yLabelFormat: "%.2f"
    property alias raSeries: included
    property alias decSeries: includedDec
    property alias excludedRaSeries: excluded
    property alias excludedDecSeries: excludedDec
    property real referenceRms: 0
    property rect bounds: Qt.rect(-1, -1, 2, 2)
    property rect dataBounds: Qt.rect(-1, -1, 2, 2)
    property var hoveredSeries: null
    antialiasing: true
    backgroundColor: AppTheme.surface
    plotAreaColor: AppTheme.surface
    plotArea: Qt.rect(66, 16, Math.max(1, width - 90), Math.max(1, height - 64))
    backgroundRoundness: 0
    legend.visible: false
    margins { left: 0; right: 0; top: 0; bottom: 0 }
    clip: true
    MouseArea {
        anchors.fill: parent
        hoverEnabled: true
        acceptedButtons: Qt.LeftButton
        cursorShape: chart.hoveredSeries ? Qt.PointingHandCursor : Qt.ArrowCursor
        onPositionChanged: mouse => chart.pointAt(mouse.x, mouse.y, false)
        onClicked: mouse => chart.pointAt(mouse.x, mouse.y, true)
        onExited: chart.clearHover()
        onCanceled: chart.clearHover()
        onWheel: wheel => wheel.accepted = false
    }
    function pointAt(x, y, select) {
        clearHover()
        var entries = [
            [selectedDec, true, true, true], [selected, true, false, true],
            [excludedDec, false, true, false], [includedDec, true, true, false],
            [excluded, false, false, false], [included, true, false, false]
        ]
        for (var entry of entries) {
            var series = entry[0]
            if (!series.visible || series.count === 0) continue
            var index = included.markerRenderer.pointIndexAt(series, Qt.point(x, y))
            if (index < 0) continue
            var point = series.at(index)
            if (series.pointHovered !== undefined) series.pointHovered = true
            updateHover(series, point, true, entry[1], entry[2], entry[3])
            if (select && !entry[3]) selectPoint(point, entry[1], entry[2])
            return
        }
    }
    function refreshSelection() {
        if (!observations) { selected.clear(); selectedDec.clear(); return }
        observations.fillSeries(selected, true, timeMode, false, true, timeOriginMjd)
        observations.fillSeries(selectedDec, true, timeMode, true, true, timeOriginMjd)
    }
    function clearHover() {
        hoveredSeries = null
        pointHint.observation = ({})
        included.pointHovered = false
        excluded.pointHovered = false
        includedDec.pointHovered = false
        excludedDec.pointHovered = false
    }
    function updateHover(series, point, state, isIncluded, declination, selectedOnly) {
        if (!observations) return
        if (!state) {
            if (hoveredSeries === series) clearHover()
            return
        }
        hoveredSeries = series
        pointHint.observation = selectedOnly ? observations.selection
            : observations.pointData(point.x + (timeMode ? timeOriginMjd : 0), point.y, timeMode, declination, isIncluded)
        pointHint.anchorPoint = mapToPosition(point, series)
    }
    function selectPoint(point, isIncluded, declination) {
        if (!observations) return
        var observation = observations.pointData(point.x + (timeMode ? timeOriginMjd : 0), point.y, timeMode, declination, isIncluded)
        if (observation.sourceRow !== undefined) observations.selectSourceRow(observation.sourceRow)
    }
    ResidualPointHint { id: pointHint; objectName: "residual-point-hint" }
    function updateBounds() {
        if (externalData) return
        var raw = dataBounds
        if (timeMode || plotArea.width <= 0 || plotArea.height <= 0) {
            bounds = raw
            return
        }
        var reference = Number.isFinite(referenceRms) && referenceRms > 0 ? referenceRms * 3.5 : 0
        var unitsPerPixel = Math.max(Math.max(raw.width / 2, reference) / plotArea.width,
                                     Math.max(raw.height / 2, reference) / plotArea.height)
        var x = unitsPerPixel * plotArea.width
        var y = unitsPerPixel * plotArea.height
        if (Math.abs(bounds.width - 2 * x) < 1e-9 && Math.abs(bounds.height - 2 * y) < 1e-9) return
        bounds = Qt.rect(-x, -y, 2 * x, 2 * y)
        if (zeroHorizontal && zeroHorizontal.count === 2) {
            zeroHorizontal.replace(0, -x, 0)
            zeroHorizontal.replace(1, x, 0)
        }
        if (zeroVertical && zeroVertical.count === 2) {
            zeroVertical.replace(0, 0, -y)
            zeroVertical.replace(1, 0, y)
        }
    }
    function refresh() {
        clearHover()
        if (externalData || !observations) { refreshSelection(); return }
        dataBounds = observations.plotBounds(timeMode)
        updateBounds()
        zeroHorizontal.clear()
        zeroHorizontal.append(bounds.x, 0)
        zeroHorizontal.append(bounds.x + bounds.width, 0)
        zeroVertical.clear()
        zeroVertical.append(0, bounds.y)
        zeroVertical.append(0, bounds.y + bounds.height)
        observations.fillSeries(included, true, timeMode, false)
        observations.fillSeries(excluded, false, timeMode, false)
        observations.fillSeries(includedDec, true, timeMode, true)
        observations.fillSeries(excludedDec, false, timeMode, true)
        refreshSelection()
    }
    function refreshReferences() {
        if (!zeroHorizontal) return
        zeroHorizontal.clear()
        zeroHorizontal.append(bounds.x, 0)
        zeroHorizontal.append(bounds.x + bounds.width, 0)
    }
    ValueAxis {
        id: horizontal
        min: chart.bounds.x
        max: chart.bounds.x + chart.bounds.width
        tickCount: 5
        titleText: chart.timeMode ? chart.timeOriginLabel.length ? qsTr("Days from %1").arg(chart.timeOriginLabel) : "MJD UTC" : "Δα cos δ, arcsec"
        labelsColor: AppTheme.chartTick
        titleBrush: AppTheme.muted
        gridLineColor: AppTheme.chartGrid
        color: AppTheme.chartGrid
        labelsFont.family: Qt.platform.os === "windows" ? "Consolas" : UI.monoFont
        labelsFont.pixelSize: 10
        titleFont.family: Qt.platform.os === "windows" ? "Segoe UI" : UI.font
        titleFont.pixelSize: 10
        labelFormat: chart.timeMode ? "%.1f" : "%.2f"
    }
    ValueAxis {
        id: vertical
        min: chart.bounds.y
        max: chart.bounds.y + chart.bounds.height
        tickCount: 5
        titleText: chart.timeMode ? chart.timeComponent === 0 ? "Δα cos δ, arcsec" : chart.timeComponent === 1 ? "Δδ, arcsec" : "Residual, arcsec" : "Δδ, arcsec"
        labelsColor: AppTheme.chartTick
        titleBrush: AppTheme.muted
        gridLineColor: AppTheme.chartGrid
        color: AppTheme.chartGrid
        labelsFont.family: Qt.platform.os === "windows" ? "Consolas" : UI.monoFont
        labelsFont.pixelSize: 10
        titleFont.family: Qt.platform.os === "windows" ? "Segoe UI" : UI.font
        titleFont.pixelSize: 10
        labelFormat: chart.yLabelFormat
    }
    LineSeries {
        id: zeroHorizontal
        axisX: horizontal; axisY: vertical; color: AppTheme.chartZero
    }
    LineSeries {
        id: zeroVertical
        visible: !chart.timeMode
        axisX: horizontal; axisY: vertical; color: AppTheme.chartZero
    }
    ResidualReferenceCircle {
        objectName: "residual-rms-circle"
        axisX: horizontal; axisY: vertical
        radius: chart.referenceRms
        visible: !chart.timeMode && radius > 0
        color: Qt.alpha(AppTheme.blue, 0.5)
    }
    ResidualReferenceCircle {
        objectName: "residual-3rms-circle"
        axisX: horizontal; axisY: vertical
        radius: chart.referenceRms * 3
        visible: !chart.timeMode && radius > 0
        color: Qt.alpha(AppTheme.amber, 0.6)
        style: Qt.DashLine
    }
    LineSeries {
        id: included
        property ResidualMarker markerRenderer: ResidualMarker { series: included }
        style: Qt.NoPen
        visible: !chart.timeMode || chart.timeComponent !== 1
        objectName: "included-residuals"
        property bool pointHovered: false
        axisX: horizontal; axisY: vertical
        color: Qt.alpha(AppTheme.accent, 0.72)
    }
    LineSeries {
        id: excluded
        property ResidualMarker markerRenderer: ResidualMarker { series: excluded }
        style: Qt.NoPen
        visible: !chart.timeMode || chart.timeComponent !== 1
        property bool pointHovered: false
        axisX: horizontal; axisY: vertical
        color: AppTheme.red
    }
    LineSeries {
        id: includedDec
        property ResidualMarker markerRenderer: ResidualMarker { series: includedDec }
        style: Qt.NoPen
        objectName: "included-dec-residuals"
        property bool pointHovered: false
        visible: chart.timeMode && chart.timeComponent !== 0
        axisX: horizontal; axisY: vertical
        color: Qt.alpha(chart.uniformColor ? AppTheme.accent : AppTheme.blue, 0.72)
    }
    LineSeries {
        id: excludedDec
        property ResidualMarker markerRenderer: ResidualMarker { series: excludedDec }
        style: Qt.NoPen
        property bool pointHovered: false
        visible: chart.timeMode && chart.timeComponent !== 0
        axisX: horizontal; axisY: vertical
        color: AppTheme.red
    }
    ScatterSeries {
        id: selected
        visible: !chart.timeMode || chart.timeComponent !== 1
        axisX: horizontal; axisY: vertical; color: "transparent"; borderColor: AppTheme.accent; markerSize: 15
    }
    ScatterSeries {
        id: selectedDec
        visible: chart.timeMode && chart.timeComponent !== 0
        axisX: horizontal; axisY: vertical; color: "transparent"; borderColor: chart.uniformColor ? AppTheme.accent : AppTheme.blue; markerSize: 15
    }
    Component.onCompleted: refresh()
    onTimeModeChanged: refresh()
    onTimeOriginMjdChanged: refreshSelection()
    onBoundsChanged: refreshReferences()
    onReferenceRmsChanged: updateBounds()
    onPlotAreaChanged: { clearHover(); updateBounds() }
    Connections {
        target: chart.observations
        function onCountChanged() { chart.refresh() }
        function onSelectionChanged() { chart.refreshSelection() }
    }
}
