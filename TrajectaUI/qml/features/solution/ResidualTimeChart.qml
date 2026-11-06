import QtQuick
import TrajectaUI

ResidualPlot {
    id: chart
    property bool declination: false
    readonly property var includedSeries: declination ? decSeries : raSeries
    readonly property var excludedSeries: declination ? excludedDecSeries : excludedRaSeries
    timeMode: true
    timeComponent: declination ? 1 : 0
    externalData: true
    uniformColor: true
    yLabelFormat: "%.1f"
    plotArea: Qt.rect(0, 0, 0, 0)
}
