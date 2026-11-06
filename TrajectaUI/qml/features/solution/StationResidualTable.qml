pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import TrajectaUI

Flickable {
    id: table
    required property StationResidualModel statistics
    signal stationRequested(string station)
    implicitHeight: content.implicitHeight + (contentWidth > width ? 12 : 0)
    contentWidth: Math.max(760, width)
    contentHeight: content.implicitHeight
    flickableDirection: Flickable.HorizontalFlick
    boundsBehavior: Flickable.StopAtBounds
    clip: true
    ScrollBar.horizontal: CustomScrollBar {}
    Column {
        id: content
        width: table.contentWidth
        StationResidualRow { width: parent.width; header: true; meanLimit: table.statistics.meanLimit }
        Repeater {
            model: table.statistics
            StationResidualRow {
                required property var model
                objectName: "station-residual-row-" + model.station
                width: content.width
                station: model.station
                includedCount: model.included
                excludedCount: model.excluded
                rms: model.rms
                meanRa: model.meanRa
                meanDec: model.meanDec
                meanLimit: table.statistics.meanLimit
                onStationRequested: function(station) { table.stationRequested(station) }
            }
        }
    }
}
