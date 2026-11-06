pragma ComponentBehavior: Bound
import QtQuick
import TrajectaUI
import "../../shared/UI.js" as UI

Rectangle {
    id: row
    property bool header: false
    property string station: ""
    property int includedCount: 0
    property int excludedCount: 0
    property real rms: NaN
    property real meanRa: NaN
    property real meanDec: NaN
    property real meanLimit: 1
    signal stationRequested(string station)
    implicitHeight: header ? 60 : 38
    color: header ? AppTheme.well : hover.hovered ? AppTheme.rowHover : "transparent"
    HoverHandler { id: hover; enabled: !row.header }
    Row {
        anchors.fill: parent
        Repeater {
            model: 3
            Item {
                required property int index
                width: row.width * [0.10, 0.08, 0.08][index]
                height: row.height
                CustomText {
                    visible: row.header || parent.index > 0
                    anchors.fill: parent
                    anchors.leftMargin: 14
                    anchors.rightMargin: 14
                    verticalAlignment: Text.AlignVCenter
                    horizontalAlignment: parent.index === 0 ? Text.AlignLeft : Text.AlignRight
                    text: row.header ? [qsTr("Station"), qsTr("Included"), qsTr("Excluded")][parent.index]
                        : [row.station, String(row.includedCount), String(row.excludedCount)][parent.index]
                    numeric: !row.header
                    font.pixelSize: row.header ? 10 : 11
                    color: row.header ? AppTheme.muted : parent.index === 2 ? row.excludedCount ? AppTheme.red : AppTheme.muted : AppTheme.ink
                }
                CustomTextLink {
                    objectName: "station-link-" + row.station
                    visible: !row.header && parent.index === 0
                    anchors.left: parent.left
                    anchors.leftMargin: 14
                    anchors.verticalCenter: parent.verticalCenter
                    text: row.station || qsTr("Unknown")
                    iconName: "arrow-up-right"
                    Accessible.name: qsTr("Observations from station %1").arg(row.station)
                    onClicked: row.stationRequested(row.station)
                }
            }
        }
        Repeater {
            model: 2
            Item {
                required property int index
                width: row.width * 0.31
                height: row.height
                ResidualOffset {
                    anchors.fill: parent
                    anchors.leftMargin: 14
                    anchors.rightMargin: 14
                    header: row.header
                    label: parent.index === 0 ? qsTr("Mean \u0394\u03b1 cos \u03b4") : qsTr("Mean \u0394\u03b4")
                    value: parent.index === 0 ? row.meanRa : row.meanDec
                    limit: row.meanLimit
                    tint: parent.index === 0 ? AppTheme.accent : AppTheme.blue
                }
            }
        }
        Item {
            width: row.width * 0.12
            height: row.height
            CustomText {
                objectName: row.header ? "station-rms-header" : "station-rms-value"
                anchors.fill: parent
                anchors.leftMargin: 14
                anchors.rightMargin: 14
                verticalAlignment: Text.AlignVCenter
                horizontalAlignment: Text.AlignRight
                text: row.header ? qsTr("RMS") : UI.number(row.rms, 3)
                numeric: !row.header
                font.pixelSize: row.header ? 10 : 11
                font.weight: row.header ? Font.Normal : Font.DemiBold
                color: row.header ? AppTheme.muted : AppTheme.ink
            }
        }
    }
    Rectangle { width: parent.width; height: 1; color: AppTheme.line; visible: row.header }
    Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: row.header ? AppTheme.line : AppTheme.lineSoft }
}
