pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

ColumnLayout {
    id: table
    required property var rows
    property string parameterTitle: qsTr("Parameter")
    property string rowPrefix: "parameter"
    property bool showSigma: false
    property bool showStatus: false
    property bool showSignificance: false
    property real columnSpacing: 12
    property real parameterShare: showSigma || showStatus ? 0.16 : 0.26
    property real unitShare: showSigma || showStatus ? 0.14 : 0.18
    property real sigmaShare: 0.20
    property real statusShare: 0.20
    property real significanceShare: 0.13
    readonly property real columnSpace: Math.max(0, width - 24 - columnSpacing * (2 + Number(showSigma) + Number(showStatus) + Number(showSignificance)))
    readonly property real parameterWidth: columnSpace * parameterShare
    readonly property real valueWidth: columnSpace - parameterWidth - sigmaWidth - unitWidth - statusWidth - significanceWidth
    readonly property real sigmaWidth: showSigma ? columnSpace * sigmaShare : 0
    readonly property real unitWidth: columnSpace * unitShare
    readonly property real statusWidth: showStatus ? columnSpace * statusShare : 0
    readonly property real significanceWidth: showSignificance ? columnSpace * significanceShare : 0
    spacing: 0

    Rectangle {
        Layout.fillWidth: true
        implicitHeight: 32
        color: AppTheme.well
        Row {
            anchors.fill: parent
            anchors.leftMargin: 12
            anchors.rightMargin: 12
            spacing: table.columnSpacing
            CustomText { text: table.parameterTitle; width: table.parameterWidth; height: parent.height; font.pixelSize: 10; secondary: true }
            CustomText { text: qsTr("Value"); width: table.valueWidth; height: parent.height; font.pixelSize: 10; secondary: true; horizontalAlignment: Text.AlignRight }
            CustomText { visible: table.showSigma; text: "\u03c3"; width: table.sigmaWidth; height: parent.height; font.pixelSize: 10; secondary: true; horizontalAlignment: Text.AlignRight }
            CustomText { visible: table.showSignificance; text: qsTr("|value| / σ"); width: table.significanceWidth; height: parent.height; font.pixelSize: 10; secondary: true; horizontalAlignment: Text.AlignRight }
            CustomText { text: qsTr("Unit"); width: table.unitWidth; height: parent.height; font.pixelSize: 10; secondary: true }
            CustomText { visible: table.showStatus; text: qsTr("Status"); width: table.statusWidth; height: parent.height; font.pixelSize: 10; secondary: true; horizontalAlignment: Text.AlignRight }
        }
        Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line }
    }
    Repeater {
        model: table.rows
        Item {
            id: row
            required property var modelData
            required property int index
            objectName: table.rowPrefix + "-row-" + modelData.label
            Layout.fillWidth: true
            implicitHeight: 36
            Rectangle { anchors.fill: parent; color: AppTheme.surfaceSecondary; visible: row.index % 2 === 1 }
            Row {
                anchors.fill: parent
                anchors.leftMargin: 12
                anchors.rightMargin: 12
                spacing: table.columnSpacing
                CustomText { text: row.modelData.label; width: table.parameterWidth; height: parent.height; font.pixelSize: 12 }
                CustomText {
                    objectName: table.rowPrefix + "-value"
                    text: row.modelData.value
                    width: table.valueWidth
                    height: parent.height
                    numeric: true
                    font.pixelSize: 12
                    font.weight: Font.DemiBold
                    horizontalAlignment: Text.AlignRight
                }
                CustomText {
                    objectName: table.rowPrefix + "-sigma"
                    visible: table.showSigma
                    text: table.showSigma ? row.modelData.sigma : ""
                    width: table.sigmaWidth
                    height: parent.height
                    numeric: true
                    secondary: true
                    font.pixelSize: 11
                    horizontalAlignment: Text.AlignRight
                }
                CustomText {
                    objectName: table.rowPrefix + "-significance"
                    visible: table.showSignificance
                    text: table.showSignificance ? row.modelData.significance : ""
                    width: table.significanceWidth
                    height: parent.height
                    numeric: true
                    secondary: true
                    font.pixelSize: 11
                    horizontalAlignment: Text.AlignRight
                }
                CustomText { text: row.modelData.unit; width: table.unitWidth; height: parent.height; secondary: true; font.pixelSize: 10 }
                CustomText {
                    objectName: table.rowPrefix + "-status"
                    visible: table.showStatus
                    text: table.showStatus ? row.modelData.status : ""
                    width: table.statusWidth
                    height: parent.height
                    font.pixelSize: 11
                    color: row.modelData.estimated ? AppTheme.accent : row.modelData.used ? AppTheme.ink : AppTheme.muted
                    horizontalAlignment: Text.AlignRight
                }
            }
            Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.lineSoft }
        }
    }
}
