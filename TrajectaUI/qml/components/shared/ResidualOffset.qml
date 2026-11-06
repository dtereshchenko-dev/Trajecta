import QtQuick
import TrajectaUI
import "../../shared/UI.js" as UI

Item {
    id: offset
    property bool header: false
    property string label: ""
    property real value: NaN
    property real limit: 1
    property color tint: AppTheme.accent
    readonly property real plotWidth: Math.max(0, width - 76)
    readonly property real position: Number.isFinite(value) && limit > 0 ? plotWidth * (0.5 + Math.max(-1, Math.min(1, value / limit)) * 0.5) : plotWidth / 2
    implicitHeight: header ? 60 : 38
    Column {
        visible: offset.header
        anchors.verticalCenter: parent.verticalCenter
        width: parent.width
        spacing: 5
        Row {
            spacing: 6
            Rectangle { width: 6; height: 6; radius: 3; color: offset.tint; anchors.verticalCenter: parent.verticalCenter }
            CustomText { text: offset.label; secondary: true; font.pixelSize: 10 }
        }
        Item {
            width: offset.plotWidth
            height: 12
            CustomText { text: "-" + offset.limit; secondary: true; numeric: true; font.pixelSize: 9 }
            CustomText { text: "0"; anchors.horizontalCenter: parent.horizontalCenter; secondary: true; numeric: true; font.pixelSize: 9 }
            CustomText { text: "+" + offset.limit; anchors.right: parent.right; secondary: true; numeric: true; font.pixelSize: 9 }
        }
    }
    Item {
        visible: !offset.header
        anchors.fill: parent
        Item {
            visible: Number.isFinite(offset.value)
            width: offset.plotWidth
            height: parent.height
            Rectangle { width: parent.width; height: 1; anchors.verticalCenter: parent.verticalCenter; color: AppTheme.lineSoft }
            Rectangle { x: parent.width / 2; width: 1; height: parent.height; color: AppTheme.controlLine }
            Rectangle {
                x: Math.min(offset.plotWidth / 2, offset.position)
                width: Math.abs(offset.position - offset.plotWidth / 2)
                height: 6
                anchors.verticalCenter: parent.verticalCenter
                color: offset.tint
                opacity: 0.42
            }
            Rectangle {
                objectName: "residual-offset-dot"
                x: offset.position - width / 2
                width: 6
                height: 6
                radius: 3
                anchors.verticalCenter: parent.verticalCenter
                color: offset.tint
            }
        }
        CustomText {
            objectName: "residual-offset-value"
            width: 64
            anchors.right: parent.right
            anchors.verticalCenter: parent.verticalCenter
            text: UI.signedNumber(offset.value, 3)
            numeric: true
            font.pixelSize: 11
            horizontalAlignment: Text.AlignRight
            color: Number.isFinite(offset.value) ? AppTheme.ink : AppTheme.muted
        }
    }
}
