pragma ComponentBehavior: Bound
import QtQuick
import TrajectaUI

Column {
    id: list
    property var rows: []
    spacing: 14
    Repeater {
        model: list.rows
        Item {
            id: row
            required property var modelData
            width: list.width
            readonly property real textBaseline: Math.max(label.baselineOffset, value.baselineOffset)
            implicitHeight: Math.max(label.y + label.implicitHeight, value.y + value.implicitHeight)
            CustomText {
                id: label
                y: row.textBaseline - baselineOffset
                width: Math.max(0, row.width - value.width - 18)
                text: row.modelData[0]
                secondary: true
                font.pixelSize: 12
                wrapMode: Text.Wrap
            }
            CustomText {
                id: value
                y: row.textBaseline - baselineOffset
                anchors.right: parent.right
                width: Math.min(implicitWidth, Math.max(0, row.width * 0.6))
                text: row.modelData[1]
                numeric: true
                font.pixelSize: 11
                font.weight: Font.DemiBold
                horizontalAlignment: Text.AlignRight
                wrapMode: Text.WrapAnywhere
                elide: Text.ElideNone
            }
        }
    }
}
