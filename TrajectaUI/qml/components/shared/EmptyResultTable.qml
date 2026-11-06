pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import TrajectaUI

Item {
    id: root
    property var columns: []
    property var columnTooltips: []
    property var columnWidths: []
    readonly property real minimumContentWidth: columns.reduce((sum, column, index) => sum + (columnWidths[index] || 150), 0)
    property string emptyText: qsTr("No results")
    property real tableHeight: 240
    property real headerHeight: 52
    implicitHeight: tableHeight
    ScrollView {
        id: scroll
        anchors.fill: parent
        contentWidth: Math.max(availableWidth, root.minimumContentWidth)
        contentHeight: availableHeight
        clip: true
        ScrollBar.horizontal: CustomScrollBar {}
        ScrollBar.vertical.policy: ScrollBar.AlwaysOff
        Item {
            width: scroll.contentWidth
            height: scroll.availableHeight
            Row {
                id: headings
                width: parent.width
                Repeater {
                    model: root.columns
                    Item {
                        id: heading
                        required property string modelData
                        required property int index
                        width: (root.columnWidths[index] || 150) + (headings.width - root.minimumContentWidth) / root.columns.length
                        height: root.headerHeight
                        HoverHandler { id: hover }
                        ToolTip.visible: hover.hovered && ToolTip.text.length > 0
                        ToolTip.delay: 500
                        ToolTip.text: root.columnTooltips[index] || ""
                        CustomText {
                            anchors.fill: parent
                            anchors.margins: 10
                            text: heading.modelData
                            secondary: true
                            font.pixelSize: 10
                            wrapMode: Text.Wrap
                            verticalAlignment: Text.AlignVCenter
                        }
                    }
                }
            }
            Rectangle { y: root.headerHeight - 1; width: parent.width; height: 1; color: AppTheme.line }
        }
    }
    CustomText {
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.verticalCenter: parent.verticalCenter
        anchors.margins: 20
        text: root.emptyText
        secondary: true
        font.pixelSize: 12
        horizontalAlignment: Text.AlignHCenter
        wrapMode: Text.Wrap
    }
}
