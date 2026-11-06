pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI

Item {
    id: root
    required property CloseApproachTableModel results
    property string emptyText
    readonly property var tooltips: [qsTr("TDB; plus/minus local linear 3-sigma time uncertainty"), qsTr("Encounter body"),
        qsTr("Nominal body-centred distance"), qsTr("Local linear 3-sigma lower bound, truncated at zero"),
        qsTr("Local linear 3-sigma upper bound"), qsTr("Relative speed at closest approach"), qsTr("Hyperbolic excess speed")]
    HorizontalHeaderView {
        id: header
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: parent.top
        syncView: table
        clip: true
        delegate: Rectangle {
            id: heading
            required property var model
            required property int column
            implicitHeight: 52
            color: AppTheme.base
            RowLayout {
                anchors.fill: parent
                anchors.margins: 6
                spacing: 3
                CustomText {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    text: heading.model.display
                    font.pixelSize: 10
                    color: AppTheme.muted
                    horizontalAlignment: heading.column >= 2 ? Text.AlignRight : Text.AlignLeft
                    wrapMode: Text.Wrap
                }
                CustomIcon {
                    size: 12
                    color: root.results.sortColumn === heading.column ? AppTheme.accent : AppTheme.muted
                    name: "chevron-right"
                    rotation: root.results.sortColumn === heading.column && root.results.descending ? 90 : -90
                    opacity: root.results.sortColumn === heading.column ? 1 : 0.35
                }
            }
            HoverHandler { id: hover }
            TapHandler { onTapped: root.results.toggleSort(heading.column) }
            ToolTip.visible: hover.hovered
            ToolTip.delay: 500
            ToolTip.text: root.tooltips[heading.column]
            Accessible.role: Accessible.Button
            Accessible.name: heading.model.display
            Accessible.onPressAction: root.results.toggleSort(heading.column)
            Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line }
        }
    }
    TableView {
        id: table
        objectName: "approaches-table"
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: header.bottom
        anchors.bottom: parent.bottom
        model: root.results
        clip: true
        boundsBehavior: Flickable.StopAtBounds
        columnWidthProvider: function(column) {
            const widths = [230, 64, 125, 125, 125, 105, 105]
            const total = widths.reduce((a, b) => a + b, 0)
            return widths[column] + Math.max(0, width - total - 14) / widths.length
        }
        rowHeightProvider: row => 34
        onWidthChanged: forceLayout()
        delegate: Rectangle {
            id: cell
            required property string display
            required property string toolTip
            required property int column
            implicitHeight: 34
            color: AppTheme.base
            CustomText {
                anchors.fill: parent
                anchors.leftMargin: 6
                anchors.rightMargin: 6
                text: cell.display
                numeric: cell.column !== 1
                font.pixelSize: 11
                horizontalAlignment: cell.column >= 2 ? Text.AlignRight : Text.AlignLeft
            }
            HoverHandler { id: cellHover }
            ToolTip.visible: cellHover.hovered && cell.toolTip.length > 0
            ToolTip.delay: 500
            ToolTip.text: cell.toolTip
            Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.lineSoft }
        }
        ScrollBar.horizontal: CustomScrollBar {}
        ScrollBar.vertical: CustomScrollBar {}
    }
    CustomText {
        anchors.centerIn: table
        visible: root.results.count === 0
        text: root.emptyText
        secondary: true
        font.pixelSize: 12
    }
}
