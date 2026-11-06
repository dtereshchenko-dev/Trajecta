import QtQuick
import TrajectaUI

CustomButton {
    id: cell
    required property var model
    required property int column
    required property ObservationTableModel observations
    objectName: "observations-header-" + column
    implicitHeight: 38
    leftPadding: 8
    rightPadding: 8
    flat: true
    accessibleLabel: qsTr("Sort by %1").arg(model.display)
    contentItem: Item {
        CustomText {
            anchors.fill: parent
            anchors.leftMargin: cell.column > 2 ? 12 : 0
            anchors.rightMargin: cell.column > 2 ? 0 : 12
            text: cell.model.display
            color: cell.observations.sortedColumn === cell.column ? AppTheme.ink : AppTheme.muted
            font.pixelSize: 10
            horizontalAlignment: cell.column > 2 ? Text.AlignRight : Text.AlignLeft
        }
        CustomIcon {
            x: cell.column > 2 ? 0 : parent.width - width
            anchors.verticalCenter: parent.verticalCenter
            visible: cell.observations.sortedColumn === cell.column
            name: "chevron-right"
            rotation: cell.observations.sortDescending ? 90 : -90
            size: 10
            color: AppTheme.accent
        }
    }
    background: Rectangle {
        color: cell.hovered ? AppTheme.rowHover : AppTheme.base
        Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line }
    }
    onClicked: cell.observations.sortColumn(column, cell.observations.sortedColumn === column ? !cell.observations.sortDescending : false)
}
