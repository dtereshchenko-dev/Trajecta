import QtQuick
import TrajectaUI

Rectangle {
    id: cell
    required property var display
    required property int sourceRow
    required property bool included
    required property bool desiredIncluded
    required property bool pending
    required property bool manual
    required property int row
    required property int column
    required property ObservationTableModel observations
    required property int hoveredSourceRow
    signal hoverChanged(int source, bool hovered)
    signal activated()
    objectName: "observation-cell-" + sourceRow + "-" + column
    implicitWidth: 120
    implicitHeight: 36
    color: observations.selectedSourceRow === sourceRow ? AppTheme.accentSoft : hoveredSourceRow === sourceRow ? AppTheme.rowHover : AppTheme.base
    CustomText {
        anchors.fill: parent
        anchors.margins: 8
        anchors.leftMargin: cell.column === 0 ? 32 : 8
        text: cell.display === undefined || cell.display === null ? "\u2014" : cell.display
        numeric: cell.column < 10
        font.pixelSize: 11
        color: cell.column === 10 && cell.pending ? AppTheme.amber : cell.column === 2 ? AppTheme.accent : cell.column === 10 && !cell.included ? AppTheme.red : cell.included ? AppTheme.ink : AppTheme.muted
        horizontalAlignment: cell.column > 2 ? Text.AlignRight : Text.AlignLeft
    }
    CustomCheckBox {
        objectName: cell.column === 0 ? "observation-included-" + cell.sourceRow : ""
        anchors.left: parent.left
        anchors.leftMargin: 6
        anchors.verticalCenter: parent.verticalCenter
        width: 24
        height: 28
        visible: cell.column === 0
        enabled: cell.observations.editingEnabled
        opacity: enabled ? 1 : 0.5
        checked: cell.desiredIncluded
        Accessible.name: qsTr("Include observation %1 in fit").arg(cell.sourceRow + 1)
        onToggled: cell.observations.setIncluded(cell.sourceRow, checked)
    }
    Rectangle {
        width: 4; height: 4; radius: 2
        anchors.left: parent.left; anchors.leftMargin: 24
        anchors.top: parent.top; anchors.topMargin: 4
        visible: cell.column === 0 && (cell.pending || cell.manual)
        color: cell.pending ? AppTheme.amber : AppTheme.accent
    }
    Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.lineSoft }
    Rectangle { width: 2; height: parent.height; color: AppTheme.accent; visible: cell.column === 0 && cell.observations.selectedSourceRow === cell.sourceRow }
    InteractionCursor {}
    HoverHandler { onHoveredChanged: cell.hoverChanged(cell.sourceRow, hovered) }
    TapHandler { enabled: cell.column !== 0; onTapped: { cell.observations.selectRow(cell.row); cell.activated() } }
    Accessible.role: Accessible.Cell
    Accessible.name: qsTr("Observation %1, %2").arg(sourceRow + 1).arg(display === undefined || display === null ? qsTr("Not available") : display)
}
