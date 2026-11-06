pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import TrajectaUI

Flickable {
    id: matrixView
    objectName: "covariance-matrix"
    required property CovarianceTableModel matrix
    readonly property real rowLabelWidth: 100
    readonly property real cellWidth: Math.max(84, (contentWidth - rowLabelWidth - matrix.dimension * 5) / Math.max(1, matrix.dimension))
    readonly property real cellHeight: 68
    implicitHeight: contentHeight + (contentWidth > width ? 14 : 0)
    contentWidth: Math.max(width, rowLabelWidth + matrix.dimension * 89)
    contentHeight: 36 + matrix.dimension * (cellHeight + 5)
    flickableDirection: Flickable.HorizontalFlick
    boundsBehavior: Flickable.StopAtBounds
    clip: true
    ScrollBar.horizontal: CustomScrollBar {}
    Repeater {
        model: matrixView.matrix.labels
        CustomText {
            required property int index
            required property string modelData
            x: matrixView.rowLabelWidth + 5 + index * (matrixView.cellWidth + 5)
            width: matrixView.cellWidth
            height: 31
            text: modelData
            secondary: true
            numeric: true
            font.pixelSize: 12
            horizontalAlignment: Text.AlignHCenter
        }
    }
    Repeater {
        model: matrixView.matrix.labels
        CustomText {
            required property int index
            required property string modelData
            y: 36 + index * (matrixView.cellHeight + 5)
            width: matrixView.rowLabelWidth
            height: matrixView.cellHeight
            text: modelData
            secondary: true
            numeric: true
            font.pixelSize: 12
            horizontalAlignment: Text.AlignHCenter
        }
    }
    TableView {
        id: table
        x: matrixView.rowLabelWidth + 5
        y: 36
        width: matrixView.contentWidth - x
        height: matrixView.contentHeight - y
        model: matrixView.matrix
        interactive: false
        columnSpacing: 5
        rowSpacing: 5
        columnWidthProvider: function(column) { return matrixView.cellWidth }
        rowHeightProvider: function(row) { return matrixView.cellHeight }
        delegate: CustomMatrixCell {
            required property int row
            required property int column
            required property var model
            objectName: "covariance-cell-" + row + "-" + column
            text: model.display === undefined || model.display === null ? "\u2014" : model.display
            coefficient: model.rho
            checked: matrixView.matrix.selectedRow === row && matrixView.matrix.selectedColumn === column
            accessibleLabel: matrixView.matrix.labels[row] + " / " + matrixView.matrix.labels[column] + ": " + text
            onClicked: matrixView.matrix.selectCell(row, column)
        }
        Connections {
            target: matrixView
            function onCellWidthChanged() { table.forceLayout() }
        }
    }
}
