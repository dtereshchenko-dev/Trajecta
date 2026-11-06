pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import TrajectaUI

FocusScope {
    id: root
    required property ObservationTableModel observations
    property int hoveredSourceRow: -1
    signal observationActivated()
    function moveSelection(delta) {
        if (!observations.count) return
        const row = observations.selectedProxyRow < 0 ? 0 : Math.max(0, Math.min(observations.count - 1, observations.selectedProxyRow + delta))
        observations.selectRow(row)
        table.positionViewAtRow(row, TableView.Contain)
    }
    Keys.onUpPressed: moveSelection(-1)
    Keys.onDownPressed: moveSelection(1)
    Keys.onPressed: event => {
        if (observations.count && (event.key === Qt.Key_Home || event.key === Qt.Key_End)) {
            const row = event.key === Qt.Key_Home ? 0 : observations.count - 1
            observations.selectRow(row)
            table.positionViewAtRow(row, TableView.Contain)
            event.accepted = true
        }
    }
    Keys.onReturnPressed: { if (observations.selectedSourceRow >= 0) observationActivated() }
    Keys.onEnterPressed: { if (observations.selectedSourceRow >= 0) observationActivated() }
    HorizontalHeaderView {
        id: header
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: parent.top
        anchors.rightMargin: verticalBar.visible ? verticalBar.width : 0
        syncView: table
        clip: true
        delegate: ObservationHeaderCell { observations: root.observations }
    }
    TableView {
        id: table
        objectName: "observations-data-view"
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: header.bottom
        anchors.bottom: parent.bottom
        clip: true
        model: root.observations
        pointerNavigationEnabled: false
        keyNavigationEnabled: false
        boundsBehavior: Flickable.StopAtBounds
        columnWidthProvider: function(column) {
            const widths = root.observations.radarMode ? [72, 182, 70, 112, 112, 106, 118, 106, 112, 86, 118]
                : [72, 182, 70, 100, 112, 106, 86, 112, 96, 86, 118]
            const total = widths.reduce((a, b) => a + b, 0)
            return widths[column] + Math.max(0, width - total - 14) / widths.length
        }
        onWidthChanged: forceLayout()
        rowHeightProvider: function(row) { return 36 }
        delegate: ObservationCell {
            observations: root.observations
            hoveredSourceRow: root.hoveredSourceRow
            onHoverChanged: (source, hovered) => { if (hovered) root.hoveredSourceRow = source; else if (root.hoveredSourceRow === source) root.hoveredSourceRow = -1 }
            onActivated: { root.forceActiveFocus(); root.observationActivated() }
        }
        ScrollBar.horizontal: CustomScrollBar { id: horizontalBar; policy: ScrollBar.AsNeeded }
        ScrollBar.vertical: CustomScrollBar { id: verticalBar; policy: ScrollBar.AsNeeded }
        bottomMargin: horizontalBar.visible ? horizontalBar.height : 0
        rightMargin: verticalBar.visible ? verticalBar.width : 0
    }
    Connections {
        target: root.observations
        function onFilterChanged() { table.contentY = 0; root.hoveredSourceRow = -1; table.forceLayout() }
        function onSortingChanged() { table.contentY = 0; root.hoveredSourceRow = -1 }
    }
}
