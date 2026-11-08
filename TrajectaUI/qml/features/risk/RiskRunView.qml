pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI

ColumnLayout {
    id: root
    required property RiskMonitorViewModel monitor
    property real tableHeight: 500
    spacing: 16
    RowLayout {
        Layout.fillWidth: true
        CustomText { text: qsTr("Events"); font.pixelSize: 14; font.weight: Font.DemiBold }
        CustomComboBox {
            objectName: "risk-run-history"
            Layout.fillWidth: true
            model: root.monitor.runs
            currentIndex: root.monitor.selectedRun
            enabled: count > 0
            onActivated: index => root.monitor.selectedRun = index
            Accessible.name: qsTr("MC run")
        }
        CustomButton {
            iconName: "file-text"
            enabled: root.monitor.archivePath.length > 0
            accessibleLabel: qsTr("Open HTML report")
            onClicked: root.monitor.openArchive(true)
        }
        CustomButton {
            iconName: "folder-open"
            enabled: root.monitor.archivePath.length > 0
            accessibleLabel: qsTr("Open numerical archive")
            onClicked: root.monitor.openArchive(false)
        }
    }
    CustomText {
        visible: root.monitor.busy
        text: root.monitor.progress
        Layout.fillWidth: true
        font.pixelSize: 11
    }
    CustomProgressBar {
        objectName: "risk-progress"
        visible: root.monitor.busy
        Layout.fillWidth: true
        value: root.monitor.progressFraction
        indeterminate: value === 0
    }
    CustomText {
        text: root.monitor.runDetails
        visible: text.length > 0
        Layout.fillWidth: true
        wrapMode: Text.Wrap
        elide: Text.ElideNone
        font.pixelSize: 11
        secondary: true
    }
    CustomText {
        text: root.monitor.summary
        visible: text.length > 0
        Layout.fillWidth: true
        wrapMode: Text.Wrap
        elide: Text.ElideNone
        font.pixelSize: 11
        secondary: true
    }
    CustomText {
        text: root.monitor.warning
        visible: text.length > 0
        Layout.fillWidth: true
        wrapMode: Text.Wrap
        elide: Text.ElideNone
        font.pixelSize: 11
        color: AppTheme.amber
    }
    CustomText {
        objectName: "risk-job-error"
        text: root.monitor.error
        visible: text.length > 0
        Layout.fillWidth: true
        wrapMode: Text.Wrap
        elide: Text.ElideNone
        font.pixelSize: 11
        color: AppTheme.red
    }
    Item {
        Layout.fillWidth: true
        Layout.preferredHeight: root.tableHeight
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
                implicitHeight: 38
                color: AppTheme.base
                CustomText {
                    anchors.fill: parent
                    anchors.margins: 8
                    text: heading.model.display
                    font.pixelSize: 10
                    color: AppTheme.muted
                    horizontalAlignment: heading.column ? Text.AlignRight : Text.AlignLeft
                }
                Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line }
            }
        }
        TableView {
            id: table
            objectName: "risk-event-table"
            anchors.fill: parent
            anchors.topMargin: header.height
            clip: true
            model: root.monitor.events
            boundsBehavior: Flickable.StopAtBounds
            columnWidthProvider: function(column) {
                const widths = [120, 216, 216, 80, 100, 104, 104]
                return widths[column] + Math.max(0, width - 940) / widths.length
            }
            onWidthChanged: forceLayout()
            rowHeightProvider: function(row) { return 36 }
            delegate: Rectangle {
                id: cell
                required property string display
                required property string detail
                required property int column
                implicitHeight: 36
                color: AppTheme.base
                CustomText {
                    anchors.fill: parent
                    anchors.leftMargin: 8
                    anchors.rightMargin: 8
                    text: cell.display
                    numeric: true
                    font.pixelSize: 11
                    horizontalAlignment: cell.column ? Text.AlignRight : Text.AlignLeft
                }
                HoverHandler { id: hover }
                ToolTip.visible: hover.hovered
                ToolTip.text: cell.detail
                ToolTip.delay: 450
                Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.lineSoft }
            }
            ScrollBar.horizontal: CustomScrollBar {}
            ScrollBar.vertical: CustomScrollBar {}
        }
        CustomText {
            anchors.centerIn: parent
            width: Math.min(parent.width, 500)
            horizontalAlignment: Text.AlignHCenter
            wrapMode: Text.Wrap
            visible: root.monitor.events.count === 0
            text: root.monitor.busy ? qsTr("Computing trajectories...")
                : root.monitor.selectedRun < 0 ? qsTr("No MC runs in this session")
                : qsTr("No boundary crossings in completed trials")
            secondary: true
        }
    }
    CustomText {
        visible: root.monitor.selectedRun >= 0
        text: qsTr("Crossing groups are separated by gaps over 24 hours. P = entries / complete trials; Wilson 95%. Zero entries is not zero risk.")
        Layout.fillWidth: true
        wrapMode: Text.Wrap
        elide: Text.ElideNone
        font.pixelSize: 10
        secondary: true
    }
}
