import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI

ColumnLayout {
    id: inspector
    required property ObservationTableModel observations
    signal closeRequested()
    spacing: 16
    RowLayout {
        Layout.fillWidth: true
        spacing: 12
        CustomText { text: inspector.observations.selectionTitle; font.pixelSize: 14; font.weight: Font.DemiBold; Layout.fillWidth: true }
        CustomButton { objectName: "observations-inspector-close"; iconName: "x"; flat: true; accessibleLabel: qsTr("Close observation details"); onClicked: inspector.closeRequested() }
    }
    ScrollView {
        id: scroll
        Layout.fillWidth: true
        Layout.fillHeight: true
        contentWidth: availableWidth
        clip: true
        ScrollBar.vertical: CustomScrollBar {
            parent: scroll
            x: scroll.width - width
            y: scroll.topPadding
            height: scroll.availableHeight
        }
        ObservationDetails { width: parent.width; observations: inspector.observations }
    }
}
