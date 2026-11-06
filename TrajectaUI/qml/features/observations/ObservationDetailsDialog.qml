import QtQuick
import QtQuick.Controls.Basic
import TrajectaUI

BaseDialog {
    id: dialog
    objectName: "observation-details-dialog"
    required property ObservationTableModel observations
    title: observations.selectionTitle
    width: Math.min(550, parent.width - 40)
    height: Math.min(implicitHeight, parent.height - 40)
    contentItem: ScrollView {
        id: scroll
        contentWidth: availableWidth
        implicitHeight: details.implicitHeight
        clip: true
        ScrollBar.vertical: CustomScrollBar {
            parent: scroll
            x: scroll.width - width
            y: scroll.topPadding
            height: scroll.availableHeight
        }
        ObservationDetails {
            id: details
            width: parent.width
            observations: dialog.observations
        }
    }
}
