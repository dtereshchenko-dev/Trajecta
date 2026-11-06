import QtQuick
import QtQuick.Layouts
import TrajectaUI

BaseDialog {
    id: dialog
    objectName: "calculation-dialog"
    required property ApplicationSession session
    title: qsTr("Orbit solution")
    width: Math.min(550, parent.width - 40)
    closeAllowed: !session.busy
    autoAccept: false
    canAccept: !session.cancelling
    acceptText: session.cancelling ? qsTr("Cancelling...") : session.busy ? qsTr("Cancel") : qsTr("Done")
    onAcceptRequested: {
        if (session.busy) session.cancel()
        else close()
    }
    contentItem: ColumnLayout {
        spacing: 22
        CustomText { text: dialog.session.inputName; font.pixelSize: 14; font.weight: Font.DemiBold; Layout.fillWidth: true }
        CustomText {
            text: dialog.session.error || dialog.session.progress
            color: dialog.session.error.length ? AppTheme.red : AppTheme.muted
            wrapMode: Text.Wrap
            elide: Text.ElideNone
            Layout.fillWidth: true
            Layout.minimumHeight: 38
        }
        CustomProgressBar {
            Layout.fillWidth: true
            indeterminate: dialog.session.busy
            value: !dialog.session.busy && !dialog.session.error.length && dialog.session.progress !== qsTr("Cancelled") ? 1 : 0
        }
    }
}
