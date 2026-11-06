import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI

Dialog {
    id: dialog
    property string acceptText: qsTr("Done")
    property string cancelText: ""
    property bool canAccept: true
    property bool acceptPrimary: false
    property bool closeAllowed: true
    property bool autoAccept: true
    readonly property var hostWindow: parent ? parent.Window.window : null
    signal acceptRequested()
    modal: true
    focus: true
    closePolicy: closeAllowed ? Popup.CloseOnEscape : Popup.NoAutoClose
    Overlay.modal: Rectangle { color: AppTheme.overlay }
    parent: Overlay.overlay
    onVisibleChanged: DialogState.update(dialog, hostWindow, visible && modal)
    onModalChanged: DialogState.update(dialog, hostWindow, visible && modal)
    onHostWindowChanged: DialogState.update(dialog, hostWindow, visible && modal)
    Component.onDestruction: DialogState.update(dialog, hostWindow, false)
    anchors.centerIn: parent
    width: Math.min(460, parent ? parent.width - 40 : 460)
    padding: 24
    spacing: 0
    background: Rectangle { color: AppTheme.surface; radius: 7; border.color: AppTheme.controlLine }
    header: Item {
        implicitHeight: 68
        RowLayout {
            anchors.fill: parent
            anchors.leftMargin: 22
            anchors.rightMargin: 22
            spacing: 16
            CustomText { text: dialog.title; font.pixelSize: 17; font.weight: Font.DemiBold; Layout.fillWidth: true }
            CustomButton { objectName: "dialog-close-button"; iconName: "x"; flat: true; enabled: dialog.closeAllowed; accessibleLabel: qsTr("Close"); onClicked: dialog.reject() }
        }
        Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line }
    }
    footer: Item {
        implicitHeight: 62
        Rectangle {
            anchors.fill: parent
            anchors.leftMargin: 1
            anchors.rightMargin: 1
            anchors.bottomMargin: 1
            bottomLeftRadius: 6
            bottomRightRadius: 6
            color: AppTheme.surfaceSecondary
        }
        Rectangle { x: 1; width: parent.width - 2; height: 1; color: AppTheme.line }
        RowLayout {
            anchors.fill: parent
            anchors.leftMargin: 22
            anchors.rightMargin: 22
            anchors.topMargin: 14
            anchors.bottomMargin: 14
            spacing: 10
            Item { Layout.fillWidth: true }
            CustomButton { text: dialog.cancelText; visible: text.length > 0; onClicked: dialog.reject() }
            CustomButton {
                objectName: "dialog-accept-button"
                text: dialog.acceptText
                enabled: dialog.canAccept
                primary: dialog.acceptPrimary
                onClicked: {
                    dialog.acceptRequested()
                    if (dialog.autoAccept) dialog.accept()
                }
            }
        }
    }
}
