import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Dialogs
import QtQuick.Layouts
import TrajectaUI

BaseDialog {
    id: dialog
    objectName: "import-observations-dialog"
    required property ApplicationSession session
    title: qsTr("Open observations")
    width: Math.min(550, parent.width - 40)
    height: Math.min(implicitHeight, parent.height - 40)
    acceptText: qsTr("Cancel")
    contentItem: ScrollView {
        id: scroll
        contentWidth: availableWidth
        implicitHeight: body.implicitHeight
        clip: true
        ScrollBar.vertical: CustomScrollBar {
            parent: scroll
            x: scroll.width - width
            y: scroll.topPadding
            height: scroll.availableHeight
        }
        ColumnLayout {
            id: body
            width: parent.width
            spacing: 22
            Rectangle {
                Layout.fillWidth: true
                implicitHeight: 210
                color: AppTheme.surfaceSecondary
                border.color: AppTheme.controlLine
                radius: 4
                ColumnLayout {
                    anchors.centerIn: parent
                    width: parent.width - 40
                    spacing: 14
                    CustomIcon { name: "folder-open"; size: 32; color: AppTheme.muted; Layout.alignment: Qt.AlignHCenter }
                    CustomText { text: qsTr("Observations file"); font.pixelSize: 14; font.weight: Font.DemiBold; Layout.alignment: Qt.AlignHCenter }
                    CustomText { text: "ADES XML / PSV · MPC80 / MPEC"; secondary: true; font.pixelSize: 12; Layout.alignment: Qt.AlignHCenter }
                    CustomButton { objectName: "choose-observations-button"; text: qsTr("Choose file"); iconName: "folder-open"; primary: true; Layout.alignment: Qt.AlignHCenter; onClicked: picker.open() }
                }
            }
            RecentObservationFiles {
                Layout.fillWidth: true
                recentFiles: dialog.session.recentFiles
                onFileRequested: function(index) { dialog.session.openRecent(index) }
            }
            CustomText { visible: text.length > 0; text: dialog.session.error; color: AppTheme.red; wrapMode: Text.Wrap; Layout.fillWidth: true }
        }
    }
    ObservationDropArea {
        objectName: "import-drop-area"
        parent: dialog.contentItem
        anchors.fill: parent
        session: dialog.session
        enabled: dialog.visible && !picker.visible && !session.busy
    }
    FileDialog {
        id: picker
        objectName: "observations-file-picker"
        title: qsTr("Choose observations file")
        nameFilters: [qsTr("Observations (*.xml *.psv *.txt *.mpc *.obs *.rwo)"), qsTr("All files (*)")]
        onAccepted: dialog.session.openInput(selectedFile)
    }
}
