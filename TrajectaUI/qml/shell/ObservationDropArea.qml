import QtQuick
import QtQuick.Layouts
import TrajectaUI

DropArea {
    id: area
    objectName: "observation-drop-area"
    required property ApplicationSession session
    enabled: !session.busy
    keys: ["text/uri-list"]
    onEntered: function(drag) {
        drag.accepted = (drag.supportedActions & Qt.CopyAction) !== 0 && session.canDropInput(drag.urls)
    }
    onDropped: function(drop) {
        if ((drop.supportedActions & Qt.CopyAction) === 0 || !session.canDropInput(drop.urls)) {
            drop.accepted = false
            return
        }
        drop.accept(Qt.CopyAction)
        session.openInput(drop.urls[0])
    }
    Rectangle {
        anchors.fill: parent
        anchors.margins: 8
        visible: area.containsDrag
        color: Qt.alpha(AppTheme.base, 0.94)
        border.color: AppTheme.accent
        border.width: 2
        radius: 4
        ColumnLayout {
            anchors.centerIn: parent
            spacing: 14
            CustomIcon { name: "folder-open"; color: AppTheme.accent; size: 32; Layout.alignment: Qt.AlignHCenter }
            CustomText { text: qsTr("Open observations"); font.pixelSize: 17; font.weight: Font.Medium }
        }
    }
}
