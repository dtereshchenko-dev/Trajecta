import QtQuick
import QtQuick.Layouts
import TrajectaUI

Item {
    id: row
    property string label: ""
    property string detail: ""
    default property alias controls: line.data
    implicitHeight: Math.max(33, caption.implicitHeight, line.implicitHeight) + 26
    Layout.fillWidth: true
    RowLayout {
        id: line
        anchors.fill: parent
        anchors.topMargin: 13
        anchors.bottomMargin: 13
        spacing: 22
        ColumnLayout {
            id: caption
            spacing: 4
            Layout.fillWidth: true
            CustomText { text: row.label; font.pixelSize: 12; Layout.fillWidth: true; wrapMode: Text.Wrap; elide: Text.ElideNone }
            CustomText { visible: text.length > 0; text: row.detail; font.pixelSize: 10; secondary: true; Layout.fillWidth: true; wrapMode: Text.Wrap; elide: Text.ElideNone }
        }
    }
    Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line; opacity: 0.5 }
}
