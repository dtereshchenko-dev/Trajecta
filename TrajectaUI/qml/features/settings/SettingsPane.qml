import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI

Flickable {
    id: pane
    property string title: ""
    default property alias fields: column.data
    contentWidth: width
    contentHeight: column.implicitHeight + 52
    clip: true
    boundsBehavior: Flickable.StopAtBounds
    ScrollBar.vertical: CustomScrollBar {}
    ColumnLayout {
        id: column
        x: 28
        y: 26
        width: pane.width - 56
        spacing: 0
        CustomText { text: pane.title; font.pixelSize: 17; font.weight: Font.DemiBold; Layout.fillWidth: true; Layout.bottomMargin: 18 }
    }
}
