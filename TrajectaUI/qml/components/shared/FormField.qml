pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

ColumnLayout {
    id: field
    property string label
    spacing: 6
    Layout.fillWidth: true
    Layout.preferredWidth: 1
    CustomText {
        text: field.label
        secondary: true
        font.pixelSize: 11
        Layout.fillWidth: true
        wrapMode: Text.Wrap
    }
}
