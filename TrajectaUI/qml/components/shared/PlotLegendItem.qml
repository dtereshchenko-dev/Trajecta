import QtQuick
import QtQuick.Layouts
import TrajectaUI

RowLayout {
    property alias label: caption.text
    property alias color: swatch.color
    spacing: 6
    Rectangle { id: swatch; implicitWidth: 7; implicitHeight: 7; radius: 3.5 }
    CustomText { id: caption; font.pixelSize: 11; secondary: true }
}
