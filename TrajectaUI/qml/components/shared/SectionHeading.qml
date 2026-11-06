import QtQuick
import QtQuick.Layouts
import TrajectaUI

RowLayout {
    property alias title: label.text
    property alias detail: detailLabel.text
    spacing: 12
    implicitHeight: 30
    CustomText { id: label; font.pixelSize: 14; font.weight: Font.DemiBold; Layout.fillWidth: true }
    CustomText { id: detailLabel; secondary: true; font.pixelSize: 11 }
}
