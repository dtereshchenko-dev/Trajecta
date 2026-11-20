import QtQuick
import QtQuick.Layouts
import TrajectaUI

RowLayout {
    id: heading
    property alias title: label.text
    property alias detail: detailLabel.text
    property int detailFontSize: 11
    spacing: 12
    implicitHeight: 30
    CustomText { id: label; font.pixelSize: 14; font.weight: Font.DemiBold; Layout.fillWidth: true }
    CustomText { id: detailLabel; secondary: true; font.pixelSize: heading.detailFontSize }
}
