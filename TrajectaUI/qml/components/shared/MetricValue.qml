import QtQuick
import QtQuick.Layouts
import TrajectaUI

ColumnLayout {
    Layout.fillWidth: false
    property alias label: caption.text
    property alias value: amount.text
    property alias context: contextLabel.text
    property alias suffix: suffixLabel.text
    property int valueSize: 21
    spacing: 5
    CustomText { id: caption; secondary: true; font.pixelSize: 11; Layout.fillWidth: true }
    RowLayout {
        spacing: 6
        CustomText { id: amount; numeric: true; font.pixelSize: valueSize; font.weight: Font.Medium }
        CustomText { id: suffixLabel; secondary: true; font.pixelSize: 11; Layout.alignment: Qt.AlignBaseline }
    }
    CustomText { id: contextLabel; secondary: true; font.pixelSize: 10; Layout.fillWidth: true }
}
