import QtQuick
import QtQuick.Layouts
import TrajectaUI

RowLayout {
    id: root
    property string label: ""
    property string value: ""
    property string unit: ""
    property int valueSize: 11
    Layout.minimumHeight: 16
    Layout.fillWidth: false
    spacing: 4
    CustomText { objectName: "inline-label"; Layout.alignment: Qt.AlignBaseline; text: root.label; visible: text.length > 0; secondary: true; font.pixelSize: 11 }
    CustomText { objectName: "inline-value"; Layout.alignment: Qt.AlignBaseline; text: root.value; numeric: true; font.pixelSize: root.valueSize; font.weight: Font.Bold }
    CustomText { objectName: "inline-unit"; Layout.alignment: Qt.AlignBaseline; text: root.unit; visible: text.length > 0; secondary: true; font.pixelSize: 11 }
}
