pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

ColumnLayout {
    id: field
    required property AtmosphereViewModel atmosphere
    required property string settingKey
    required property string label
    Layout.fillWidth: true
    Layout.preferredWidth: 1
    spacing: 6
    CustomText { text: field.label; secondary: true; font.pixelSize: 11; Layout.fillWidth: true; wrapMode: Text.Wrap }
    RowLayout {
        Layout.fillWidth: true
        spacing: 8
        AtmosphereNumberField { atmosphere: field.atmosphere; settingKey: field.settingKey + "Min"; placeholder: field.label + qsTr(" minimum") }
        CustomText { text: qsTr("to"); secondary: true; font.pixelSize: 11 }
        AtmosphereNumberField { atmosphere: field.atmosphere; settingKey: field.settingKey + "Max"; placeholder: field.label + qsTr(" maximum") }
    }
}
