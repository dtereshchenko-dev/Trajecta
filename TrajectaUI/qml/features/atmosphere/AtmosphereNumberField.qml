pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

ColumnLayout {
    id: field
    required property AtmosphereViewModel atmosphere
    required property string settingKey
    property string label
    property string placeholder
    Layout.fillWidth: true
    Layout.preferredWidth: 1
    spacing: 6
    CustomText {
        text: field.label
        visible: text.length > 0
        secondary: true
        font.pixelSize: 11
        Layout.fillWidth: true
        wrapMode: Text.Wrap
    }
    CustomTextField {
        objectName: "atmosphere-" + field.settingKey
        Layout.fillWidth: true
        Layout.minimumWidth: 65
        text: String(field.atmosphere.values[field.settingKey] ?? "")
        placeholderText: field.placeholder
        font.family: UI.monoFont
        font.pixelSize: 12
        Accessible.name: field.label || field.placeholder
        onTextEdited: field.atmosphere.setValue(field.settingKey, text)
    }
}
