pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Controls.Basic
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

BaseDialog {
    id: dialog
    objectName: "ephemeris-settings-dialog"
    required property EphemerisViewModel ephemerides
    title: qsTr("Ephemeris settings")
    width: Math.min(600, parent ? parent.width - 40 : 600)
    height: Math.min(690, parent ? parent.height - 40 : 690)
    contentItem: ScrollView {
        id: settingsScroll
        clip: true
        contentWidth: availableWidth
        ScrollBar.vertical: CustomScrollBar { policy: ScrollBar.AsNeeded }
        ColumnLayout {
            width: settingsScroll.availableWidth
            spacing: 14
            GridLayout {
                Layout.fillWidth: true
                columns: 2
                columnSpacing: 16
                rowSpacing: 6
                CustomText { text: qsTr("Time scale"); secondary: true; font.pixelSize: 10 }
                CustomText { text: qsTr("Time"); secondary: true; font.pixelSize: 10 }
                CustomComboBox {
                    objectName: "ephemeris-time-scale"
                    Layout.fillWidth: true
                    model: ["UTC", "TT", "TDB"]
                    currentIndex: dialog.ephemerides.scale
                    Accessible.name: qsTr("Time scale")
                    onActivated: dialog.ephemerides.scale = currentIndex
                }
                CustomTextField {
                    objectName: "ephemeris-start-time"
                    Layout.fillWidth: true
                    font.family: UI.monoFont
                    font.pixelSize: 12
                    text: dialog.ephemerides.startTime
                    color: dialog.ephemerides.timeValid ? AppTheme.ink : AppTheme.red
                    placeholderText: "00:00:00"
                    Accessible.name: qsTr("Time of the first epoch")
                    onTextEdited: dialog.ephemerides.startTime = text
                }
            }
            CustomCheckBox {
                objectName: "ephemeris-light-time"
                text: qsTr("Light time")
                checked: dialog.ephemerides.lightTime
                onToggled: dialog.ephemerides.lightTime = checked
            }
            CustomText { text: qsTr("Additional columns"); font.pixelSize: 12; font.weight: Font.DemiBold }
            GridLayout {
                Layout.fillWidth: true
                columns: settingsScroll.availableWidth < 480 ? 1 : 2
                columnSpacing: 16
                rowSpacing: 8
                Repeater {
                    model: dialog.ephemerides.table.columnOptions
                    delegate: CustomCheckBox {
                        required property var modelData
                        objectName: "ephemeris-column-" + modelData.key
                        Layout.fillWidth: true
                        text: modelData.label
                        checked: modelData.enabled
                        onToggled: dialog.ephemerides.table.setColumnEnabled(modelData.id, checked)
                    }
                }
            }
            CustomText {
                Layout.fillWidth: true
                text: qsTr("Dates use the selected time scale. Geocentric and topocentric positions include light-time correction when enabled; otherwise they are geometric. Heliocentric positions are always geometric. Stellar aberration is not applied to coordinates. The solution's force model is used, including its non-gravitational parameters.")
                secondary: true
                font.pixelSize: 11
                wrapMode: Text.Wrap
                elide: Text.ElideNone
            }
        }
    }
}
