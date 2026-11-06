import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../../shared/UI.js" as UI

BaseDialog {
    id: dialog
    objectName: "ephemeris-settings-dialog"
    required property EphemerisViewModel ephemerides
    title: qsTr("Ephemeris settings")
    contentItem: ColumnLayout {
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
        CustomText {
            Layout.fillWidth: true
            text: qsTr("Dates on the page and in the table use the selected scale. Astrometric positions include light time only, without stellar aberration; heliocentric positions are geometric. The solution's force model is used, including its non-gravitational parameters.")
            secondary: true
            font.pixelSize: 11
            wrapMode: Text.Wrap
            elide: Text.ElideNone
        }
    }
}
