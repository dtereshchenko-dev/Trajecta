pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

SettingsPane {
    id: pane
    required property CalculationSettingsModel settings
    readonly property int defaultMask: settings.nonGravDefaultMask
    readonly property var masks: settings.nonGrav.model === 4 ? [0, 2] : [0, 3, 7, 15]
    title: qsTr("Non-gravitational forces")
    SettingRow {
        label: qsTr("Include in fit")
        CustomSwitch { objectName: "settings-nongrav"; checked: pane.settings.nonGrav.enabled; Accessible.name: qsTr("Include non-gravitational forces in fit"); onToggled: pane.settings.nonGrav.enabled = checked }
    }
    ColumnLayout {
        Layout.fillWidth: true
        spacing: 0
        enabled: pane.settings.nonGrav.enabled
        SettingRow {
            label: qsTr("Model g(r)")
            CustomComboBox {
                objectName: "settings-nongrav-model"
                model: ["Marsden H2O", "Marsden CO", "Inverse square", "Yarkovsky"]
                currentIndex: pane.settings.nonGrav.model - 1
                Accessible.name: qsTr("Non-gravitational model")
                onActivated: {
                    pane.settings.nonGrav.model = currentIndex + 1
                    pane.settings.nonGrav.solveMask = 0
                }
            }
        }
        SettingRow {
            label: qsTr("Estimated coefficients")
            CustomComboBox {
                objectName: "settings-nongrav-coefficients"
                model: pane.settings.nonGrav.model === 4 ? [qsTr("Model default"), "A2"] : [qsTr("Model default"), "A1 + A2", "A1 + A2 + A3", "A1 + A2 + A3 + DT"]
                currentIndex: pane.masks.indexOf(pane.settings.nonGrav.solveMask)
                Accessible.name: qsTr("Estimated coefficients")
                onActivated: pane.settings.nonGrav.solveMask = pane.masks[currentIndex]
            }
        }
        GridLayout {
            Layout.fillWidth: true
            Layout.topMargin: 12
            Layout.bottomMargin: 12
            columns: 2
            columnSpacing: 20
            Repeater {
                model: ["A1", "A2", "A3", "DT"]
                CustomCheckBox {
                    required property int index
                    required property string modelData
                    Layout.fillWidth: true
                    text: qsTr("Estimate %1").arg(modelData)
                    enabled: false
                    opacity: enabled ? 1 : 0.45
                    checked: ((pane.settings.nonGrav.solveMask || pane.defaultMask) & (1 << index)) !== 0
                }
            }
        }
    }
}
