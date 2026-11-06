pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

SettingsPane {
    id: pane
    required property CalculationSettingsModel settings
    required property DataSourcesModel sources
    title: qsTr("Perturbing bodies")
    GridLayout {
        Layout.fillWidth: true
        columns: 2
        rowSpacing: 8
        columnSpacing: 20
        enabled: pane.settings.perturbations.enabled
        opacity: enabled ? 1 : 0.45
        Repeater {
            model: BodyPalette.planets
            CustomCheckBox {
                required property var modelData
                Layout.fillWidth: true
                text: modelData.full
                checked: (pane.settings.perturbations.planetMask & (1 << modelData.bit)) !== 0
                onToggled: pane.settings.perturbations.planetMask = checked ? pane.settings.perturbations.planetMask | (1 << modelData.bit) : pane.settings.perturbations.planetMask & ~(1 << modelData.bit)
            }
        }
    }
    SettingRow {
        Layout.topMargin: 12
        label: qsTr("Small bodies SB441")
        detail: pane.sources.asteroidsAvailable ? qsTr("16 massive asteroids") : qsTr("SB441 unavailable")
        CustomSwitch { objectName: "settings-sb441"; checked: pane.settings.perturbations.useSb441; enabled: pane.settings.perturbations.enabled && (pane.sources.asteroidsAvailable || checked); Accessible.name: qsTr("Small bodies SB441"); onToggled: pane.settings.perturbations.useSb441 = checked }
    }
    GridLayout {
        Layout.fillWidth: true
        Layout.topMargin: 16
        columns: 2
        columnSpacing: 20
        enabled: pane.settings.perturbations.enabled && !pane.settings.perturbations.useSb441
        opacity: enabled ? 1 : 0.45
        Repeater {
            model: BodyPalette.asteroids
            CustomCheckBox {
                required property var modelData
                Layout.fillWidth: true
                text: modelData.full
                checked: pane.settings.perturbations.useSb441 || (pane.settings.perturbations.planetMask & (1 << modelData.bit)) !== 0
                onToggled: pane.settings.perturbations.planetMask = checked ? pane.settings.perturbations.planetMask | (1 << modelData.bit) : pane.settings.perturbations.planetMask & ~(1 << modelData.bit)
            }
        }
    }
}
