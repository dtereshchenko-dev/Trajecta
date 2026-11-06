import QtQuick
import TrajectaUI

SettingsPane {
    id: pane
    required property CalculationSettingsModel settings
    required property DataSourcesModel sources
    title: qsTr("Integrator")
    SettingRow {
        label: qsTr("N-body propagation")
        detail: pane.sources.planetaryAvailable ? "DE440" : qsTr("DE440 unavailable")
        CustomSwitch { objectName: "settings-nbody"; checked: pane.settings.perturbations.enabled; enabled: pane.sources.planetaryAvailable || checked; Accessible.name: qsTr("N-body propagation"); onToggled: pane.settings.perturbations.enabled = checked }
    }
    SettingRow { label: qsTr("Method"); CustomText { text: pane.settings.perturbations.enabled || pane.settings.nonGrav.enabled ? "RA15 / Gauss-Radau" : "Kepler / two-body"; numeric: true; font.pixelSize: 12 } }
    SettingRow {
        label: qsTr("Initial step, days")
        CustomNumberField { objectName: "settings-initial-step"; value: pane.settings.perturbations.integratorInitialStepDays; minimumValue: 1e-6; enabled: pane.settings.perturbations.enabled || pane.settings.nonGrav.enabled; Accessible.name: qsTr("Initial step in days"); onValueEdited: number => pane.settings.perturbations.integratorInitialStepDays = number }
    }
    SettingRow { label: qsTr("Accuracy"); CustomText { text: qsTr("Core-controlled"); secondary: true; font.pixelSize: 12 } }
    SettingRow { label: qsTr("Light-time correction"); CustomText { text: qsTr("Enabled"); secondary: true; font.pixelSize: 12 } }
}
