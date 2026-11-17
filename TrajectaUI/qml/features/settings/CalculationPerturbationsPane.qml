pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

SettingsPane {
    id: pane
    required property CalculationSettingsModel settings
    required property DataSourcesModel sources
    title: qsTr("Forces")
    CustomText { text: qsTr("Perturbing bodies"); font.weight: Font.DemiBold; Layout.topMargin: 8; Layout.bottomMargin: 12 }
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
    SettingRow {
        Layout.topMargin: 20
        label: qsTr("Relativity")
        CustomComboBox {
            objectName: "settings-relativity"
            model: [qsTr("Full EIH"), qsTr("EIH: Sun outer sum"), qsTr("Off")]
            currentIndex: pane.settings.perturbations.relativityMode
            enabled: pane.settings.perturbations.enabled
            Accessible.name: qsTr("Relativity model")
            onActivated: pane.settings.perturbations.relativityMode = currentIndex
        }
    }
    CustomText { text: qsTr("Earth gravity"); font.weight: Font.DemiBold; Layout.topMargin: 20; Layout.bottomMargin: 8 }
    GridLayout {
        Layout.fillWidth: true
        columns: 3
        enabled: pane.settings.perturbations.enabled && (pane.settings.perturbations.planetMask & (1 << 2)) !== 0
        opacity: enabled ? 1 : 0.45
        CustomCheckBox { objectName: "settings-earth-j2"; text: "J2"; checked: pane.settings.perturbations.earthJ2; onToggled: pane.settings.perturbations.earthJ2 = checked }
        CustomCheckBox { objectName: "settings-earth-j3"; text: "J3"; checked: pane.settings.perturbations.earthJ3; onToggled: pane.settings.perturbations.earthJ3 = checked }
        CustomCheckBox { objectName: "settings-earth-j4"; text: "J4"; checked: pane.settings.perturbations.earthJ4; onToggled: pane.settings.perturbations.earthJ4 = checked }
    }
    SettingRow {
        label: qsTr("Earth pole")
        CustomComboBox {
            objectName: "settings-earth-pole"
            model: [qsTr("Time-varying"), qsTr("Fixed J2000")]
            currentIndex: pane.settings.perturbations.rotatingEarthPole ? 0 : 1
            enabled: pane.settings.perturbations.enabled && (pane.settings.perturbations.planetMask & (1 << 2)) !== 0
            Accessible.name: qsTr("Earth pole")
            onActivated: pane.settings.perturbations.rotatingEarthPole = currentIndex === 0
        }
    }
}
