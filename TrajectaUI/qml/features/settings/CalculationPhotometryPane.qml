import QtQuick
import QtQuick.Layouts
import TrajectaUI

SettingsPane {
    id: pane
    required property CalculationSettingsModel settings
    title: qsTr("Photometry")
    SettingRow {
        label: qsTr("Asteroid magnitude slope G")
        CustomNumberField { value: pane.settings.photometry.asteroidG; Accessible.name: qsTr("Asteroid magnitude slope G"); onValueEdited: number => pane.settings.photometry.asteroidG = number }
    }
    SettingRow {
        label: qsTr("Comet total magnitude slope K1")
        CustomNumberField { value: pane.settings.photometry.cometK1; Accessible.name: qsTr("Comet total magnitude slope K1"); onValueEdited: number => pane.settings.photometry.cometK1 = number }
    }
    SettingRow {
        label: qsTr("Comet nuclear magnitude slope K2")
        CustomNumberField { value: pane.settings.photometry.cometK2; Accessible.name: qsTr("Comet nuclear magnitude slope K2"); onValueEdited: number => pane.settings.photometry.cometK2 = number }
    }
    SettingRow {
        label: qsTr("Comet nuclear phase coefficient")
        detail: qsTr("PHCOF · mag/degree")
        CustomNumberField { value: pane.settings.photometry.nuclearPhaseCoefficient; Accessible.name: qsTr("Comet nuclear phase coefficient"); onValueEdited: number => pane.settings.photometry.nuclearPhaseCoefficient = number }
    }
}
