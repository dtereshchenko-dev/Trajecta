import QtQuick
import QtQuick.Layouts
import TrajectaUI

SettingsPane {
    id: pane
    required property CalculationSettingsModel settings
    title: qsTr("Fit and observations")
    SettingRow {
        label: qsTr("Astrometric catalogue debias")
        detail: "EFCC18"
        CustomSwitch { objectName: "settings-debias"; checked: pane.settings.weighting.debiasEnabled; Accessible.name: qsTr("Astrometric catalogue debias"); onToggled: pane.settings.weighting.debiasEnabled = checked }
    }
    SettingRow {
        label: qsTr("Observation weighting")
        CustomSwitch { objectName: "settings-weighting"; checked: pane.settings.weighting.useWeights; Accessible.name: qsTr("Observation weighting"); onToggled: pane.settings.weighting.useWeights = checked }
    }
    SettingRow {
        label: qsTr("Weighting profile")
        CustomComboBox { objectName: "settings-weight-profile"; model: ["vfcc17.rules", "Uniform 0.5 arcsec"]; currentIndex: pane.settings.weighting.weightProfile - 1; enabled: pane.settings.weighting.useWeights; Accessible.name: qsTr("Weighting profile"); onActivated: pane.settings.weighting.weightProfile = currentIndex + 1 }
    }
    SettingRow {
        label: qsTr("Prefer observer-supplied sigmas")
        detail: "ADES rmsRA / rmsDec"
        CustomSwitch { checked: pane.settings.weighting.observerSigmas; enabled: pane.settings.weighting.useWeights; Accessible.name: qsTr("Prefer observer-supplied sigmas"); onToggled: pane.settings.weighting.observerSigmas = checked }
    }
    SettingRow {
        label: qsTr("Timing uncertainty")
        CustomSwitch { checked: pane.settings.weighting.timingEnabled; enabled: pane.settings.weighting.useWeights; Accessible.name: qsTr("Timing uncertainty"); onToggled: pane.settings.weighting.timingEnabled = checked }
    }
    SettingRow {
        label: qsTr("Station-cluster reweighting")
        detail: "vfcc17"
        CustomText { text: pane.settings.weighting.useWeights && pane.settings.weighting.weightProfile === 1 ? qsTr("Automatic") : qsTr("Off"); secondary: true; font.pixelSize: 12 }
    }
    SettingRow {
        label: qsTr("Sigma-clip rejection")
        CustomSwitch { objectName: "settings-clipping"; checked: pane.settings.sigmaClip.enabled; Accessible.name: qsTr("Sigma-clip rejection"); onToggled: pane.settings.sigmaClip.enabled = checked }
    }
    SettingRow {
        label: qsTr("Rejection threshold, sigma")
        CustomNumberField { objectName: "settings-clip-threshold"; value: pane.settings.sigmaClip.threshold; minimumValue: 0; exclusiveMinimum: true; enabled: pane.settings.sigmaClip.enabled; Accessible.name: qsTr("Rejection threshold"); onValueEdited: number => pane.settings.sigmaClip.threshold = number }
    }
}
