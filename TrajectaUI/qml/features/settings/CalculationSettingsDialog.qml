pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

BaseDialog {
    id: dialog
    objectName: "calculation-settings-dialog"
    required property ApplicationSession session
    property int currentSection: 0
    readonly property CalculationSettingsModel settings: session.calculationSettings
    title: qsTr("Settings")
    width: Math.min(880, parent ? parent.width - 40 : 880)
    height: Math.min(820, parent ? parent.height - 40 : 820)
    padding: 1
    onAboutToShow: settings.reset()
    onClosed: settings.reset()
    contentItem: RowLayout {
        spacing: 0
        Rectangle {
            Layout.preferredWidth: 185
            Layout.fillHeight: true
            color: AppTheme.base
            ColumnLayout {
                anchors.fill: parent
                anchors.margins: 12
                anchors.topMargin: 17
                spacing: 4
                Repeater {
                    model: [
                        { name: qsTr("Appearance"), icon: "moon" },
                        { name: qsTr("Fit"), icon: "scan-line" },
                        { name: qsTr("Integrator"), icon: "activity" },
                        { name: qsTr("Perturbations"), icon: "orbit" },
                        { name: qsTr("Non-gravitational"), icon: "sigma" },
                        { name: qsTr("Photometry"), icon: "sun" },
                        { name: qsTr("About"), icon: "info" }
                    ]
                    NavigationItem {
                        required property int index
                        required property var modelData
                        objectName: "settings-section-" + index
                        Layout.fillWidth: true
                        text: modelData.name
                        iconName: modelData.icon
                        checked: dialog.currentSection === index
                        labelSize: 12
                        leftPadding: 10
                        rightPadding: 4
                        onClicked: dialog.currentSection = index
                    }
                }
                Item { Layout.fillHeight: true }
            }
            Rectangle { anchors.right: parent.right; width: 1; height: parent.height; color: AppTheme.line }
        }
        StackLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            currentIndex: dialog.currentSection
            SettingsAppearancePane { session: dialog.session }
            CalculationFitPane { settings: dialog.settings; enabled: !dialog.session.busy }
            CalculationIntegratorPane { settings: dialog.settings; sources: dialog.session.dataSources; enabled: !dialog.session.busy }
            CalculationPerturbationsPane { settings: dialog.settings; sources: dialog.session.dataSources; enabled: !dialog.session.busy }
            CalculationNonGravPane { settings: dialog.settings; enabled: !dialog.session.busy }
            CalculationPhotometryPane { settings: dialog.settings; enabled: !dialog.session.busy }
            CalculationAboutPane { settings: dialog.settings }
        }
    }
    footer: ColumnLayout {
        spacing: 0
        CustomText {
            objectName: "settings-validation-error"
            visible: !dialog.settings.valid
            text: dialog.settings.validationError
            color: AppTheme.red
            font.pixelSize: 11
            Layout.fillWidth: true
            Layout.margins: visible ? 14 : 0
            wrapMode: Text.Wrap
            elide: Text.ElideNone
        }
        Rectangle { Layout.fillWidth: true; implicitHeight: 1; color: AppTheme.line }
        RowLayout {
            Layout.fillWidth: true
            Layout.leftMargin: 22
            Layout.rightMargin: 22
            Layout.topMargin: 14
            Layout.bottomMargin: 14
            spacing: 10
            CustomText { text: qsTr("For the next calculation"); secondary: true; font.pixelSize: 10; Layout.fillWidth: true }
            CustomButton { objectName: "settings-cancel"; text: qsTr("Cancel"); onClicked: dialog.reject() }
            CustomButton { objectName: "settings-apply"; text: qsTr("Apply"); enabled: dialog.settings.valid && !dialog.session.busy; onClicked: if (dialog.session.applySettings()) dialog.accept() }
            CustomButton { objectName: "settings-recalculate"; text: qsTr("Apply and re-solve"); iconName: "play"; primary: true; enabled: dialog.settings.valid && dialog.session.canRun; onClicked: if (dialog.session.applySettings(true)) dialog.accept() }
        }
    }
}
