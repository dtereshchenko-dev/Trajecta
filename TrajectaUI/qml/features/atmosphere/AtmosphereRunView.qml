pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import QtQuick.Controls.Basic
import TrajectaUI

ColumnLayout {
    id: view
    objectName: "atmosphere-run-view"
    required property AtmosphereViewModel atmosphere
    spacing: 20
    RowLayout {
        Layout.fillWidth: true
        CustomText { text: qsTr("Atmospheric MC"); font.pixelSize: 14; font.weight: Font.DemiBold }
        CustomComboBox {
            objectName: "atmosphere-run-history"
            Layout.fillWidth: true
            model: view.atmosphere.runs
            currentIndex: view.atmosphere.selectedRun
            enabled: count > 0
            onActivated: view.atmosphere.selectedRun = currentIndex
            Accessible.name: qsTr("Atmospheric run")
        }
        CustomButton {
            objectName: "atmosphere-report"
            iconName: "file-text"
            enabled: view.atmosphere.reportAvailable
            accessibleLabel: qsTr("Open HTML report")
            onClicked: view.atmosphere.openResult(true)
            ToolTip.visible: hovered
            ToolTip.text: qsTr("Open HTML report")
        }
        CustomButton {
            iconName: "folder-open"
            enabled: view.atmosphere.resultPath.length > 0
            accessibleLabel: qsTr("Open numerical archive")
            onClicked: view.atmosphere.openResult()
            ToolTip.visible: hovered
            ToolTip.text: qsTr("Open numerical archive")
        }
    }
    CustomText {
        text: view.atmosphere.progress
        visible: view.atmosphere.busy
        Layout.fillWidth: true
        wrapMode: Text.Wrap
        elide: Text.ElideNone
        font.pixelSize: 12
    }
    CustomProgressBar {
        objectName: "atmosphere-progress"
        Layout.fillWidth: true
        visible: view.atmosphere.busy
        value: view.atmosphere.progressFraction
        indeterminate: value === 0
    }
    CustomText {
        text: view.atmosphere.summary
        Layout.fillWidth: true
        wrapMode: Text.Wrap
        elide: Text.ElideNone
        font.pixelSize: 12
    }
    CustomText {
        text: qsTr("Conditional on the selected atmospheric entry states")
        visible: view.atmosphere.selectedRun >= 0
        Layout.fillWidth: true
        secondary: true
        font.pixelSize: 11
    }
    CustomText {
        objectName: "atmosphere-run-error"
        text: view.atmosphere.runError
        color: AppTheme.amber
        visible: text.length > 0
        Layout.fillWidth: true
        wrapMode: Text.Wrap
        elide: Text.ElideNone
        font.pixelSize: 12
    }
    CustomText {
        text: view.atmosphere.resultPath
        Layout.fillWidth: true
        wrapMode: Text.WrapAnywhere
        elide: Text.ElideNone
        secondary: true
        font.pixelSize: 11
    }
}
