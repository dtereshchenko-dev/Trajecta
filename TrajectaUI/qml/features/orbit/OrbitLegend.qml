pragma ComponentBehavior: Bound
import QtQuick
import QtQuick.Layouts
import TrajectaUI

Flow {
    id: legend
    objectName: "orbit-legend"
    required property var entries
    required property string selectedTarget
    signal targetSelected(string key)
    spacing: 4
    Repeater {
        model: legend.entries
        CustomButton {
            id: entry
            required property var modelData
            objectName: "orbit-follow-" + modelData.key
            text: modelData.label
            accessibleLabel: qsTr("Follow %1").arg(text)
            flat: true
            checked: legend.selectedTarget === modelData.key
            implicitHeight: 28
            implicitWidth: Math.min(implicitContentWidth + 16, legend.width)
            leftPadding: 8
            rightPadding: 8
            topPadding: 4
            bottomPadding: 4
            contentItem: RowLayout {
                spacing: 6
                Rectangle { implicitWidth: 6; implicitHeight: 6; radius: 3; color: entry.modelData.color }
                CustomText {
                    Layout.fillWidth: true
                    text: entry.text
                    color: entry.checked || entry.hovered ? AppTheme.ink : AppTheme.muted
                    font.pixelSize: 11
                    elide: Text.ElideRight
                }
            }
            onClicked: legend.targetSelected(modelData.key)
        }
    }
}
