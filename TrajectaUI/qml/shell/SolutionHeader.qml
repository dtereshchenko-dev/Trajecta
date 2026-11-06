import QtQuick
import QtQuick.Layouts
import TrajectaUI
import "../shared/UI.js" as UI

Item {
    id: header
    required property SolutionViewModel solution
    property real viewportWidth: 1440
    implicitHeight: viewportWidth >= 1700 ? 87 : viewportWidth <= 1200 ? 73 : 78
    RowLayout {
        anchors.fill: parent
        anchors.leftMargin: header.viewportWidth <= 1200 ? 20 : 30
        anchors.rightMargin: anchors.leftMargin
        spacing: 13
        CustomText { text: header.solution.available ? header.solution.objectTitle : qsTr("Solution"); font.pixelSize: UI.headingSize; font.weight: Font.DemiBold; Layout.maximumWidth: header.width * 0.55 }
        Repeater {
            model: header.solution.available ? header.solution.orbit.orbitalClasses : []
            Rectangle {
                required property string modelData
                objectName: "solution-orbital-class"
                implicitWidth: classification.implicitWidth + 14
                implicitHeight: 20
                color: "transparent"
                radius: 3
                border.color: AppTheme.controlLine
                CustomText { id: classification; anchors.centerIn: parent; text: parent.modelData; secondary: true; font.pixelSize: 11 }
            }
        }
        CustomText { text: header.solution.designation !== header.solution.objectTitle ? header.solution.designation : ""; secondary: true; numeric: true; font.pixelSize: 11 }
        Item { Layout.fillWidth: true }
        Rectangle { width: 6; height: 6; radius: 3; color: header.solution.converged ? AppTheme.accent : header.solution.available ? AppTheme.amber : AppTheme.muted }
        CustomText { text: header.solution.status; font.pixelSize: 11 }
    }
}
