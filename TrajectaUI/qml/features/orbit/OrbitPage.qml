import QtQuick
import QtQuick.Layouts
import TrajectaUI

Item {
    id: page
    required property SolutionViewModel solution
    property bool suspended: false
    signal importRequested()
    ColumnLayout {
        anchors.fill: parent
        spacing: 0
        visible: page.solution.available
        Rectangle {
            Layout.fillWidth: true
            implicitHeight: 52
            color: AppTheme.well
            RowLayout {
                anchors.fill: parent
                anchors.leftMargin: 24
                anchors.rightMargin: 24
                spacing: 16
                CustomText { text: qsTr("Heliocentric \u00b7 ECLM J2000"); secondary: true; font.pixelSize: 11 }
                CustomText { text: page.solution.scene.canAnimate ? qsTr("Osculating two-body") : qsTr("Static osculating geometry"); secondary: true; font.pixelSize: 11; Layout.fillWidth: true }
                CustomButton { objectName: "orbit-layers"; iconName: "sliders-horizontal"; accessibleLabel: qsTr("Orbit display"); onClicked: layers.open() }
            }
            Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line }
        }
        Loader {
            Layout.fillWidth: true
            Layout.fillHeight: true
            active: page.visible && page.solution.available
            sourceComponent: OrbitScene {
                solution: page.solution
                showPlanets: layers.planets
                showGrid: layers.grid
                showLabels: layers.labels
                showPerihelion: layers.perihelion
            }
        }
        OrbitPlayback { Layout.fillWidth: true; scene: page.solution.scene; suspended: page.suspended || !page.visible }
    }
    ColumnLayout {
        visible: !page.solution.available
        anchors.centerIn: parent
        spacing: 16
        CustomIcon { Layout.alignment: Qt.AlignHCenter; name: "orbit"; size: 34; color: AppTheme.muted }
        CustomText { text: qsTr("No orbital solution"); Layout.alignment: Qt.AlignHCenter }
        CustomButton { text: qsTr("Open observations"); iconName: "folder-open"; primary: true; Layout.alignment: Qt.AlignHCenter; onClicked: page.importRequested() }
    }
    OrbitLayerDialog { id: layers }
}
