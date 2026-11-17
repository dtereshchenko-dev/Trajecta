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
                CustomButton { text: qsTr("Heliocentric"); checkable: true; checked: !page.solution.scene.relativeMode; onClicked: page.solution.scene.relativeMode = false }
                CustomButton { objectName: "orbit-relative-mode"; text: qsTr("Relative motion"); checkable: true; checked: page.solution.scene.relativeMode; onClicked: page.solution.scene.relativeMode = true }
                CustomText { text: page.solution.scene.relativeMode ? (page.solution.scene.relative.rotatingFrame ? qsTr("Sun–planet rotating") : qsTr("Planet-centered inertial")) : page.solution.scene.canAnimate ? qsTr("Osculating two-body · ECLM J2000") : qsTr("Static osculating geometry"); secondary: true; font.pixelSize: 11; Layout.fillWidth: true }
                CustomButton { objectName: "orbit-layers"; visible: !page.solution.scene.relativeMode; iconName: "sliders-horizontal"; accessibleLabel: qsTr("Orbit display"); onClicked: layers.open() }
            }
            Rectangle { anchors.bottom: parent.bottom; width: parent.width; height: 1; color: AppTheme.line }
        }
        RowLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: 0
            Loader {
                id: sceneLoader
                Layout.fillWidth: true
                Layout.fillHeight: true
                active: page.visible && page.solution.available
                sourceComponent: page.solution.scene.relativeMode ? relativeComponent : heliocentricComponent
            }
            RelativeMotionPanel {
                visible: page.solution.scene.relativeMode
                Layout.fillHeight: true
                model: page.solution.scene.relative
                onCaptureSelected: index => {
                    if (page.solution.scene.seekCaptureStart(index)) {
                        playback.playing = false
                        if (sceneLoader.item) sceneLoader.item.focusOnObject()
                    }
                }
            }
        }
        OrbitPlayback { id: playback; Layout.fillWidth: true; scene: page.solution.scene; suspended: page.suspended || !page.visible }
    }
    Component { id: heliocentricComponent; OrbitScene { solution: page.solution; showPlanets: layers.planets; showGrid: layers.grid; showLabels: layers.labels; showPerihelion: layers.perihelion } }
    Component { id: relativeComponent; RelativeOrbitScene { solution: page.solution } }
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
