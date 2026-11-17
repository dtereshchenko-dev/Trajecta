import QtQuick
import QtQuick.Layouts
import TrajectaUI

RowLayout {
    id: controls
    property string objectPrefix: "orbit"
    signal zoomIn()
    signal zoomOut()
    signal resetView()
    signal topView()

    spacing: 4
    CustomButton { objectName: controls.objectPrefix + "-zoom-in"; iconName: "plus"; accessibleLabel: qsTr("Zoom in"); onClicked: controls.zoomIn() }
    CustomButton { objectName: controls.objectPrefix + "-zoom-out"; iconName: "minus"; accessibleLabel: qsTr("Zoom out"); onClicked: controls.zoomOut() }
    CustomButton { objectName: controls.objectPrefix + "-reset-view"; iconName: "rotate-ccw"; accessibleLabel: qsTr("Reset view"); onClicked: controls.resetView() }
    CustomButton { objectName: controls.objectPrefix + "-top-view"; iconName: "view"; accessibleLabel: qsTr("Ecliptic view"); onClicked: controls.topView() }
}
