import QtQuick
import QtQuick.Controls.Basic
import TrajectaUI

ProgressBar {
    id: control
    implicitHeight: 4
    implicitWidth: 200
    padding: 0
    readonly property bool animating: indeterminate && visible && enabled
        && width > 0 && height > 0 && Window.window !== null && Window.window.visible
    background: Rectangle {
        radius: height / 2
        color: AppTheme.lineSoft
    }
    contentItem: Item {
        id: track
        property real phase: 0
        Rectangle {
            id: indicator
            objectName: "progress-indicator"
            x: control.indeterminate ? (parent.width - width) * track.phase : 0
            height: parent.height
            width: control.indeterminate ? Math.min(180, parent.width * 0.24) : parent.width * control.position
            radius: height / 2
            color: AppTheme.accent
            Behavior on width {
                enabled: !control.indeterminate && control.visible
                NumberAnimation { duration: 160; easing.type: Easing.OutCubic }
            }
        }
        SequentialAnimation on phase {
            running: control.animating
            loops: Animation.Infinite
            NumberAnimation { from: 0; to: 1; duration: 1500; easing.type: Easing.InOutSine }
            NumberAnimation { from: 1; to: 0; duration: 1500; easing.type: Easing.InOutSine }
        }
        transform: Scale {
            origin.x: track.width / 2
            xScale: control.mirrored ? -1 : 1
        }
    }
}
