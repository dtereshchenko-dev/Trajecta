import QtQuick
import QtQuick.Layouts
import TrajectaUI

Rectangle {
    id: bar
    required property OrbitSceneModel scene
    property bool suspended: false
    property bool playing: false
    readonly property real daysPerSecond: [10, 50, 100][speed.currentIndex]
    implicitHeight: 58
    color: AppTheme.surface
    onSuspendedChanged: if (suspended) playing = false
    onVisibleChanged: if (!visible) playing = false
    Connections { target: bar.scene; function onChanged() { bar.playing = false } }
    Timer {
        property double lastTick: 0
        interval: 33
        repeat: true
        running: bar.playing && bar.visible && !bar.suspended
        onRunningChanged: lastTick = Date.now()
        onTriggered: {
            const now = Date.now()
            const elapsed = Math.min(0.25, (now - lastTick) / 1000)
            lastTick = now
            const previous = bar.scene.position
            bar.scene.position = Math.min(1, bar.scene.position + elapsed * bar.daysPerSecond / bar.scene.spanDays)
            if (bar.scene.position >= 1 || bar.scene.position === previous) bar.playing = false
        }
    }
    Rectangle { width: parent.width; height: 1; color: AppTheme.line }
    RowLayout {
        anchors.fill: parent
        anchors.leftMargin: 24
        anchors.rightMargin: 24
        spacing: 14
        CustomButton {
            objectName: "orbit-play"
            iconName: bar.playing ? "pause" : "play"
            accessibleLabel: bar.playing ? qsTr("Pause") : qsTr("Play")
            enabled: bar.scene.canAnimate && !bar.suspended
            onClicked: {
                if (bar.scene.position >= 1) bar.scene.position = 0
                bar.playing = !bar.playing
            }
        }
        CustomButton { objectName: "orbit-reset-time"; iconName: "rotate-ccw"; accessibleLabel: qsTr("Solution epoch"); enabled: bar.scene.canAnimate; onClicked: { bar.playing = false; bar.scene.position = 0 } }
        InlineValue { objectName: "orbit-utc"; label: "UTC"; value: bar.scene.epochUtc || qsTr("Unavailable"); Layout.preferredWidth: 235 }
        CustomSlider {
            objectName: "orbit-time"
            Layout.fillWidth: true
            from: 0
            to: 1
            value: bar.scene.position
            enabled: bar.scene.canAnimate && !bar.suspended
            Accessible.name: qsTr("Orbit epoch")
            onMoved: { bar.playing = false; bar.scene.position = value }
        }
        CustomComboBox {
            id: speed
            objectName: "orbit-speed"
            Layout.preferredWidth: 106
            model: ["10 d/s", "50 d/s", "100 d/s"]
            enabled: bar.scene.canAnimate
            Accessible.name: qsTr("Playback speed")
        }
    }
}
