import QtQuick
import QtCharts

LineSeries {
    property real radius: 0
    width: 1
    function rebuild() {
        clear()
        if (!Number.isFinite(radius) || radius <= 0) return
        for (var i = 0; i <= 128; ++i) {
            var angle = i * 2 * Math.PI / 128
            append(radius * Math.cos(angle), radius * Math.sin(angle))
        }
    }
    onRadiusChanged: rebuild()
    Component.onCompleted: rebuild()
}
