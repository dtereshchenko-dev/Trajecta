import QtQuick
import QtQuick.Controls.Basic
import TrajectaUI
import "../shared/UI.js" as UI

ToolButton {
    id: symbol
    property string name: ""
    property color color: AppTheme.muted
    property int size: 17
    implicitWidth: size
    implicitHeight: size
    padding: 0
    focusPolicy: Qt.NoFocus
    hoverEnabled: false
    enabled: false
    icon.source: name.length ? UI.getSVG(name) : ""
    icon.color: color
    icon.width: size
    icon.height: size
    palette.mid: color
    palette.buttonText: color
    palette.disabled.buttonText: color
    background: Item {}
}
